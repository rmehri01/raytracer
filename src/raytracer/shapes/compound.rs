use std::collections::BTreeSet;

use crate::{
    core::matrix::Transformation,
    raytracer::{
        bounds::{Bounded, Bounds},
        intersection::Intersections,
        ray::Ray,
    },
};

use super::{HasProperties, Intersect, Properties, Shape, Single};

#[derive(Debug, PartialEq, Clone)]
pub struct Compound {
    properties: Properties,
    kind: CompoundKind,
}

impl Compound {
    fn new(kind: CompoundKind) -> Self {
        Self {
            properties: Properties::default(),
            kind,
        }
    }

    pub fn new_group(children: Vec<Shape>) -> Self {
        Self::new(CompoundKind::Group(children))
    }

    pub fn new_csg(operation: Operation, left: Shape, right: Shape) -> Self {
        Self::new(CompoundKind::Csg {
            operation,
            left: Box::new(left),
            right: Box::new(right),
        })
    }

    pub fn includes(&self, other: &Single) -> bool {
        match &self.kind {
            CompoundKind::Group(children) => children.iter().any(|child| child.includes(other)),
            CompoundKind::Csg { left, right, .. } => left.includes(other) || right.includes(other),
        }
    }

    pub fn as_shape(self) -> Shape {
        Shape::Compound(self)
    }
}

impl HasProperties for Compound {
    fn properties(&self) -> &Properties {
        &self.properties
    }

    fn properties_mut(&mut self) -> &mut Properties {
        &mut self.properties
    }
}

impl Intersect for Compound {
    fn local_intersect(&self, ray: &Ray, trail: &im_rc::Vector<Transformation>) -> Intersections {
        let mut new_trail = trail.clone();
        new_trail.push_front(self.properties.transform);

        match &self.kind {
            CompoundKind::Group(children) => {
                if self.bounds().intersects(ray) {
                    let intersections = children.iter().fold(BTreeSet::new(), |mut acc, child| {
                        acc.append(&mut child.intersect(ray, &new_trail).0);
                        acc
                    });

                    Intersections(intersections)
                } else {
                    Intersections(BTreeSet::new())
                }
            }
            CompoundKind::Csg {
                operation,
                left,
                right,
            } => {
                let mut left_intersections = left.intersect(ray, &new_trail).0;
                left_intersections.append(&mut right.intersect(ray, &new_trail).0);

                operation.filter_intersections(left, Intersections(left_intersections))
            }
        }
    }
}

impl Bounded for Compound {
    fn bounds(&self) -> Bounds {
        match &self.kind {
            // TODO: recursively combine
            CompoundKind::Group(children) => {
                children
                    .iter()
                    .fold(Bounds::default(), |mut bounds, child| {
                        bounds.transform(&child.properties().transform);
                        bounds
                    })
            }
            CompoundKind::Csg { left, right, .. } => {
                let mut bounds = Bounds::default();

                bounds.transform(&left.properties().transform);
                bounds.transform(&right.properties().transform);

                bounds
            }
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
enum CompoundKind {
    /// A collection of shapes that are transformed as a unit.
    Group(Vec<Shape>),
    /// A constructive solid geometry shape that is composed of an operation
    /// and two operand shapes.
    Csg {
        operation: Operation,
        left: Box<Shape>,
        right: Box<Shape>,
    },
}

/// The possible CSG operations.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Operation {
    Union,
    Intersection,
    Difference,
}

impl Operation {
    fn filter_intersections<'shape>(
        self,
        left: &Shape,
        intersections: Intersections<'shape>,
    ) -> Intersections<'shape> {
        let mut in_l = false;
        let mut in_r = false;

        let intersections = intersections
            .0
            .into_iter()
            .filter(|i| {
                let l_hit = left.includes(i.shape);
                let allowed = self.intersection_allowed(l_hit, in_l, in_r);

                if l_hit {
                    in_l = !in_l;
                } else {
                    in_r = !in_r;
                }

                allowed
            })
            .collect();

        Intersections(intersections)
    }

    /// Determines which intersections should be preserved.
    ///
    /// * `l_hit`: True if the left shape was hit, and false if the right shape was hit.
    /// * `in_l`: True if the hit occurs inside the left shape.
    /// * `in_r`: True if the hit occurs inside the right shape.
    fn intersection_allowed(self, l_hit: bool, in_l: bool, in_r: bool) -> bool {
        match self {
            Self::Union => (l_hit && !in_r) || (!l_hit && !in_l),
            Self::Intersection => (l_hit && in_r) || (!l_hit && in_l),
            Self::Difference => (l_hit && !in_r) || (!l_hit && in_l),
        }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use crate::{
        core::{matrix::Matrix, point::Point, vector::Vector},
        raytracer::{intersection::Intersection, shapes::SetProperties},
    };

    use super::*;

    #[test]
    fn create_group() {
        let group = Compound::new_group(Vec::new());

        assert_eq!(group.properties.transform, Matrix::identity());
        if let CompoundKind::Group(children) = group.kind {
            assert!(children.is_empty());
        } else {
            panic!("expected a group");
        }
    }

    #[test]
    fn create_group_with_children() {
        let sphere = Single::new_sphere().as_shape();
        let group = Compound::new_group(vec![sphere.clone()]);

        assert_eq!(group.properties.transform, Matrix::identity());
        if let CompoundKind::Group(children) = group.kind {
            assert_eq!(children.len(), 1);
            assert_eq!(children[0], sphere);
        } else {
            panic!("expected a group");
        }
    }

    #[test]
    fn intersect_ray_with_empty_group() {
        let group = Compound::new_group(Vec::new());

        let r = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));
        let xs = group.intersect(&r, &im_rc::Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn intersect_ray_with_group_of_shapes() {
        let s1 = Single::new_sphere();
        let s2 = Single::new_sphere().with_transform(Matrix::translation(0.0, 0.0, -3.0));
        let s3 = Single::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));
        let group = Compound::new_group(vec![
            s1.clone().as_shape(),
            s2.clone().as_shape(),
            s3.as_shape(),
        ]);

        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let xs = group
            .intersect(&r, &im_rc::Vector::new())
            .0
            .iter()
            .map(|i| i.shape)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 4);
        assert_eq!(xs[0], &s2);
        assert_eq!(xs[1], &s2);
        assert_eq!(xs[2], &s1);
        assert_eq!(xs[3], &s1);
    }

    #[test]
    fn intersect_ray_with_transformed_group() {
        let s = Single::new_sphere()
            .with_transform(Matrix::translation(5.0, 0.0, 0.0))
            .as_shape();
        let group = Compound::new_group(vec![s]).with_transform(Matrix::scaling(2.0, 2.0, 2.0));

        let r = Ray::new(Point::new(10.0, 0.0, -10.0), Vector::new(0.0, 0.0, 1.0));
        let xs = group.intersect(&r, &im_rc::Vector::new());

        assert_eq!(xs.0.len(), 2);
    }

    #[test]
    fn evaluate_rule_for_csg_operation() {
        let scenarios = [
            (Operation::Union, true, true, true, false),
            (Operation::Union, true, true, false, true),
            (Operation::Union, true, false, true, false),
            (Operation::Union, true, false, false, true),
            (Operation::Union, false, true, true, false),
            (Operation::Union, false, true, false, false),
            (Operation::Union, false, false, true, true),
            (Operation::Union, false, false, false, true),
            (Operation::Intersection, true, true, true, true),
            (Operation::Intersection, true, true, false, false),
            (Operation::Intersection, true, false, true, true),
            (Operation::Intersection, true, false, false, false),
            (Operation::Intersection, false, true, true, true),
            (Operation::Intersection, false, true, false, true),
            (Operation::Intersection, false, false, true, false),
            (Operation::Intersection, false, false, false, false),
            (Operation::Difference, true, true, true, false),
            (Operation::Difference, true, true, false, true),
            (Operation::Difference, true, false, true, false),
            (Operation::Difference, true, false, false, true),
            (Operation::Difference, false, true, true, true),
            (Operation::Difference, false, true, false, true),
            (Operation::Difference, false, false, true, false),
            (Operation::Difference, false, false, false, false),
        ];

        for (op, l_hit, in_l, in_r, result) in scenarios {
            assert_eq!(op.intersection_allowed(l_hit, in_l, in_r), result);
        }
    }

    #[test]
    fn filtering_intersections() {
        let scenarios = [
            (Operation::Union, 0, 3),
            (Operation::Intersection, 1, 2),
            (Operation::Difference, 0, 1),
        ];

        for (operation, x0, x1) in scenarios {
            let s1 = Single::new_sphere();
            let s2 = Single::new_cube();
            let xs = Intersections::new([
                Intersection::new(1.0, &s1, im_rc::Vector::new()),
                Intersection::new(2.0, &s2, im_rc::Vector::new()),
                Intersection::new(3.0, &s1, im_rc::Vector::new()),
                Intersection::new(4.0, &s2, im_rc::Vector::new()),
            ]);

            let result = operation
                .filter_intersections(&s1.clone().as_shape(), xs.clone())
                .0
                .iter()
                .map(|i| i.t)
                .collect::<Vec<_>>();

            assert_eq!(result.len(), 2);
            assert_abs_diff_eq!(result[0], xs.0.iter().nth(x0).unwrap().t);
            assert_abs_diff_eq!(result[1], xs.0.iter().nth(x1).unwrap().t);
        }
    }

    #[test]
    fn ray_misses_csg_object() {
        let c = Compound::new_csg(
            Operation::Union,
            Single::new_sphere().as_shape(),
            Single::new_cube().as_shape(),
        );
        let r = Ray::new(Point::new(0.0, 2.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let xs = c.intersect(&r, &im_rc::Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_hits_csg_object() {
        let s1 = Single::new_sphere();
        let s2 = Single::new_sphere().with_transform(Matrix::translation(0.0, 0.0, 0.5));

        let c = Compound::new_csg(
            Operation::Union,
            s1.clone().as_shape(),
            s2.clone().as_shape(),
        );
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let intersections = c.intersect(&r, &im_rc::Vector::new());
        let xs = intersections.0.iter().collect::<Vec<_>>();

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0].t, 4.0);
        assert_eq!(*xs[0].shape, s1);
        assert_abs_diff_eq!(xs[1].t, 6.5);
        assert_eq!(*xs[1].shape, s2);
    }
}
