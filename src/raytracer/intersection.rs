use std::collections::BTreeSet;

use approx::{AbsDiffEq, UlpsEq};
use im::Vector;

use crate::core::{matrix::Matrix, tuple::Tuple};

use super::{ray::Ray, shape::Shape};

#[derive(Debug, Clone)]
pub struct Intersection<'shape> {
    pub t: f64,
    pub object: &'shape Shape,
    // TODO: maybe another type for triangle intersection
    pub u: f64,
    pub v: f64,
    trail: Vector<Matrix<4>>,
}

impl<'shape> Intersection<'shape> {
    pub fn new(t: f64, object: &'shape Shape, trail: Vector<Matrix<4>>) -> Self {
        Self {
            t,
            object,
            trail,
            u: 0.0,
            v: 0.0,
        }
    }

    pub fn new_with_uv(
        t: f64,
        object: &'shape Shape,
        trail: Vector<Matrix<4>>,
        u: f64,
        v: f64,
    ) -> Self {
        Self {
            t,
            object,
            u,
            v,
            trail,
        }
    }

    pub fn prepare_computations(&self, ray: &Ray, xs: &Intersections<'shape>) -> Computations {
        // TODO: nicer way to do the containers that doesn't mix with what happens last
        // TODO: need type annotation?
        let mut containers: Vec<&'shape Shape> = Vec::new();

        let mut n1 = None;
        let mut n2 = None;

        for i in &xs.0 {
            if i == self {
                n1 = containers
                    .last()
                    .map(|shape| shape.material.refractive_index);
            }

            match containers.iter().position(|shape| *shape == i.object) {
                Some(shape_idx) => {
                    containers.remove(shape_idx);
                }
                None => containers.push(i.object),
            }

            if i == self {
                n2 = containers
                    .last()
                    .map(|shape| shape.material.refractive_index);
                break;
            }
        }

        let point = ray.position(self.t);
        let eye_v = -ray.direction;

        let normal_v = self.object.normal_at(&point, self, &self.trail);
        let inside = normal_v.dot(&eye_v) < 0.0;
        let normal_v = if inside { -normal_v } else { normal_v };

        Computations {
            t: self.t,
            shape: self.object,
            trail: self.trail.clone(),
            point,
            over_point: point + normal_v * Tuple::default_epsilon(),
            under_point: point - normal_v * Tuple::default_epsilon(),
            eye_v,
            normal_v,
            reflect_v: ray.direction.reflect(&normal_v),
            inside,
            n1: n1.unwrap_or(1.0),
            n2: n2.unwrap_or(1.0),
        }
    }
}

impl PartialEq for Intersection<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.t
            .ulps_eq(&other.t, Self::default_epsilon(), Self::default_max_ulps())
    }
}

impl Eq for Intersection<'_> {}

impl PartialOrd for Intersection<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.t.partial_cmp(&other.t)
    }
}

impl Ord for Intersection<'_> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.t.partial_cmp(&other.t).expect("t is not NaN")
    }
}

impl AbsDiffEq for Intersection<'_> {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-5
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.t.abs_diff_eq(&other.t, epsilon)
    }
}

impl UlpsEq for Intersection<'_> {
    fn default_max_ulps() -> u32 {
        4
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        self.t.ulps_eq(&other.t, epsilon, max_ulps)
    }
}

// TODO: need to keep duplicates?
#[derive(Debug, Clone)]
pub struct Intersections<'shape>(pub BTreeSet<Intersection<'shape>>);

impl<'shape> Intersections<'shape> {
    pub fn new<const N: usize>(xs: [Intersection<'shape>; N]) -> Self {
        Self(BTreeSet::from(xs))
    }

    pub fn hit(&self) -> Option<Intersection> {
        self.0.iter().find(|i| i.t >= 0.0).cloned()
    }
}

/// Encapsulates precomputed information for an intersection.
pub struct Computations<'shape> {
    pub t: f64,
    pub shape: &'shape Shape,
    pub trail: Vector<Matrix<4>>,
    pub point: Tuple,
    pub over_point: Tuple,
    pub under_point: Tuple,
    pub eye_v: Tuple,
    pub normal_v: Tuple,
    pub reflect_v: Tuple,
    pub inside: bool,
    pub n1: f64,
    pub n2: f64,
}

impl Computations<'_> {
    /// The Schlick approximation for the Fresnel reflectance.
    /// Computes the reflectance, which is a number between 0 and 1
    /// representing the fraction of light reflected.
    pub fn schlick(&self) -> f64 {
        let mut cos = self.eye_v.dot(&self.normal_v);

        if self.n1 > self.n2 {
            let n = self.n1 / self.n2;
            let sin2_t = n.powi(2) * (1.0 - cos.powi(2));

            if sin2_t > 1.0 {
                return 1.0;
            }

            let cos_t = (1.0 - sin2_t).sqrt();
            cos = cos_t;
        }

        let r0 = ((self.n1 - self.n2) / (self.n1 + self.n2)).powi(2);

        r0 + (1.0 - r0) * (1.0 - cos).powi(5)
    }
}

#[cfg(test)]
mod tests {
    use approx::{assert_abs_diff_eq, assert_relative_eq};

    use crate::raytracer::material::Material;

    use super::*;

    #[test]
    fn test_intersection_new() {
        let shape = Shape::new_sphere();
        let intersection = Intersection::new(3.5, &shape, Vector::new());

        assert_relative_eq!(intersection.t, 3.5);
        assert_eq!(intersection.object, &shape);
    }

    #[test]
    fn hit_all_positive() {
        let shape = Shape::new_sphere();
        let i1 = Intersection::new(1.0, &shape, Vector::new());
        let i2 = Intersection::new(2.0, &shape, Vector::new());

        let xs = Intersections::new([i1.clone(), i2]);
        let hit = xs.hit().expect("valid hit");

        assert_abs_diff_eq!(hit, i1);
    }

    #[test]
    fn hit_some_negative() {
        let shape = Shape::new_sphere();
        let i1 = Intersection::new(-1.0, &shape, Vector::new());
        let i2 = Intersection::new(1.0, &shape, Vector::new());

        let xs = Intersections::new([i1, i2.clone()]);
        let hit = xs.hit().expect("valid hit");

        assert_abs_diff_eq!(hit, i2);
    }

    #[test]
    fn hit_all_negative() {
        let shape = Shape::new_sphere();
        let i1 = Intersection::new(-2.0, &shape, Vector::new());
        let i2 = Intersection::new(-1.0, &shape, Vector::new());

        let xs = Intersections::new([i1, i2]);
        let hit = xs.hit();

        assert_eq!(hit, None);
    }

    #[test]
    fn hit_lowest_nonnegative() {
        let shape = Shape::new_sphere();
        let i1 = Intersection::new(5.0, &shape, Vector::new());
        let i2 = Intersection::new(7.0, &shape, Vector::new());
        let i3 = Intersection::new(-3.0, &shape, Vector::new());
        let i4 = Intersection::new(2.0, &shape, Vector::new());

        let xs = Intersections::new([i1, i2, i3, i4.clone()]);
        let hit = xs.hit().expect("valid hit");

        assert_abs_diff_eq!(hit, i4);
    }

    #[test]
    fn precomputing_state_of_intersection() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Shape::new_sphere();
        let i = Intersection::new(4.0, &shape, Vector::new());

        let comps = i.prepare_computations(&r, &Intersections::new([i.clone()]));

        assert_abs_diff_eq!(comps.t, i.t);
        assert_eq!(comps.shape, i.object);
        assert_abs_diff_eq!(comps.point, Tuple::point(0.0, 0.0, -1.0));
        assert_abs_diff_eq!(comps.eye_v, Tuple::vector(0.0, 0.0, -1.0));
        assert_abs_diff_eq!(comps.normal_v, Tuple::vector(0.0, 0.0, -1.0));
    }

    #[test]
    fn precomputing_reflection_vector() {
        let shape = Shape::new_plane();
        let r = Ray::new(
            Tuple::point(0.0, 1.0, -1.0),
            Tuple::vector(0.0, -(2.0_f64.sqrt()) / 2.0, 2.0_f64.sqrt() / 2.0),
        );
        let i = Intersection::new(2.0_f64.sqrt(), &shape, Vector::new());

        let comps = i.prepare_computations(&r, &Intersections::new([i.clone()]));

        assert_abs_diff_eq!(
            comps.reflect_v,
            Tuple::vector(0.0, 2.0_f64.sqrt() / 2.0, 2.0_f64.sqrt() / 2.0)
        );
    }

    #[test]
    fn hit_when_intersection_outside() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Shape::new_sphere();
        let i = Intersection::new(4.0, &shape, Vector::new());

        let comps = i.prepare_computations(&r, &Intersections::new([i.clone()]));

        assert!(!comps.inside);
    }

    #[test]
    fn hit_when_intersection_inside() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Shape::new_sphere();
        let i = Intersection::new(1.0, &shape, Vector::new());

        let comps = i.prepare_computations(&r, &Intersections::new([i.clone()]));

        assert_abs_diff_eq!(comps.point, Tuple::point(0.0, 0.0, 1.0));
        assert_abs_diff_eq!(comps.eye_v, Tuple::vector(0.0, 0.0, -1.0));
        assert!(comps.inside);
        assert_abs_diff_eq!(comps.normal_v, Tuple::vector(0.0, 0.0, -1.0));
    }

    #[test]
    fn hit_should_offset_the_point() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Shape::new_sphere().with_transform(Matrix::translation(0.0, 0.0, 1.0));

        let i = Intersection::new(5.0, &shape, Vector::new());
        let comps = i.prepare_computations(&r, &Intersections::new([i.clone()]));
        let epsilon = Tuple::default_epsilon();

        assert!(comps.over_point.z < epsilon / 2.0);
        assert!(comps.point.z > comps.over_point.z);
    }

    #[test]
    fn finding_n1_and_n2_at_various_intersections() {
        let a = Shape::new_glass_sphere()
            .with_transform(Matrix::scaling(2.0, 2.0, 2.0))
            .with_material(Material {
                refractive_index: 1.5,
                ..Material::default()
            });

        let b = Shape::new_glass_sphere()
            .with_transform(Matrix::translation(0.0, 0.0, -0.25))
            .with_material(Material {
                refractive_index: 2.0,
                ..Material::default()
            });

        let c = Shape::new_glass_sphere()
            .with_transform(Matrix::translation(0.0, 0.0, 0.25))
            .with_material(Material {
                refractive_index: 2.5,
                ..Material::default()
            });

        let r = Ray::new(Tuple::point(0.0, 0.0, -4.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = Intersections::new([
            Intersection::new(2.0, &a, Vector::new()),
            Intersection::new(2.75, &b, Vector::new()),
            Intersection::new(3.25, &c, Vector::new()),
            Intersection::new(4.75, &b, Vector::new()),
            Intersection::new(5.25, &c, Vector::new()),
            Intersection::new(6.0, &a, Vector::new()),
        ]);

        let scenarios = [
            (0, 1.0, 1.5),
            (1, 1.5, 2.0),
            (2, 2.0, 2.5),
            (3, 2.5, 2.5),
            (4, 2.5, 1.5),
            (5, 1.5, 1.0),
        ];

        for (index, n1, n2) in scenarios {
            let comps =
                xs.0.iter()
                    .nth(index)
                    .unwrap()
                    .prepare_computations(&r, &xs);
            assert_abs_diff_eq!(comps.n1, n1);
            assert_abs_diff_eq!(comps.n2, n2);
        }
    }

    #[test]
    fn under_point_offset_below_surface() {
        let shape = Shape::new_glass_sphere().with_transform(Matrix::translation(0.0, 0.0, -1.0));

        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let i = Intersection::new(5.0, &shape, Vector::new());

        let xs = Intersections::new([i.clone()]);
        let comps = i.prepare_computations(&r, &xs);

        assert!(comps.under_point.z > Tuple::default_epsilon() / 2.0);
        assert!(comps.point.z < comps.under_point.z);
    }

    #[test]
    fn schlick_approximation_for_total_internal_reflection() {
        let shape = Shape::new_glass_sphere();

        let r = Ray::new(
            Tuple::point(0.0, 0.0, 2.0_f64.sqrt() / 2.0),
            Tuple::vector(0.0, 1.0, 0.0),
        );
        let xs = Intersections::new([
            Intersection::new(-(2.0_f64.sqrt()) / 2.0, &shape, Vector::new()),
            Intersection::new(2.0_f64.sqrt() / 2.0, &shape, Vector::new()),
        ]);

        let comps = xs.0.iter().nth(1).unwrap().prepare_computations(&r, &xs);

        assert_abs_diff_eq!(comps.schlick(), 1.0);
    }

    #[test]
    fn schlick_approximation_perpendicular_viewing_angle() {
        let shape = Shape::new_glass_sphere();

        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));
        let xs = Intersections::new([
            Intersection::new(-1.0, &shape, Vector::new()),
            Intersection::new(1.0, &shape, Vector::new()),
        ]);

        let comps = xs.0.iter().nth(1).unwrap().prepare_computations(&r, &xs);

        assert_abs_diff_eq!(comps.schlick(), 0.04);
    }

    #[test]
    fn schlick_approximation_small_angle_n2_greater_than_n1() {
        let shape = Shape::new_glass_sphere();

        let r = Ray::new(Tuple::point(0.0, 0.99, -2.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = Intersections::new([Intersection::new(1.8589, &shape, Vector::new())]);

        let comps = xs.0.iter().next().unwrap().prepare_computations(&r, &xs);

        assert_abs_diff_eq!(comps.schlick(), 0.48873, epsilon = Tuple::default_epsilon());
    }

    #[test]
    fn prepare_normal_on_smooth_triangle() {
        let p1 = Tuple::point(0.0, 1.0, 0.0);
        let p2 = Tuple::point(-1.0, 0.0, 0.0);
        let p3 = Tuple::point(1.0, 0.0, 0.0);
        let n1 = Tuple::vector(0.0, 1.0, 0.0);
        let n2 = Tuple::vector(-1.0, 0.0, 0.0);
        let n3 = Tuple::vector(1.0, 0.0, 0.0);

        let tri = Shape::new_smooth_triangle(p1, p2, p3, n1, n2, n3);
        let i = Intersection::new_with_uv(1.0, &tri, Vector::new(), 0.45, 0.25);
        let r = Ray::new(Tuple::point(-0.2, 0.3, -2.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = Intersections::new([i.clone()]);

        let comps = i.prepare_computations(&r, &xs);
        assert_eq!(comps.normal_v, Tuple::vector(-0.5547, 0.83205, 0.0));
    }
}
