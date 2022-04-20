use std::collections::BTreeSet;

use approx::AbsDiffEq;

use crate::core::{
    matrix::{Matrix, Transformation},
    point::Point,
    vector::Vector,
};

use super::{
    bounds::Bounds,
    intersection::{Intersection, Intersections},
    material::Material,
    ray::Ray,
};

#[derive(Debug, PartialEq, Clone)]
pub struct Shape {
    pub transform: Transformation,
    transform_inversed: Transformation,
    pub material: Material,
    pub kind: ShapeKind,
}

impl Shape {
    fn new(kind: ShapeKind) -> Self {
        Self {
            transform: Matrix::identity(),
            transform_inversed: Matrix::identity(),
            material: Material::default(),
            kind,
        }
    }

    pub fn new_sphere() -> Self {
        Self::new(ShapeKind::Single(SingleKind::Sphere))
    }

    pub fn new_glass_sphere() -> Self {
        Self::new_sphere().with_material(Material {
            transparency: 1.0,
            refractive_index: 1.5,
            ..Material::default()
        })
    }

    pub fn new_plane() -> Self {
        Self::new(ShapeKind::Single(SingleKind::Plane))
    }

    pub fn new_cube() -> Self {
        Self::new(ShapeKind::Single(SingleKind::Cube))
    }

    pub fn new_cylinder(conic: Conic) -> Self {
        Self::new(ShapeKind::Single(SingleKind::Cylinder(conic)))
    }

    pub fn new_cone(conic: Conic) -> Self {
        Self::new(ShapeKind::Single(SingleKind::Cone(conic)))
    }

    pub fn new_group(children: Vec<Self>) -> Self {
        Self::new(ShapeKind::Group(children))
    }

    pub fn new_triangle(p1: Point, p2: Point, p3: Point) -> Self {
        let triangular = Triangular::new(p1, p2, p3);
        let normal = triangular.e1.cross(&triangular.e2).normalize();
        Self::new(ShapeKind::Single(SingleKind::Triangle {
            triangular,
            normal,
        }))
    }

    pub fn new_smooth_triangle(
        p1: Point,
        p2: Point,
        p3: Point,
        n1: Vector,
        n2: Vector,
        n3: Vector,
    ) -> Self {
        let triangular = Triangular::new(p1, p2, p3);

        Self::new(ShapeKind::SmoothTriangle {
            triangular,
            n1,
            n2,
            n3,
        })
    }

    pub fn new_csg(operation: Operation, left: Self, right: Self) -> Self {
        Self::new(ShapeKind::CSG {
            operation,
            left: Box::new(left),
            right: Box::new(right),
        })
    }

    pub fn with_transform(mut self, transform: Transformation) -> Self {
        self.transform = transform;
        self.transform_inversed = transform.inverse();
        self
    }

    pub fn with_material(mut self, material: Material) -> Self {
        self.material = material;
        self
    }

    pub fn intersect(&self, ray: &Ray, mut trail: im::Vector<Transformation>) -> Intersections {
        let local_ray = ray.transform(&self.transform_inversed);

        match &self.kind {
            ShapeKind::Single(single) => {
                let ts = single.intersect(&local_ray);
                let intersections = ts
                    .iter()
                    .map(|t| Intersection::new(*t, self, trail.clone()))
                    .collect();

                Intersections(intersections)
            }
            ShapeKind::Group(children) => {
                trail.push_front(self.transform_inversed);
                if self.bounds().intersects(&local_ray) {
                    let intersections = children
                        .iter()
                        .flat_map(|child| child.intersect(&local_ray, trail.clone()).0)
                        .collect();

                    Intersections(intersections)
                } else {
                    Intersections(BTreeSet::new())
                }
            }
            ShapeKind::SmoothTriangle { triangular, .. } => {
                // TODO: kind of duplicated with single case
                let (ts, u, v) = triangular.intersect_uv(&local_ray);
                let intersections = ts
                    .iter()
                    .map(|t| Intersection::new_with_uv(*t, self, trail.clone(), u, v))
                    .collect();

                Intersections(intersections)
            }
            ShapeKind::CSG {
                operation,
                left,
                right,
            } => {
                // TODO: could be cleaner?
                let mut left_intersections = left.intersect(&local_ray, trail.clone()).0;
                let mut right_intersections = right.intersect(&local_ray, trail.clone()).0;
                left_intersections.append(&mut right_intersections);

                operation.filter_intersections(left, Intersections(left_intersections))
            }
        }
    }

    fn bounds(&self) -> Bounds {
        match &self.kind {
            ShapeKind::Single(single) => single.bounds(),
            ShapeKind::Group(children) => {
                children.iter().fold(Bounds::default(), |bounds, child| {
                    bounds.transform(&child.transform)
                })
            }
            ShapeKind::SmoothTriangle { triangular, .. } => triangular.bounds(),
            ShapeKind::CSG { left, right, .. } => Bounds::default()
                .transform(&left.transform)
                .transform(&right.transform),
        }
    }

    // TODO: hit only used for smooth triangles
    pub fn normal_at(
        &self,
        point: &Point,
        hit: &Intersection,
        trail: &im::Vector<Transformation>,
    ) -> Vector {
        match &self.kind {
            ShapeKind::Single(single) => {
                let local_point = self.world_to_object(point, trail);
                let local_normal = single.normal_at(&local_point);

                self.normal_to_world(&local_normal, trail)
            }
            // TODO: duplicated with Single
            &ShapeKind::SmoothTriangle { n1, n2, n3, .. } => self.normal_to_world(
                &(n2 * hit.u + n3 * hit.v + n1 * (1.0 - hit.u - hit.v)),
                trail,
            ),
            ShapeKind::Group(_) => panic!("groups should not have normal vectors"),
            ShapeKind::CSG { .. } => panic!("csg objects should not have normal vectors"),
        }
    }

    pub fn world_to_object(
        &self,
        world_point: &Point,
        trail: &im::Vector<Transformation>,
    ) -> Point {
        let trail_point = trail.iter().rev().fold(*world_point, |acc, &mat| mat * acc);

        self.transform_inversed * trail_point
    }

    fn normal_to_world(&self, normal: &Vector, trail: &im::Vector<Transformation>) -> Vector {
        let normal = self.transform_inversed.transpose() * *normal;
        let normal = normal.normalize();

        trail.iter().fold(normal, |acc, mat| {
            let normal = mat.transpose() * acc;
            normal.normalize()
        })
    }

    fn includes(&self, other: &Self) -> bool {
        match &self.kind {
            ShapeKind::Group(children) => children.iter().any(|child| child.includes(other)),
            ShapeKind::CSG { left, right, .. } => left.includes(other) || right.includes(other),
            _ => self == other,
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub enum ShapeKind {
    /// A single shape.
    Single(SingleKind),
    // TODO: might move group and smooth triangle up so that they can return different types
    /// A collection of shapes that are transformed as a unit.
    Group(Vec<Shape>),
    /// A triangle with vertices at p1, p2, and p3 and a normal at each of
    /// the vertices. Uses normal interpolation to calculate the normal at
    /// any point on the triangle.
    SmoothTriangle {
        triangular: Triangular,
        n1: Vector,
        n2: Vector,
        n3: Vector,
    },
    /// A constructive solid geometry shape that is composed of an operation
    /// and two operand shapes.
    // TODO: subtype of group?
    CSG {
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
                let l_hit = left.includes(i.object);
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

#[derive(Debug, PartialEq, Clone)]
pub enum SingleKind {
    /// A unit sphere with center at the origin.
    Sphere,
    /// A perfectly flat surface that extends infinitely in x and z.
    Plane,
    /// An axis-aligned bounding box that is centered at the origin and
    /// extends from -1 to 1 along each axis.
    Cube,
    /// A cylinder that is centered at the origin with radius 1 and extends
    /// from minimum to maximum exclusive along the y axis.
    Cylinder(Conic),
    /// A double-napped cone that is centered at the origin and extends
    /// from minimum to maximum exclusive along the y axis.
    Cone(Conic),
    /// A triangle with vertices at p1, p2, and p3, along with two edge
    /// vectors and a normal vector to optimize intersection calculations.
    Triangle {
        triangular: Triangular,
        normal: Vector,
    },
}

impl SingleKind {
    pub fn intersect(&self, ray: &Ray) -> Vec<f64> {
        match self {
            Self::Sphere => Self::sphere_intersect(ray),
            Self::Plane => Self::plane_intersect(ray),
            Self::Cube => Self::cube_intersect(ray),
            Self::Cylinder(conic) => Self::cylinder_intersect(ray, conic),
            Self::Cone(conic) => Self::cone_intersect(ray, conic),
            Self::Triangle { triangular, .. } => triangular.intersect(ray),
        }
    }

    fn sphere_intersect(ray: &Ray) -> Vec<f64> {
        let sphere_to_ray = ray.origin - Point::new(0.0, 0.0, 0.0);

        let a = ray.direction.dot(&ray.direction);
        let b = 2.0 * sphere_to_ray.dot(&ray.direction);
        let c = sphere_to_ray.dot(&sphere_to_ray) - 1.0;

        let discriminant = b.powi(2) - 4.0 * a * c;

        if discriminant < 0.0 {
            Vec::new()
        } else {
            let t1 = (-b - discriminant.sqrt()) / (2.0 * a);
            let t2 = (-b + discriminant.sqrt()) / (2.0 * a);

            vec![t1, t2]
        }
    }

    fn plane_intersect(ray: &Ray) -> Vec<f64> {
        if ray.direction.y.abs() < Vector::default_epsilon() {
            Vec::new()
        } else {
            vec![-ray.origin.y / ray.direction.y]
        }
    }

    fn cube_intersect(ray: &Ray) -> Vec<f64> {
        let (x_t_min, x_t_max) = Self::check_axis(ray.origin.x, ray.direction.x);
        let (y_t_min, y_t_max) = Self::check_axis(ray.origin.y, ray.direction.y);
        let (z_t_min, z_t_max) = Self::check_axis(ray.origin.z, ray.direction.z);

        let t_min = x_t_min.max(y_t_min).max(z_t_min);
        let t_max = x_t_max.min(y_t_max).min(z_t_max);

        if t_min > t_max {
            Vec::new()
        } else {
            vec![t_min, t_max]
        }
    }

    fn check_axis(origin: f64, direction: f64) -> (f64, f64) {
        let tmin_numerator = -1.0 - origin;
        let tmax_numerator = 1.0 - origin;

        let tmin = tmin_numerator / direction;
        let tmax = tmax_numerator / direction;

        if tmin > tmax {
            (tmax, tmin)
        } else {
            (tmin, tmax)
        }
    }

    fn cylinder_intersect(ray: &Ray, conic: &Conic) -> Vec<f64> {
        let a = ray.direction.x.powi(2) + ray.direction.z.powi(2);
        let b = 2.0 * ray.origin.x * ray.direction.x + 2.0 * ray.origin.z * ray.direction.z;
        let c = ray.origin.x.powi(2) + ray.origin.z.powi(2) - 1.0;

        conic.intersect(ray, a, b, c, |_| 1.0)
    }

    fn cone_intersect(ray: &Ray, conic: &Conic) -> Vec<f64> {
        let a = ray.direction.x.powi(2) - ray.direction.y.powi(2) + ray.direction.z.powi(2);
        let b = 2.0 * ray.origin.x * ray.direction.x - 2.0 * ray.origin.y * ray.direction.y
            + 2.0 * ray.origin.z * ray.direction.z;
        let c = ray.origin.x.powi(2) - ray.origin.y.powi(2) + ray.origin.z.powi(2);

        conic.intersect(ray, a, b, c, f64::abs)
    }

    pub fn normal_at(&self, point: &Point) -> Vector {
        match self {
            Self::Sphere => Self::sphere_normal_at(point),
            Self::Plane => Self::plane_normal_at(),
            Self::Cube => Self::cube_normal_at(point),
            Self::Cylinder(conic) => Self::cylinder_normal_at(point, conic),
            Self::Cone(conic) => Self::cone_normal_at(point, conic),
            Self::Triangle { normal, .. } => *normal,
        }
    }

    fn sphere_normal_at(object_point: &Point) -> Vector {
        *object_point - Point::new(0.0, 0.0, 0.0)
    }

    fn plane_normal_at() -> Vector {
        Vector::new(0.0, 1.0, 0.0)
    }

    fn cube_normal_at(object_point: &Point) -> Vector {
        let max_c = object_point
            .x
            .abs()
            .max(object_point.y.abs())
            .max(object_point.z.abs());

        if max_c.abs_diff_eq(&object_point.x.abs(), Point::default_epsilon()) {
            Vector::new(object_point.x, 0.0, 0.0)
        } else if max_c.abs_diff_eq(&object_point.y.abs(), Point::default_epsilon()) {
            Vector::new(0.0, object_point.y, 0.0)
        } else {
            Vector::new(0.0, 0.0, object_point.z)
        }
    }

    fn cylinder_normal_at(object_point: &Point, conic: &Conic) -> Vector {
        conic.normal_at(object_point, 0.0)
    }

    fn cone_normal_at(object_point: &Point, conic: &Conic) -> Vector {
        let y = object_point.x.hypot(object_point.z);
        let y = if object_point.y > 0.0 { -y } else { y };

        conic.normal_at(object_point, y)
    }

    fn bounds(&self) -> Bounds {
        match self {
            Self::Cube | Self::Sphere => {
                Bounds::new(Point::new(-1.0, -1.0, -1.0), Point::new(1.0, 1.0, 1.0))
            }
            Self::Plane => Bounds::new(
                Point::new(f64::NEG_INFINITY, 0.0, f64::NEG_INFINITY),
                Point::new(f64::INFINITY, 0.0, f64::INFINITY),
            ),
            Self::Cylinder(conic) | Self::Cone(conic) => conic.bounds(),
            Self::Triangle { triangular, .. } => triangular.bounds(),
        }
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Triangular {
    pub p1: Point,
    pub p2: Point,
    pub p3: Point,
    e1: Vector,
    e2: Vector,
}

impl Triangular {
    pub fn new(p1: Point, p2: Point, p3: Point) -> Self {
        let e1 = p2 - p1;
        let e2 = p3 - p1;

        Self { p1, p2, p3, e1, e2 }
    }

    fn intersect(&self, ray: &Ray) -> Vec<f64> {
        let (ts, _, _) = self.intersect_uv(ray);
        ts
    }

    fn intersect_uv(&self, ray: &Ray) -> (Vec<f64>, f64, f64) {
        let dir_cross_e2 = ray.direction.cross(&self.e2);
        let det = self.e1.dot(&dir_cross_e2);

        let f = 1.0 / det;
        let p1_to_origin = ray.origin - self.p1;
        let u = f * p1_to_origin.dot(&dir_cross_e2);

        let origin_cross_e1 = p1_to_origin.cross(&self.e1);
        let v = f * ray.direction.dot(&origin_cross_e1);

        let ts = if det.abs().abs_diff_eq(&0.0, Point::default_epsilon())
            || u < 0.0
            || u > 1.0
            || v < 0.0
            || (u + v) > 1.0
        {
            Vec::new()
        } else {
            let t = f * self.e2.dot(&origin_cross_e1);

            vec![t]
        };

        (ts, u, v)
    }

    fn bounds(&self) -> Bounds {
        let mut bounds = Bounds::default();

        bounds.add_point(&self.p1);
        bounds.add_point(&self.p2);
        bounds.add_point(&self.p3);

        bounds
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Conic {
    minimum: f64,
    maximum: f64,
    closed: bool,
}

impl Conic {
    pub fn new(minimum: f64, maximum: f64, closed: bool) -> Self {
        Self {
            minimum,
            maximum,
            closed,
        }
    }

    fn normal_at(&self, object_point: &Point, y: f64) -> Vector {
        let dist = object_point.x.powi(2) + object_point.z.powi(2);

        if dist < 1.0 && object_point.y >= self.maximum - Point::default_epsilon() {
            Vector::new(0.0, 1.0, 0.0)
        } else if dist < 1.0 && object_point.y <= self.minimum + Point::default_epsilon() {
            Vector::new(0.0, -1.0, 0.0)
        } else {
            Vector::new(object_point.x, y, object_point.z)
        }
    }

    fn intersect(&self, ray: &Ray, a: f64, b: f64, c: f64, radius_at: fn(f64) -> f64) -> Vec<f64> {
        let discriminant = b.powi(2) - 4.0 * a * c;
        let mut xs = self.intersect_caps(ray, radius_at);

        if discriminant >= 0.0 {
            if a.abs_diff_eq(&0.0, Point::default_epsilon()) {
                if !b.abs_diff_eq(&0.0, Point::default_epsilon()) {
                    let t = -c / (2.0 * b);
                    xs.push(t);
                }
            } else {
                let t0 = (-b - discriminant.sqrt()) / (2.0 * a);
                let t1 = (-b + discriminant.sqrt()) / (2.0 * a);

                let y0 = ray.origin.y + t0 * ray.direction.y;
                if self.minimum < y0 && y0 < self.maximum {
                    xs.push(t0);
                }

                let y1 = ray.origin.y + t1 * ray.direction.y;
                if self.minimum < y1 && y1 < self.maximum {
                    xs.push(t1);
                }
            }
        }

        xs
    }

    fn intersect_caps(&self, ray: &Ray, radius_at: fn(f64) -> f64) -> Vec<f64> {
        let mut xs = Vec::new();

        if self.closed && !ray.direction.y.abs_diff_eq(&0.0, Vector::default_epsilon()) {
            let t = (self.minimum - ray.origin.y) / ray.direction.y;
            if Self::check_cap(ray, t, radius_at(self.minimum)) {
                xs.push(t);
            }

            let t = (self.maximum - ray.origin.y) / ray.direction.y;
            if Self::check_cap(ray, t, radius_at(self.maximum)) {
                xs.push(t);
            }
        }

        xs
    }

    fn check_cap(ray: &Ray, t: f64, radius: f64) -> bool {
        let x = ray.origin.x + t * ray.direction.x;
        let z = ray.origin.z + t * ray.direction.z;

        (x.powi(2) + z.powi(2)) <= radius.powi(2)
    }

    fn bounds(&self) -> Bounds {
        Bounds::new(
            Point::new(-1.0, self.minimum, -1.0),
            Point::new(1.0, self.maximum, 1.0),
        )
    }
}

impl Default for Conic {
    fn default() -> Self {
        Self {
            minimum: f64::NEG_INFINITY,
            maximum: f64::INFINITY,
            closed: false,
        }
    }
}

impl PartialEq for Conic {
    fn eq(&self, other: &Self) -> bool {
        self.abs_diff_eq(other, Self::default_epsilon())
    }
}

impl AbsDiffEq for Conic {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-5
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.minimum.abs_diff_eq(&other.minimum, epsilon)
            && self.maximum.abs_diff_eq(&other.maximum, epsilon)
            && self.closed == other.closed
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, FRAC_PI_4};

    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn default_transformation() {
        let shape = Shape::new_sphere();

        assert_eq!(shape.transform, Matrix::identity());
    }

    #[test]
    fn change_transformation() {
        let shape = Shape::new_sphere().with_transform(Matrix::translation(2.0, 3.0, 4.0));

        assert_eq!(shape.transform, Matrix::translation(2.0, 3.0, 4.0));
    }

    #[test]
    fn default_material() {
        let shape = Shape::new_sphere();

        assert_eq!(shape.material, Material::default());
    }

    #[test]
    fn assign_material() {
        let material = Material {
            ambient: 1.0,
            ..Material::default()
        };

        let shape = Shape::new_sphere().with_material(material);

        assert_eq!(shape.material, material);
    }

    #[test]
    fn glass_sphere() {
        let shape = Shape::new_glass_sphere();

        assert_eq!(shape.transform, Matrix::identity());
        assert_abs_diff_eq!(shape.material.transparency, 1.0);
        assert_abs_diff_eq!(shape.material.refractive_index, 1.5);
    }

    #[test]
    fn create_group() {
        let group = Shape::new_group(Vec::new());

        assert_eq!(group.transform, Matrix::identity());
        if let ShapeKind::Group(children) = group.kind {
            assert!(children.is_empty());
        } else {
            panic!("expected a group");
        }
    }

    #[test]
    fn create_group_with_children() {
        let sphere = Shape::new_sphere();
        let group = Shape::new_group(vec![sphere.clone()]);

        assert_eq!(group.transform, Matrix::identity());
        if let ShapeKind::Group(children) = group.kind {
            assert_eq!(children.len(), 1);
            assert_eq!(children[0], sphere);
        } else {
            panic!("expected a group");
        }
    }

    #[test]
    fn constructing_a_triangle() {
        let p1 = Point::new(0.0, 1.0, 0.0);
        let p2 = Point::new(-1.0, 0.0, 0.0);
        let p3 = Point::new(1.0, 0.0, 0.0);
        let triangle = Shape::new_triangle(p1, p2, p3);

        if let ShapeKind::Single(SingleKind::Triangle { triangular, normal }) = triangle.kind {
            assert_eq!(p1, triangular.p1);
            assert_eq!(p2, triangular.p2);
            assert_eq!(p3, triangular.p3);
            assert_eq!(triangular.e1, Vector::new(-1.0, -1.0, 0.0));
            assert_eq!(triangular.e2, Vector::new(1.0, -1.0, 0.0));
            assert_eq!(normal, Vector::new(0.0, 0.0, 1.0));
        } else {
            panic!("expected a triangle");
        }
    }

    #[test]
    fn intersect_scaled_sphere_with_ray() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let shape = Shape::new_sphere().with_transform(Matrix::scaling(2.0, 2.0, 2.0));

        let xs = shape
            .intersect(&r, im::Vector::new())
            .0
            .iter()
            .map(|i| i.t)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], 3.0);
        assert_abs_diff_eq!(xs[1], 7.0);
    }

    #[test]
    fn intersect_translated_sphere_with_ray() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let shape = Shape::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));

        assert!(shape.intersect(&r, im::Vector::new()).0.is_empty());
    }

    #[test]
    fn ray_intersects_sphere() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let s = SingleKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], 4.0);
        assert_abs_diff_eq!(xs[1], 6.0);
    }

    #[test]
    fn ray_tangent_to_sphere() {
        let r = Ray::new(Point::new(0.0, 1.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let s = SingleKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], 5.0);
        assert_abs_diff_eq!(xs[1], 5.0);
    }

    #[test]
    fn ray_misses_sphere() {
        let r = Ray::new(Point::new(0.0, 2.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let s = SingleKind::Sphere;

        let intersects = s.intersect(&r);

        assert!(intersects.is_empty());
    }

    #[test]
    fn ray_originates_inside_sphere() {
        let r = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));
        let s = SingleKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], -1.0);
        assert_abs_diff_eq!(xs[1], 1.0);
    }

    #[test]
    fn sphere_behind_ray() {
        let r = Ray::new(Point::new(0.0, 0.0, 5.0), Vector::new(0.0, 0.0, 1.0));
        let s = SingleKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], -6.0);
        assert_abs_diff_eq!(xs[1], -4.0);
    }

    #[test]
    fn intersecting_a_ray_parallel_to_the_plane() {
        let plane = SingleKind::Plane;
        let r = Ray::new(Point::new(0.0, 10.0, 0.0), Vector::new(0.0, 0.0, 1.0));

        let xs = plane.intersect(&r);

        assert!(xs.is_empty());
    }

    #[test]
    fn intersecting_a_coplanar_ray() {
        let plane = SingleKind::Plane;
        let r = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));

        let xs = plane.intersect(&r);

        assert!(xs.is_empty());
    }

    #[test]
    fn ray_intersecting_from_above() {
        let plane = SingleKind::Plane;
        let r = Ray::new(Point::new(0.0, 1.0, 0.0), Vector::new(0.0, -1.0, 0.0));

        let xs = plane.intersect(&r);

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 1.0);
    }

    #[test]
    fn ray_intersecting_from_below() {
        let plane = SingleKind::Plane;
        let r = Ray::new(Point::new(0.0, -1.0, 0.0), Vector::new(0.0, 1.0, 0.0));

        let xs = plane.intersect(&r);

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 1.0);
    }

    #[test]
    fn intersect_sets_shape() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let shape = Shape::new_sphere();

        let shapes = shape
            .intersect(&r, im::Vector::new())
            .0
            .iter()
            .map(|i| i.object)
            .collect::<Vec<_>>();

        assert_eq!(shapes[0], &shape);
        assert_eq!(shapes[1], &shape);
    }

    #[test]
    fn ray_intersects_cube() {
        let scenarios = [
            (
                Point::new(5.0, 0.5, 0.0),
                Vector::new(-1.0, 0.0, 0.0),
                4.0,
                6.0,
            ),
            (
                Point::new(-5.0, 0.5, 0.0),
                Vector::new(1.0, 0.0, 0.0),
                4.0,
                6.0,
            ),
            (
                Point::new(0.5, 5.0, 0.0),
                Vector::new(0.0, -1.0, 0.0),
                4.0,
                6.0,
            ),
            (
                Point::new(0.5, -5.0, 0.0),
                Vector::new(0.0, 1.0, 0.0),
                4.0,
                6.0,
            ),
            (
                Point::new(0.5, 0.0, 5.0),
                Vector::new(0.0, 0.0, -1.0),
                4.0,
                6.0,
            ),
            (
                Point::new(0.5, 0.0, -5.0),
                Vector::new(0.0, 0.0, 1.0),
                4.0,
                6.0,
            ),
            (
                Point::new(0.0, 0.5, 0.0),
                Vector::new(0.0, 0.0, 1.0),
                -1.0,
                1.0,
            ),
        ];

        for (origin, direction, t0, t1) in scenarios {
            let c = SingleKind::Cube;
            let r = Ray::new(origin, direction);

            let xs = c.intersect(&r);

            assert_eq!(xs.len(), 2);
            assert_abs_diff_eq!(xs[0], t0);
            assert_abs_diff_eq!(xs[1], t1);
        }
    }

    #[test]
    fn ray_misses_cube() {
        let scenarios = [
            (
                Point::new(-2.0, 0.0, 0.0),
                Vector::new(0.2673, 0.5345, 0.8018),
            ),
            (
                Point::new(0.0, -2.0, 0.0),
                Vector::new(0.8018, 0.2673, 0.5345),
            ),
            (
                Point::new(0.0, 0.0, -2.0),
                Vector::new(0.5345, 0.8018, 0.2673),
            ),
            (Point::new(2.0, 0.0, 2.0), Vector::new(0.0, 0.0, -1.0)),
            (Point::new(0.0, 2.0, 2.0), Vector::new(0.0, -1.0, 0.0)),
            (Point::new(2.0, 2.0, 0.0), Vector::new(-1.0, 0.0, 0.0)),
        ];

        for (origin, direction) in scenarios {
            let c = SingleKind::Cube;
            let r = Ray::new(origin, direction);

            let xs = c.intersect(&r);

            assert!(xs.is_empty());
        }
    }

    #[test]
    fn ray_misses_cylinder() {
        let scenarios = [
            (Point::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
            (Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
            (Point::new(0.0, 0.0, -5.0), Vector::new(1.0, 1.0, 1.0)),
        ];

        for (origin, direction) in scenarios {
            let c = SingleKind::Cylinder(Conic::default());
            let r = Ray::new(origin, direction);

            let xs = c.intersect(&r);

            assert!(xs.is_empty());
        }
    }

    #[test]
    fn ray_intersects_cylinder() {
        let scenarios = [
            (
                Point::new(1.0, 0.0, -5.0),
                Vector::new(0.0, 0.0, 1.0),
                5.0,
                5.0,
            ),
            (
                Point::new(0.0, 0.0, -5.0),
                Vector::new(0.0, 0.0, 1.0),
                4.0,
                6.0,
            ),
            (
                Point::new(0.5, 0.0, -5.0),
                Vector::new(0.1, 1.0, 1.0),
                6.80798,
                7.08872,
            ),
        ];

        for (origin, direction, t0, t1) in scenarios {
            let c = SingleKind::Cylinder(Conic::default());
            let r = Ray::new(origin, direction.normalize());

            let xs = c.intersect(&r);

            assert_eq!(xs.len(), 2);
            assert_abs_diff_eq!(xs[0], t0, epsilon = Point::default_epsilon());
            assert_abs_diff_eq!(xs[1], t1, epsilon = Point::default_epsilon());
        }
    }

    #[test]
    fn intersecting_constrained_cylinder() {
        let scenarios = [
            (Point::new(0.0, 1.5, 0.0), Vector::new(0.1, 1.0, 0.0), 0),
            (Point::new(0.0, 3.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
            (Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
            (Point::new(0.0, 2.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
            (Point::new(0.0, 1.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
            (Point::new(0.0, 1.5, -2.0), Vector::new(0.0, 0.0, 1.0), 2),
        ];

        for (point, direction, count) in scenarios {
            let cyl = SingleKind::Cylinder(Conic::new(1.0, 2.0, false));

            let r = Ray::new(point, direction.normalize());
            let xs = cyl.intersect(&r);

            assert_eq!(xs.len(), count);
        }
    }

    #[test]
    fn intersecting_caps_of_closed_cylinder() {
        let scenarios = [
            (Point::new(0.0, 3.0, 0.0), Vector::new(0.0, -1.0, 0.0), 2),
            (Point::new(0.0, 3.0, -2.0), Vector::new(0.0, -1.0, 2.0), 2),
            (Point::new(0.0, 4.0, -2.0), Vector::new(0.0, -1.0, 1.0), 2),
            (Point::new(0.0, 0.0, -2.0), Vector::new(0.0, 1.0, 2.0), 2),
            (Point::new(0.0, -1.0, -2.0), Vector::new(0.0, 1.0, 1.0), 2),
        ];

        for (point, direction, count) in scenarios {
            let cyl = SingleKind::Cylinder(Conic::new(1.0, 2.0, true));

            let r = Ray::new(point, direction.normalize());
            let xs = cyl.intersect(&r);

            assert_eq!(xs.len(), count);
        }
    }

    #[test]
    fn intersect_cone_with_ray() {
        let scenarios = [
            (
                Point::new(0.0, 0.0, -5.0),
                Vector::new(0.0, 0.0, 1.0),
                5.0,
                5.0,
            ),
            (
                Point::new(0.0, 0.0, -5.0),
                Vector::new(1.0, 1.0, 1.0),
                8.66025,
                8.66025,
            ),
            (
                Point::new(1.0, 1.0, -5.0),
                Vector::new(-0.5, -1.0, 1.0),
                4.55006,
                49.44994,
            ),
        ];

        for (origin, direction, t0, t1) in scenarios {
            let shape = SingleKind::Cone(Conic::default());

            let r = Ray::new(origin, direction.normalize());
            let xs = shape.intersect(&r);

            assert_eq!(xs.len(), 2);
            assert_abs_diff_eq!(xs[0], t0, epsilon = Point::default_epsilon());
            assert_abs_diff_eq!(xs[1], t1, epsilon = Point::default_epsilon());
        }
    }

    #[test]
    fn intersect_cone_with_ray_parallel_to_one_half() {
        let shape = SingleKind::Cone(Conic::default());

        let direction = Vector::new(0.0, 1.0, 1.0).normalize();
        let r = Ray::new(Point::new(0.0, 0.0, -1.0), direction);

        let xs = shape.intersect(&r);

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 0.35355, epsilon = Point::default_epsilon());
    }

    #[test]
    fn intersect_cone_end_caps() {
        let scenarios = [
            (Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 1.0, 0.0), 0),
            (Point::new(0.0, 0.0, -0.25), Vector::new(0.0, 1.0, 1.0), 2),
            (Point::new(0.0, 0.0, -0.25), Vector::new(0.0, 1.0, 0.0), 4),
        ];

        for (origin, direction, count) in scenarios {
            let shape = SingleKind::Cone(Conic::new(-0.5, 0.5, true));

            let r = Ray::new(origin, direction.normalize());
            let xs = shape.intersect(&r);

            assert_eq!(xs.len(), count);
        }
    }

    #[test]
    fn intersect_ray_with_empty_group() {
        let group = Shape::new_group(Vec::new());

        let r = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));
        let xs = group.intersect(&r, im::Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn intersect_ray_with_group_of_shapes() {
        let s1 = Shape::new_sphere();
        let s2 = Shape::new_sphere().with_transform(Matrix::translation(0.0, 0.0, -3.0));
        let s3 = Shape::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));
        let group = Shape::new_group(vec![s1.clone(), s2.clone(), s3]);

        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let xs = group
            .intersect(&r, im::Vector::new())
            .0
            .iter()
            .map(|i| i.object)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 4);
        assert_eq!(xs[0], &s2);
        assert_eq!(xs[1], &s2);
        assert_eq!(xs[2], &s1);
        assert_eq!(xs[3], &s1);
    }

    #[test]
    fn intersect_ray_with_transformed_group() {
        let s = Shape::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));
        let group = Shape::new_group(vec![s]).with_transform(Matrix::scaling(2.0, 2.0, 2.0));

        let r = Ray::new(Point::new(10.0, 0.0, -10.0), Vector::new(0.0, 0.0, 1.0));
        let xs = group.intersect(&r, im::Vector::new());

        assert_eq!(xs.0.len(), 2);
    }

    #[test]
    fn intersect_ray_parallel_to_triangle() {
        let t = Shape::new_triangle(
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Point::new(0.0, -1.0, -2.0), Vector::new(0.0, 1.0, 0.0));

        let xs = t.intersect(&r, im::Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_misses_p1_p3_edge() {
        let t = Shape::new_triangle(
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Point::new(1.0, 1.0, -2.0), Vector::new(0.0, 0.0, 1.0));

        let xs = t.intersect(&r, im::Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_misses_p1_p2_edge() {
        let t = Shape::new_triangle(
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Point::new(-1.0, 1.0, -2.0), Vector::new(0.0, 0.0, 1.0));

        let xs = t.intersect(&r, im::Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_misses_p2_p3_edge() {
        let t = Shape::new_triangle(
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Point::new(0.0, -1.0, -2.0), Vector::new(0.0, 0.0, 1.0));

        let xs = t.intersect(&r, im::Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_strikes_triangle() {
        let t = Shape::new_triangle(
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Point::new(0.0, 0.5, -2.0), Vector::new(0.0, 0.0, 1.0));

        let xs = t
            .intersect(&r, im::Vector::new())
            .0
            .iter()
            .map(|i| i.t)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 2.0);
    }

    #[test]
    fn intersection_with_smooth_triangle_stores_uv() {
        let p1 = Point::new(0.0, 1.0, 0.0);
        let p2 = Point::new(-1.0, 0.0, 0.0);
        let p3 = Point::new(1.0, 0.0, 0.0);
        let n1 = Vector::new(0.0, 1.0, 0.0);
        let n2 = Vector::new(-1.0, 0.0, 0.0);
        let n3 = Vector::new(1.0, 0.0, 0.0);

        let tri = Shape::new_smooth_triangle(p1, p2, p3, n1, n2, n3);
        let r = Ray::new(Point::new(-0.2, 0.3, -2.0), Vector::new(0.0, 0.0, 1.0));
        let intersections = tri.intersect(&r, im::Vector::new());
        let xs = intersections.0.iter().collect::<Vec<_>>();

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0].u, 0.45);
        assert_abs_diff_eq!(xs[0].v, 0.25);
    }

    #[test]
    fn normal_on_translated_sphere() {
        let shape = Shape::new_sphere().with_transform(Matrix::translation(0.0, 1.0, 0.0));
        let n = shape.normal_at(
            &Point::new(0.0, 1.70711, -FRAC_1_SQRT_2),
            &Intersection::new(0.0, &shape, im::Vector::new()),
            &im::Vector::new(),
        );

        assert_abs_diff_eq!(n, Vector::new(0.0, FRAC_1_SQRT_2, -FRAC_1_SQRT_2));
    }

    #[test]
    fn normal_on_transformed_sphere() {
        let shape = Shape::new_sphere()
            .with_transform(Matrix::scaling(1.0, 0.5, 1.0) * Matrix::rotation_z(FRAC_PI_4));
        let n = shape.normal_at(
            &Point::new(0.0, 2.0_f64.sqrt() / 2.0, -(2.0_f64.sqrt()) / 2.0),
            &Intersection::new(0.0, &shape, im::Vector::new()),
            &im::Vector::new(),
        );

        assert_abs_diff_eq!(n, Vector::new(0.0, 0.97014, -0.24254));
    }

    #[test]
    fn sphere_normal_x_axis() {
        let s = SingleKind::Sphere;
        let n = s.normal_at(&Point::new(1.0, 0.0, 0.0));

        assert_abs_diff_eq!(n, Vector::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn sphere_normal_y_axis() {
        let s = SingleKind::Sphere;
        let n = s.normal_at(&Point::new(0.0, 1.0, 0.0));

        assert_abs_diff_eq!(n, Vector::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn sphere_normal_z_axis() {
        let s = SingleKind::Sphere;
        let n = s.normal_at(&Point::new(0.0, 0.0, 1.0));

        assert_abs_diff_eq!(n, Vector::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn sphere_normal_non_axial_point() {
        let s = SingleKind::Sphere;
        let n = s.normal_at(&Point::new(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
        ));

        assert_abs_diff_eq!(
            n,
            Vector::new(
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0
            )
        );
    }

    #[test]
    fn normal_is_normalized() {
        let s = SingleKind::Sphere;
        let n = s.normal_at(&Point::new(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
        ));

        assert_abs_diff_eq!(n, n.normalize());
    }

    #[test]
    fn normal_of_plane_is_constant() {
        let plane = SingleKind::Plane;
        let n1 = plane.normal_at(&Point::new(0.0, 0.0, 0.0));
        let n2 = plane.normal_at(&Point::new(10.0, 0.0, -10.0));
        let n3 = plane.normal_at(&Point::new(-5.0, 0.0, 150.0));

        assert_abs_diff_eq!(n1, Vector::new(0.0, 1.0, 0.0));
        assert_abs_diff_eq!(n2, Vector::new(0.0, 1.0, 0.0));
        assert_abs_diff_eq!(n3, Vector::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn normal_on_surface_of_cube() {
        let scenarios = [
            (Point::new(1.0, 0.5, -0.8), Vector::new(1.0, 0.0, 0.0)),
            (Point::new(-1.0, -0.2, 0.9), Vector::new(-1.0, 0.0, 0.0)),
            (Point::new(-0.4, 1.0, -0.1), Vector::new(0.0, 1.0, 0.0)),
            (Point::new(0.3, -1.0, -0.7), Vector::new(0.0, -1.0, 0.0)),
            (Point::new(-0.6, 0.3, 1.0), Vector::new(0.0, 0.0, 1.0)),
            (Point::new(0.4, 0.4, -1.0), Vector::new(0.0, 0.0, -1.0)),
            (Point::new(1.0, 1.0, 1.0), Vector::new(1.0, 0.0, 0.0)),
            (Point::new(-1.0, -1.0, -1.0), Vector::new(-1.0, 0.0, 0.0)),
        ];

        for (point, normal) in scenarios {
            let c = SingleKind::Cube;
            let n = c.normal_at(&point);

            assert_abs_diff_eq!(n, normal);
        }
    }

    #[test]
    fn normal_on_cylinder() {
        let scenarios = [
            (Point::new(1.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
            (Point::new(0.0, 5.0, -1.0), Vector::new(0.0, 0.0, -1.0)),
            (Point::new(0.0, -2.0, 1.0), Vector::new(0.0, 0.0, 1.0)),
            (Point::new(-1.0, 1.0, 0.0), Vector::new(-1.0, 0.0, 0.0)),
        ];

        for (point, normal) in scenarios {
            let c = SingleKind::Cylinder(Conic::default());
            let n = c.normal_at(&point);

            assert_abs_diff_eq!(n, normal);
        }
    }

    #[test]
    fn normal_on_cylinder_end_caps() {
        let scenarios = [
            (Point::new(0.0, 1.0, 0.03), Vector::new(0.0, -1.0, 0.0)),
            (Point::new(0.5, 1.0, 0.0), Vector::new(0.0, -1.0, 0.0)),
            (Point::new(0.0, 1.0, 0.5), Vector::new(0.0, -1.0, 0.0)),
            (Point::new(0.0, 2.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
            (Point::new(0.5, 2.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
            (Point::new(0.0, 2.0, 0.5), Vector::new(0.0, 1.0, 0.0)),
        ];

        for (point, normal) in scenarios {
            let c = SingleKind::Cylinder(Conic::new(1.0, 2.0, true));
            let n = c.normal_at(&point);

            assert_abs_diff_eq!(n, normal);
        }
    }

    #[test]
    fn normal_on_cone() {
        let scenarios = [
            (Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 0.0)),
            (
                Point::new(1.0, 1.0, 1.0),
                Vector::new(1.0, -(2.0_f64.sqrt()), 1.0),
            ),
            (Point::new(-1.0, -1.0, 0.0), Vector::new(-1.0, 1.0, 0.0)),
        ];

        for (point, normal) in scenarios {
            let c = SingleKind::Cone(Conic::default());
            let n = c.normal_at(&point);

            assert_abs_diff_eq!(n, normal);
        }
    }

    #[test]
    fn normal_on_child_object() {
        let s = Shape::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));
        let g2 = Shape::new_group(vec![s.clone()]).with_transform(Matrix::scaling(1.0, 2.0, 3.0));
        let g1 = Shape::new_group(vec![g2.clone()]).with_transform(Matrix::rotation_y(FRAC_PI_2));

        let n = s.normal_at(
            &Point::new(1.7321, 1.1547, -5.5774),
            &Intersection::new(0.0, &s, im::Vector::new()),
            &im::vector![g2.transform_inversed, g1.transform_inversed],
        );

        assert_abs_diff_eq!(n, Vector::new(0.2857, 0.4286, -0.8571));
    }

    #[test]
    fn normal_on_triangle() {
        let t = Shape::new_triangle(
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        );

        assert_eq!(
            t.normal_at(
                &Point::new(0.0, 0.5, 0.0),
                &Intersection::new(0.0, &t, im::Vector::new()),
                &im::Vector::new()
            ),
            Vector::new(0.0, 0.0, 1.0)
        );
        assert_eq!(
            t.normal_at(
                &Point::new(-0.5, 0.75, 0.0),
                &Intersection::new(0.0, &t, im::Vector::new()),
                &im::Vector::new()
            ),
            Vector::new(0.0, 0.0, 1.0)
        );
        assert_eq!(
            t.normal_at(
                &Point::new(0.5, 0.25, 0.0),
                &Intersection::new(0.0, &t, im::Vector::new()),
                &im::Vector::new()
            ),
            Vector::new(0.0, 0.0, 1.0)
        );
    }

    #[test]
    fn smooth_triangle_uses_uv_to_interpolate_normal() {
        let p1 = Point::new(0.0, 1.0, 0.0);
        let p2 = Point::new(-1.0, 0.0, 0.0);
        let p3 = Point::new(1.0, 0.0, 0.0);
        let n1 = Vector::new(0.0, 1.0, 0.0);
        let n2 = Vector::new(-1.0, 0.0, 0.0);
        let n3 = Vector::new(1.0, 0.0, 0.0);

        let tri = Shape::new_smooth_triangle(p1, p2, p3, n1, n2, n3);
        let i = Intersection::new_with_uv(1.0, &tri, im::Vector::new(), 0.45, 0.25);
        let n = tri.normal_at(&Point::new(0.0, 0.0, 0.0), &i, &im::Vector::new());

        assert_eq!(n, Vector::new(-0.5547, 0.83205, 0.0));
    }

    #[test]
    fn convert_point_world_to_object() {
        let s = Shape::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));
        let g2 = Shape::new_group(vec![s.clone()]).with_transform(Matrix::scaling(2.0, 2.0, 2.0));
        let g1 = Shape::new_group(vec![g2.clone()]).with_transform(Matrix::rotation_y(FRAC_PI_2));

        assert_eq!(
            s.world_to_object(
                &Point::new(-2.0, 0.0, -10.0),
                &im::vector![g2.transform_inversed, g1.transform_inversed]
            ),
            Point::new(0.0, 0.0, -1.0)
        );
    }

    #[test]
    fn convert_normal_object_to_world() {
        let s = Shape::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));
        let g2 = Shape::new_group(vec![s.clone()]).with_transform(Matrix::scaling(1.0, 2.0, 3.0));
        let g1 = Shape::new_group(vec![g2.clone()]).with_transform(Matrix::rotation_y(FRAC_PI_2));

        let n = s.normal_to_world(
            &Vector::new(
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0,
            ),
            &im::vector![g2.transform_inversed, g1.transform_inversed],
        );

        assert_eq!(n, Vector::new(0.2857, 0.4286, -0.8571));
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
            let s1 = Shape::new_sphere();
            let s2 = Shape::new_cube();
            let xs = Intersections::new([
                Intersection::new(1.0, &s1, im::Vector::new()),
                Intersection::new(2.0, &s2, im::Vector::new()),
                Intersection::new(3.0, &s1, im::Vector::new()),
                Intersection::new(4.0, &s2, im::Vector::new()),
            ]);

            let result = operation
                .filter_intersections(&s1, xs.clone())
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
        let c = Shape::new_csg(Operation::Union, Shape::new_sphere(), Shape::new_cube());
        let r = Ray::new(Point::new(0.0, 2.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let xs = c.intersect(&r, im::Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_hits_csg_object() {
        let s1 = Shape::new_sphere();
        let s2 = Shape::new_sphere().with_transform(Matrix::translation(0.0, 0.0, 0.5));

        let c = Shape::new_csg(Operation::Union, s1.clone(), s2.clone());
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let intersections = c.intersect(&r, im::Vector::new());
        let xs = intersections.0.iter().collect::<Vec<_>>();

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0].t, 4.0);
        assert_eq!(*xs[0].object, s1);
        assert_abs_diff_eq!(xs[1].t, 6.5);
        assert_eq!(*xs[1].object, s2);
    }
}
