use std::collections::BTreeSet;

use approx::AbsDiffEq;
use im::Vector;

use crate::core::{matrix::Matrix, tuple::Tuple};

use super::{
    bounds::Bounds,
    intersection::{Intersection, Intersections},
    material::Material,
    ray::Ray,
};

#[derive(Debug, PartialEq, Clone)]
pub struct Shape {
    pub transform: Matrix<4>,
    transform_inversed: Matrix<4>,
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

    pub fn new_triangle(p1: Tuple, p2: Tuple, p3: Tuple) -> Self {
        let triangular = Triangular::new(p1, p2, p3);
        let normal = triangular.e1.cross(&triangular.e2).normalize();
        Self::new(ShapeKind::Single(SingleKind::Triangle {
            triangular,
            normal,
        }))
    }

    pub fn new_smooth_triangle(
        p1: Tuple,
        p2: Tuple,
        p3: Tuple,
        n1: Tuple,
        n2: Tuple,
        n3: Tuple,
    ) -> Self {
        let triangular = Triangular::new(p1, p2, p3);

        Self::new(ShapeKind::SmoothTriangle {
            triangular,
            n1,
            n2,
            n3,
        })
    }

    pub fn with_transform(mut self, transform: Matrix<4>) -> Self {
        self.transform = transform;
        self.transform_inversed = transform.inverse();
        self
    }

    pub fn with_material(mut self, material: Material) -> Self {
        self.material = material;
        self
    }

    pub fn intersect(&self, ray: &Ray, mut trail: Vector<Matrix<4>>) -> Intersections {
        let local_ray = ray.transform(&self.transform_inversed);

        match &self.kind {
            ShapeKind::Single(single) => {
                let ts = single.intersect(&local_ray);
                let intersections = ts
                    .iter()
                    .map(|t| Intersection::new(*t, self, trail.clone()));

                Intersections(BTreeSet::from_iter(intersections))
            }
            ShapeKind::Group(children) => {
                trail.push_front(self.transform);
                if self.bounds().intersects(&local_ray) {
                    let intersections = children
                        .iter()
                        .flat_map(|child| child.intersect(&local_ray, trail.clone()).0);

                    Intersections(BTreeSet::from_iter(intersections))
                } else {
                    Intersections(BTreeSet::new())
                }
            }
            ShapeKind::SmoothTriangle { triangular, .. } => {
                // TODO: kind of duplicated with single case
                let (ts, u, v) = triangular.intersect_uv(&local_ray);
                let intersections = ts
                    .iter()
                    .map(|t| Intersection::new_with_uv(*t, self, trail.clone(), u, v));

                Intersections(BTreeSet::from_iter(intersections))
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
        }
    }

    // TODO: hit only used for smooth triangles
    pub fn normal_at(&self, point: &Tuple, hit: &Intersection, trail: &Vector<Matrix<4>>) -> Tuple {
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
            ShapeKind::Group(_) => {
                panic!("groups should not have normal vectors");
            }
        }
    }

    // TODO: should mutate vector to avoid duplication?
    pub fn world_to_object(&self, world_point: &Tuple, trail: &Vector<Matrix<4>>) -> Tuple {
        let trail_point = trail
            .iter()
            .rev()
            .fold(*world_point, |acc, mat| mat.inverse() * acc);

        self.transform_inversed * trail_point
    }

    fn normal_to_world(&self, normal: &Tuple, trail: &Vector<Matrix<4>>) -> Tuple {
        let mut normal = self.transform_inversed.transpose() * *normal;
        normal.w = 0.0;
        let normal = normal.normalize();

        trail.iter().fold(normal, |acc, mat| {
            let mut normal = mat.inverse().transpose() * acc;
            normal.w = 0.0;
            normal.normalize()
        })
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
        n1: Tuple,
        n2: Tuple,
        n3: Tuple,
    },
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
        normal: Tuple,
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
        let sphere_to_ray = ray.origin - Tuple::point(0.0, 0.0, 0.0);

        let a = ray.direction.dot(&ray.direction);
        let b = 2.0 * sphere_to_ray.dot(&ray.direction);
        let c = sphere_to_ray.dot(&sphere_to_ray) - 1.0;

        let discriminant = b * b - 4.0 * a * c;

        if discriminant < 0.0 {
            Vec::new()
        } else {
            let t1 = (-b - discriminant.sqrt()) / (2.0 * a);
            let t2 = (-b + discriminant.sqrt()) / (2.0 * a);

            vec![t1, t2]
        }
    }

    fn plane_intersect(ray: &Ray) -> Vec<f64> {
        if ray.direction.y.abs() < Tuple::default_epsilon() {
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

    pub fn normal_at(&self, point: &Tuple) -> Tuple {
        match self {
            Self::Sphere => Self::sphere_normal_at(point),
            Self::Plane => Self::plane_normal_at(),
            Self::Cube => Self::cube_normal_at(point),
            Self::Cylinder(conic) => Self::cylinder_normal_at(point, conic),
            Self::Cone(conic) => Self::cone_normal_at(point, conic),
            Self::Triangle { normal, .. } => *normal,
        }
    }

    fn sphere_normal_at(object_point: &Tuple) -> Tuple {
        *object_point - Tuple::point(0.0, 0.0, 0.0)
    }

    fn plane_normal_at() -> Tuple {
        Tuple::vector(0.0, 1.0, 0.0)
    }

    fn cube_normal_at(object_point: &Tuple) -> Tuple {
        let max_c = object_point
            .x
            .abs()
            .max(object_point.y.abs())
            .max(object_point.z.abs());

        if max_c == object_point.x.abs() {
            Tuple::vector(object_point.x, 0.0, 0.0)
        } else if max_c == object_point.y.abs() {
            Tuple::vector(0.0, object_point.y, 0.0)
        } else {
            Tuple::vector(0.0, 0.0, object_point.z)
        }
    }

    fn cylinder_normal_at(object_point: &Tuple, conic: &Conic) -> Tuple {
        conic.normal_at(object_point, 0.0)
    }

    fn cone_normal_at(object_point: &Tuple, conic: &Conic) -> Tuple {
        let y = (object_point.x.powi(2) + object_point.z.powi(2)).sqrt();
        let y = if object_point.y > 0.0 { -y } else { y };

        conic.normal_at(object_point, y)
    }

    fn bounds(&self) -> Bounds {
        match self {
            Self::Sphere => {
                Bounds::new(Tuple::point(-1.0, -1.0, -1.0), Tuple::point(1.0, 1.0, 1.0))
            }
            Self::Plane => Bounds::new(
                Tuple::point(f64::NEG_INFINITY, 0.0, f64::NEG_INFINITY),
                Tuple::point(f64::INFINITY, 0.0, f64::INFINITY),
            ),
            Self::Cube => Bounds::new(Tuple::point(-1.0, -1.0, -1.0), Tuple::point(1.0, 1.0, 1.0)),
            Self::Cylinder(conic) => conic.bounds(),
            Self::Cone(conic) => conic.bounds(),
            Self::Triangle { triangular, .. } => triangular.bounds(),
        }
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Triangular {
    pub p1: Tuple,
    pub p2: Tuple,
    pub p3: Tuple,
    e1: Tuple,
    e2: Tuple,
}

impl Triangular {
    pub fn new(p1: Tuple, p2: Tuple, p3: Tuple) -> Self {
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

        let ts = if det.abs().abs_diff_eq(&0.0, Tuple::default_epsilon())
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

    fn normal_at(&self, object_point: &Tuple, y: f64) -> Tuple {
        let dist = object_point.x.powi(2) + object_point.z.powi(2);

        if dist < 1.0 && object_point.y >= self.maximum - Tuple::default_epsilon() {
            Tuple::vector(0.0, 1.0, 0.0)
        } else if dist < 1.0 && object_point.y <= self.minimum + Tuple::default_epsilon() {
            Tuple::vector(0.0, -1.0, 0.0)
        } else {
            Tuple::vector(object_point.x, y, object_point.z)
        }
    }

    fn intersect(&self, ray: &Ray, a: f64, b: f64, c: f64, radius_at: fn(f64) -> f64) -> Vec<f64> {
        let discriminant = b.powi(2) - 4.0 * a * c;
        let mut xs = self.intersect_caps(ray, radius_at);

        if discriminant >= 0.0 {
            if a.abs_diff_eq(&0.0, Tuple::default_epsilon()) {
                if !b.abs_diff_eq(&0.0, Tuple::default_epsilon()) {
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

        if self.closed && !ray.direction.y.abs_diff_eq(&0.0, Tuple::default_epsilon()) {
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
            Tuple::point(-1.0, self.minimum, -1.0),
            Tuple::point(1.0, self.maximum, 1.0),
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
    use im::vector;

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
        assert_eq!(shape.material.transparency, 1.0);
        assert_eq!(shape.material.refractive_index, 1.5);
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
        let p1 = Tuple::point(0.0, 1.0, 0.0);
        let p2 = Tuple::point(-1.0, 0.0, 0.0);
        let p3 = Tuple::point(1.0, 0.0, 0.0);
        let triangle = Shape::new_triangle(p1, p2, p3);

        if let ShapeKind::Single(SingleKind::Triangle { triangular, normal }) = triangle.kind {
            assert_eq!(p1, triangular.p1);
            assert_eq!(p2, triangular.p2);
            assert_eq!(p3, triangular.p3);
            assert_eq!(triangular.e1, Tuple::vector(-1.0, -1.0, 0.0));
            assert_eq!(triangular.e2, Tuple::vector(1.0, -1.0, 0.0));
            assert_eq!(normal, Tuple::vector(0.0, 0.0, 1.0));
        } else {
            panic!("expected a triangle");
        }
    }

    #[test]
    fn intersect_scaled_sphere_with_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Shape::new_sphere().with_transform(Matrix::scaling(2.0, 2.0, 2.0));

        let xs = shape
            .intersect(&r, Vector::new())
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
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Shape::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));

        assert!(shape.intersect(&r, Vector::new()).0.is_empty());
    }

    #[test]
    fn ray_intersects_sphere() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = SingleKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], 4.0);
        assert_abs_diff_eq!(xs[1], 6.0);
    }

    #[test]
    fn ray_tangent_to_sphere() {
        let r = Ray::new(Tuple::point(0.0, 1.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = SingleKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], 5.0);
        assert_abs_diff_eq!(xs[1], 5.0);
    }

    #[test]
    fn ray_misses_sphere() {
        let r = Ray::new(Tuple::point(0.0, 2.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = SingleKind::Sphere;

        let intersects = s.intersect(&r);

        assert!(intersects.is_empty());
    }

    #[test]
    fn ray_originates_inside_sphere() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = SingleKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], -1.0);
        assert_abs_diff_eq!(xs[1], 1.0);
    }

    #[test]
    fn sphere_behind_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = SingleKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], -6.0);
        assert_abs_diff_eq!(xs[1], -4.0);
    }

    #[test]
    fn intersecting_a_ray_parallel_to_the_plane() {
        let plane = SingleKind::Plane;
        let r = Ray::new(Tuple::point(0.0, 10.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = plane.intersect(&r);

        assert!(xs.is_empty());
    }

    #[test]
    fn intersecting_a_coplanar_ray() {
        let plane = SingleKind::Plane;
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = plane.intersect(&r);

        assert!(xs.is_empty());
    }

    #[test]
    fn ray_intersecting_from_above() {
        let plane = SingleKind::Plane;
        let r = Ray::new(Tuple::point(0.0, 1.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));

        let xs = plane.intersect(&r);

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 1.0);
    }

    #[test]
    fn ray_intersecting_from_below() {
        let plane = SingleKind::Plane;
        let r = Ray::new(Tuple::point(0.0, -1.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));

        let xs = plane.intersect(&r);

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 1.0);
    }

    #[test]
    fn intersect_sets_shape() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Shape::new_sphere();

        let shapes = shape
            .intersect(&r, Vector::new())
            .0
            .iter()
            .map(|i| i.shape)
            .collect::<Vec<_>>();

        assert_eq!(shapes[0], &shape);
        assert_eq!(shapes[1], &shape);
    }

    #[test]
    fn ray_intersects_cube() {
        let scenarios = [
            (
                Tuple::point(5.0, 0.5, 0.0),
                Tuple::vector(-1.0, 0.0, 0.0),
                4.0,
                6.0,
            ),
            (
                Tuple::point(-5.0, 0.5, 0.0),
                Tuple::vector(1.0, 0.0, 0.0),
                4.0,
                6.0,
            ),
            (
                Tuple::point(0.5, 5.0, 0.0),
                Tuple::vector(0.0, -1.0, 0.0),
                4.0,
                6.0,
            ),
            (
                Tuple::point(0.5, -5.0, 0.0),
                Tuple::vector(0.0, 1.0, 0.0),
                4.0,
                6.0,
            ),
            (
                Tuple::point(0.5, 0.0, 5.0),
                Tuple::vector(0.0, 0.0, -1.0),
                4.0,
                6.0,
            ),
            (
                Tuple::point(0.5, 0.0, -5.0),
                Tuple::vector(0.0, 0.0, 1.0),
                4.0,
                6.0,
            ),
            (
                Tuple::point(0.0, 0.5, 0.0),
                Tuple::vector(0.0, 0.0, 1.0),
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
                Tuple::point(-2.0, 0.0, 0.0),
                Tuple::vector(0.2673, 0.5345, 0.8018),
            ),
            (
                Tuple::point(0.0, -2.0, 0.0),
                Tuple::vector(0.8018, 0.2673, 0.5345),
            ),
            (
                Tuple::point(0.0, 0.0, -2.0),
                Tuple::vector(0.5345, 0.8018, 0.2673),
            ),
            (Tuple::point(2.0, 0.0, 2.0), Tuple::vector(0.0, 0.0, -1.0)),
            (Tuple::point(0.0, 2.0, 2.0), Tuple::vector(0.0, -1.0, 0.0)),
            (Tuple::point(2.0, 2.0, 0.0), Tuple::vector(-1.0, 0.0, 0.0)),
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
            (Tuple::point(1.0, 0.0, 0.0), Tuple::vector(0.0, 1.0, 0.0)),
            (Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 1.0, 0.0)),
            (Tuple::point(0.0, 0.0, -5.0), Tuple::vector(1.0, 1.0, 1.0)),
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
                Tuple::point(1.0, 0.0, -5.0),
                Tuple::vector(0.0, 0.0, 1.0),
                5.0,
                5.0,
            ),
            (
                Tuple::point(0.0, 0.0, -5.0),
                Tuple::vector(0.0, 0.0, 1.0),
                4.0,
                6.0,
            ),
            (
                Tuple::point(0.5, 0.0, -5.0),
                Tuple::vector(0.1, 1.0, 1.0),
                6.80798,
                7.08872,
            ),
        ];

        for (origin, direction, t0, t1) in scenarios {
            let c = SingleKind::Cylinder(Conic::default());
            let r = Ray::new(origin, direction.normalize());

            let xs = c.intersect(&r);

            assert_eq!(xs.len(), 2);
            assert_abs_diff_eq!(xs[0], t0, epsilon = Tuple::default_epsilon());
            assert_abs_diff_eq!(xs[1], t1, epsilon = Tuple::default_epsilon());
        }
    }

    #[test]
    fn intersecting_constrained_cylinder() {
        let scenarios = [
            (Tuple::point(0.0, 1.5, 0.0), Tuple::vector(0.1, 1.0, 0.0), 0),
            (
                Tuple::point(0.0, 3.0, -5.0),
                Tuple::vector(0.0, 0.0, 1.0),
                0,
            ),
            (
                Tuple::point(0.0, 0.0, -5.0),
                Tuple::vector(0.0, 0.0, 1.0),
                0,
            ),
            (
                Tuple::point(0.0, 2.0, -5.0),
                Tuple::vector(0.0, 0.0, 1.0),
                0,
            ),
            (
                Tuple::point(0.0, 1.0, -5.0),
                Tuple::vector(0.0, 0.0, 1.0),
                0,
            ),
            (
                Tuple::point(0.0, 1.5, -2.0),
                Tuple::vector(0.0, 0.0, 1.0),
                2,
            ),
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
            (
                Tuple::point(0.0, 3.0, 0.0),
                Tuple::vector(0.0, -1.0, 0.0),
                2,
            ),
            (
                Tuple::point(0.0, 3.0, -2.0),
                Tuple::vector(0.0, -1.0, 2.0),
                2,
            ),
            (
                Tuple::point(0.0, 4.0, -2.0),
                Tuple::vector(0.0, -1.0, 1.0),
                2,
            ),
            (
                Tuple::point(0.0, 0.0, -2.0),
                Tuple::vector(0.0, 1.0, 2.0),
                2,
            ),
            (
                Tuple::point(0.0, -1.0, -2.0),
                Tuple::vector(0.0, 1.0, 1.0),
                2,
            ),
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
                Tuple::point(0.0, 0.0, -5.0),
                Tuple::vector(0.0, 0.0, 1.0),
                5.0,
                5.0,
            ),
            (
                Tuple::point(0.0, 0.0, -5.0),
                Tuple::vector(1.0, 1.0, 1.0),
                8.66025,
                8.66025,
            ),
            (
                Tuple::point(1.0, 1.0, -5.0),
                Tuple::vector(-0.5, -1.0, 1.0),
                4.55006,
                49.44994,
            ),
        ];

        for (origin, direction, t0, t1) in scenarios {
            let shape = SingleKind::Cone(Conic::default());

            let r = Ray::new(origin, direction.normalize());
            let xs = shape.intersect(&r);

            assert_eq!(xs.len(), 2);
            assert_abs_diff_eq!(xs[0], t0, epsilon = Tuple::default_epsilon());
            assert_abs_diff_eq!(xs[1], t1, epsilon = Tuple::default_epsilon());
        }
    }

    #[test]
    fn intersect_cone_with_ray_parallel_to_one_half() {
        let shape = SingleKind::Cone(Conic::default());

        let direction = Tuple::vector(0.0, 1.0, 1.0).normalize();
        let r = Ray::new(Tuple::point(0.0, 0.0, -1.0), direction);

        let xs = shape.intersect(&r);

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 0.35355, epsilon = Tuple::default_epsilon());
    }

    #[test]
    fn intersect_cone_end_caps() {
        let scenarios = [
            (
                Tuple::point(0.0, 0.0, -5.0),
                Tuple::vector(0.0, 1.0, 0.0),
                0,
            ),
            (
                Tuple::point(0.0, 0.0, -0.25),
                Tuple::vector(0.0, 1.0, 1.0),
                2,
            ),
            (
                Tuple::point(0.0, 0.0, -0.25),
                Tuple::vector(0.0, 1.0, 0.0),
                4,
            ),
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

        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = group.intersect(&r, Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn intersect_ray_with_group_of_shapes() {
        let s1 = Shape::new_sphere();
        let s2 = Shape::new_sphere().with_transform(Matrix::translation(0.0, 0.0, -3.0));
        let s3 = Shape::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));
        let group = Shape::new_group(vec![s1.clone(), s2.clone(), s3]);

        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = group
            .intersect(&r, Vector::new())
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
        let s = Shape::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));
        let group = Shape::new_group(vec![s]).with_transform(Matrix::scaling(2.0, 2.0, 2.0));

        let r = Ray::new(Tuple::point(10.0, 0.0, -10.0), Tuple::vector(0.0, 0.0, 1.0));
        let xs = group.intersect(&r, Vector::new());

        assert_eq!(xs.0.len(), 2);
    }

    #[test]
    fn intersect_ray_parallel_to_triangle() {
        let t = Shape::new_triangle(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Tuple::point(0.0, -1.0, -2.0), Tuple::vector(0.0, 1.0, 0.0));

        let xs = t.intersect(&r, Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_misses_p1_p3_edge() {
        let t = Shape::new_triangle(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Tuple::point(1.0, 1.0, -2.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = t.intersect(&r, Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_misses_p1_p2_edge() {
        let t = Shape::new_triangle(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Tuple::point(-1.0, 1.0, -2.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = t.intersect(&r, Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_misses_p2_p3_edge() {
        let t = Shape::new_triangle(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Tuple::point(0.0, -1.0, -2.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = t.intersect(&r, Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_strikes_triangle() {
        let t = Shape::new_triangle(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Tuple::point(0.0, 0.5, -2.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = t
            .intersect(&r, Vector::new())
            .0
            .iter()
            .map(|i| i.t)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 1);
        assert_eq!(xs[0], 2.0);
    }

    #[test]
    fn intersection_with_smooth_triangle_stores_uv() {
        let p1 = Tuple::point(0.0, 1.0, 0.0);
        let p2 = Tuple::point(-1.0, 0.0, 0.0);
        let p3 = Tuple::point(1.0, 0.0, 0.0);
        let n1 = Tuple::vector(0.0, 1.0, 0.0);
        let n2 = Tuple::vector(-1.0, 0.0, 0.0);
        let n3 = Tuple::vector(1.0, 0.0, 0.0);

        let tri = Shape::new_smooth_triangle(p1, p2, p3, n1, n2, n3);
        let r = Ray::new(Tuple::point(-0.2, 0.3, -2.0), Tuple::vector(0.0, 0.0, 1.0));
        let intersections = tri.intersect(&r, Vector::new());
        let xs = intersections.0.iter().collect::<Vec<_>>();

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0].u, 0.45);
        assert_abs_diff_eq!(xs[0].v, 0.25);
    }

    #[test]
    fn normal_on_translated_sphere() {
        let shape = Shape::new_sphere().with_transform(Matrix::translation(0.0, 1.0, 0.0));
        let n = shape.normal_at(
            &Tuple::point(0.0, 1.70711, -FRAC_1_SQRT_2),
            &Intersection::new(0.0, &shape, Vector::new()),
            &Vector::new(),
        );

        assert_abs_diff_eq!(n, Tuple::vector(0.0, FRAC_1_SQRT_2, -FRAC_1_SQRT_2));
    }

    #[test]
    fn normal_on_transformed_sphere() {
        let shape = Shape::new_sphere()
            .with_transform(Matrix::scaling(1.0, 0.5, 1.0) * Matrix::rotation_z(FRAC_PI_4));
        let n = shape.normal_at(
            &Tuple::point(0.0, 2.0_f64.sqrt() / 2.0, -(2.0_f64.sqrt()) / 2.0),
            &Intersection::new(0.0, &shape, Vector::new()),
            &Vector::new(),
        );

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 0.97014, -0.24254));
    }

    #[test]
    fn sphere_normal_x_axis() {
        let s = SingleKind::Sphere;
        let n = s.normal_at(&Tuple::point(1.0, 0.0, 0.0));

        assert_abs_diff_eq!(n, Tuple::vector(1.0, 0.0, 0.0));
    }

    #[test]
    fn sphere_normal_y_axis() {
        let s = SingleKind::Sphere;
        let n = s.normal_at(&Tuple::point(0.0, 1.0, 0.0));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 1.0, 0.0));
    }

    #[test]
    fn sphere_normal_z_axis() {
        let s = SingleKind::Sphere;
        let n = s.normal_at(&Tuple::point(0.0, 0.0, 1.0));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 0.0, 1.0));
    }

    #[test]
    fn sphere_normal_non_axial_point() {
        let s = SingleKind::Sphere;
        let n = s.normal_at(&Tuple::point(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
        ));

        assert_abs_diff_eq!(
            n,
            Tuple::vector(
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0
            )
        );
    }

    #[test]
    fn normal_is_normalized() {
        let s = SingleKind::Sphere;
        let n = s.normal_at(&Tuple::point(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
        ));

        assert_abs_diff_eq!(n, n.normalize());
    }

    #[test]
    fn normal_of_plane_is_constant() {
        let plane = SingleKind::Plane;
        let n1 = plane.normal_at(&Tuple::point(0.0, 0.0, 0.0));
        let n2 = plane.normal_at(&Tuple::point(10.0, 0.0, -10.0));
        let n3 = plane.normal_at(&Tuple::point(-5.0, 0.0, 150.0));

        assert_abs_diff_eq!(n1, Tuple::vector(0.0, 1.0, 0.0));
        assert_abs_diff_eq!(n2, Tuple::vector(0.0, 1.0, 0.0));
        assert_abs_diff_eq!(n3, Tuple::vector(0.0, 1.0, 0.0));
    }

    #[test]
    fn normal_on_surface_of_cube() {
        let scenarios = [
            (Tuple::point(1.0, 0.5, -0.8), Tuple::vector(1.0, 0.0, 0.0)),
            (Tuple::point(-1.0, -0.2, 0.9), Tuple::vector(-1.0, 0.0, 0.0)),
            (Tuple::point(-0.4, 1.0, -0.1), Tuple::vector(0.0, 1.0, 0.0)),
            (Tuple::point(0.3, -1.0, -0.7), Tuple::vector(0.0, -1.0, 0.0)),
            (Tuple::point(-0.6, 0.3, 1.0), Tuple::vector(0.0, 0.0, 1.0)),
            (Tuple::point(0.4, 0.4, -1.0), Tuple::vector(0.0, 0.0, -1.0)),
            (Tuple::point(1.0, 1.0, 1.0), Tuple::vector(1.0, 0.0, 0.0)),
            (
                Tuple::point(-1.0, -1.0, -1.0),
                Tuple::vector(-1.0, 0.0, 0.0),
            ),
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
            (Tuple::point(1.0, 0.0, 0.0), Tuple::vector(1.0, 0.0, 0.0)),
            (Tuple::point(0.0, 5.0, -1.0), Tuple::vector(0.0, 0.0, -1.0)),
            (Tuple::point(0.0, -2.0, 1.0), Tuple::vector(0.0, 0.0, 1.0)),
            (Tuple::point(-1.0, 1.0, 0.0), Tuple::vector(-1.0, 0.0, 0.0)),
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
            (Tuple::point(0.0, 1.0, 0.03), Tuple::vector(0.0, -1.0, 0.0)),
            (Tuple::point(0.5, 1.0, 0.0), Tuple::vector(0.0, -1.0, 0.0)),
            (Tuple::point(0.0, 1.0, 0.5), Tuple::vector(0.0, -1.0, 0.0)),
            (Tuple::point(0.0, 2.0, 0.0), Tuple::vector(0.0, 1.0, 0.0)),
            (Tuple::point(0.5, 2.0, 0.0), Tuple::vector(0.0, 1.0, 0.0)),
            (Tuple::point(0.0, 2.0, 0.5), Tuple::vector(0.0, 1.0, 0.0)),
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
            (Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 0.0)),
            (
                Tuple::point(1.0, 1.0, 1.0),
                Tuple::vector(1.0, -(2.0_f64.sqrt()), 1.0),
            ),
            (Tuple::point(-1.0, -1.0, 0.0), Tuple::vector(-1.0, 1.0, 0.0)),
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
            &Tuple::point(1.7321, 1.1547, -5.5774),
            &Intersection::new(0.0, &s, Vector::new()),
            &vector![g2.transform, g1.transform],
        );

        assert_abs_diff_eq!(n, Tuple::vector(0.2857, 0.4286, -0.8571));
    }

    #[test]
    fn normal_on_triangle() {
        let t = Shape::new_triangle(
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::point(-1.0, 0.0, 0.0),
            Tuple::point(1.0, 0.0, 0.0),
        );

        assert_eq!(
            t.normal_at(
                &Tuple::point(0.0, 0.5, 0.0),
                &Intersection::new(0.0, &t, Vector::new()),
                &Vector::new()
            ),
            Tuple::vector(0.0, 0.0, 1.0)
        );
        assert_eq!(
            t.normal_at(
                &Tuple::point(-0.5, 0.75, 0.0),
                &Intersection::new(0.0, &t, Vector::new()),
                &Vector::new()
            ),
            Tuple::vector(0.0, 0.0, 1.0)
        );
        assert_eq!(
            t.normal_at(
                &Tuple::point(0.5, 0.25, 0.0),
                &Intersection::new(0.0, &t, Vector::new()),
                &Vector::new()
            ),
            Tuple::vector(0.0, 0.0, 1.0)
        );
    }

    #[test]
    fn smooth_triangle_uses_uv_to_interpolate_normal() {
        let p1 = Tuple::point(0.0, 1.0, 0.0);
        let p2 = Tuple::point(-1.0, 0.0, 0.0);
        let p3 = Tuple::point(1.0, 0.0, 0.0);
        let n1 = Tuple::vector(0.0, 1.0, 0.0);
        let n2 = Tuple::vector(-1.0, 0.0, 0.0);
        let n3 = Tuple::vector(1.0, 0.0, 0.0);

        let tri = Shape::new_smooth_triangle(p1, p2, p3, n1, n2, n3);
        let i = Intersection::new_with_uv(1.0, &tri, Vector::new(), 0.45, 0.25);
        let n = tri.normal_at(&Tuple::point(0.0, 0.0, 0.0), &i, &Vector::new());

        assert_eq!(n, Tuple::vector(-0.5547, 0.83205, 0.0));
    }

    #[test]
    fn convert_point_world_to_object() {
        let s = Shape::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));
        let g2 = Shape::new_group(vec![s.clone()]).with_transform(Matrix::scaling(2.0, 2.0, 2.0));
        let g1 = Shape::new_group(vec![g2.clone()]).with_transform(Matrix::rotation_y(FRAC_PI_2));

        assert_eq!(
            s.world_to_object(
                &Tuple::point(-2.0, 0.0, -10.0),
                &vector![g2.transform, g1.transform]
            ),
            Tuple::point(0.0, 0.0, -1.0)
        );
    }

    #[test]
    fn convert_normal_object_to_world() {
        let s = Shape::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));
        let g2 = Shape::new_group(vec![s.clone()]).with_transform(Matrix::scaling(1.0, 2.0, 3.0));
        let g1 = Shape::new_group(vec![g2.clone()]).with_transform(Matrix::rotation_y(FRAC_PI_2));

        let n = s.normal_to_world(
            &Tuple::vector(
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0,
            ),
            &vector![g2.transform, g1.transform],
        );

        assert_eq!(n, Tuple::vector(0.2857, 0.4286, -0.8571));
    }
}
