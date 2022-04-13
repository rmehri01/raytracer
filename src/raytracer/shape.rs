use std::collections::BTreeSet;

use approx::AbsDiffEq;

use crate::core::{matrix::Matrix, tuple::Tuple};

use super::{
    intersection::{Intersection, Intersections},
    material::Material,
    ray::Ray,
};

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Shape {
    // TODO: better way of setting/using this
    pub transform: Matrix<4>,
    pub material: Material,
    kind: ShapeKind,
}

impl Shape {
    pub fn new(kind: ShapeKind) -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
            kind,
        }
    }

    pub fn new_sphere() -> Self {
        Self::new(ShapeKind::Sphere)
    }

    pub fn new_glass_sphere() -> Self {
        let mut sphere = Self::new_sphere();
        sphere.material.transparency = 1.0;
        sphere.material.refractive_index = 1.5;

        sphere
    }

    pub fn new_plane() -> Self {
        Self::new(ShapeKind::Plane)
    }

    pub fn new_cube() -> Self {
        Self::new(ShapeKind::Cube)
    }

    pub fn new_cylinder() -> Self {
        Self::new(ShapeKind::Cylinder(Conic::default()))
    }

    pub fn new_cone() -> Self {
        Self::new(ShapeKind::Cone(Conic::default()))
    }

    pub fn intersect(&self, ray: &Ray) -> Intersections {
        let local_ray = ray.transform(&self.transform.inverse());
        let ts = self.kind.intersect(&local_ray);
        let intersections = ts.iter().map(|t| Intersection::new(*t, *self));

        Intersections(BTreeSet::from_iter(intersections))
    }

    pub fn normal_at(&self, point: &Tuple) -> Tuple {
        let local_point = self.transform.inverse() * *point;
        let local_normal = self.kind.normal_at(&local_point);

        let mut world_normal = self.transform.inverse().transpose() * local_normal;
        world_normal.w = 0.0;

        world_normal.normalize()
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum ShapeKind {
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
}

impl ShapeKind {
    pub fn intersect(&self, ray: &Ray) -> Vec<f64> {
        match self {
            Self::Sphere => Self::sphere_intersect(ray),
            Self::Plane => Self::plane_intersect(ray),
            Self::Cube => Self::cube_intersect(ray),
            Self::Cylinder(conic) => Self::cylinder_intersect(ray, conic),
            Self::Cone(conic) => Self::cone_intersect(ray, conic),
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
}

#[derive(Debug, Clone, Copy)]
pub struct Conic {
    minimum: f64,
    maximum: f64,
    closed: bool,
}

impl Conic {
    fn new(minimum: f64, maximum: f64, closed: bool) -> Self {
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
    use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_4};

    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn default_transformation() {
        let shape = Shape::new_sphere();

        assert_eq!(shape.transform, Matrix::identity());
    }

    #[test]
    fn change_transformation() {
        let mut shape = Shape::new_sphere();
        shape.transform = Matrix::translation(2.0, 3.0, 4.0);

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

        let mut shape = Shape::new_sphere();
        shape.material = material;

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
    fn intersect_scaled_sphere_with_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut shape = Shape::new_sphere();
        shape.transform = Matrix::scaling(2.0, 2.0, 2.0);

        let xs = shape
            .intersect(&r)
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
        let mut shape = Shape::new_sphere();
        shape.transform = Matrix::translation(5.0, 0.0, 0.0);

        assert!(shape.intersect(&r).0.is_empty());
    }

    #[test]
    fn ray_intersects_sphere() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = ShapeKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], 4.0);
        assert_abs_diff_eq!(xs[1], 6.0);
    }

    #[test]
    fn ray_tangent_to_sphere() {
        let r = Ray::new(Tuple::point(0.0, 1.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = ShapeKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], 5.0);
        assert_abs_diff_eq!(xs[1], 5.0);
    }

    #[test]
    fn ray_misses_sphere() {
        let r = Ray::new(Tuple::point(0.0, 2.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = ShapeKind::Sphere;

        let intersects = s.intersect(&r);

        assert!(intersects.is_empty());
    }

    #[test]
    fn ray_originates_inside_sphere() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = ShapeKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], -1.0);
        assert_abs_diff_eq!(xs[1], 1.0);
    }

    #[test]
    fn sphere_behind_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = ShapeKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], -6.0);
        assert_abs_diff_eq!(xs[1], -4.0);
    }

    #[test]
    fn intersecting_a_ray_parallel_to_the_plane() {
        let plane = ShapeKind::Plane;
        let r = Ray::new(Tuple::point(0.0, 10.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = plane.intersect(&r);

        assert!(xs.is_empty());
    }

    #[test]
    fn intersecting_a_coplanar_ray() {
        let plane = ShapeKind::Plane;
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = plane.intersect(&r);

        assert!(xs.is_empty());
    }

    #[test]
    fn ray_intersecting_from_above() {
        let plane = ShapeKind::Plane;
        let r = Ray::new(Tuple::point(0.0, 1.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));

        let xs = plane.intersect(&r);

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 1.0);
    }

    #[test]
    fn ray_intersecting_from_below() {
        let plane = ShapeKind::Plane;
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
            .intersect(&r)
            .0
            .iter()
            .map(|i| i.shape)
            .collect::<Vec<_>>();

        assert_eq!(shapes[0], shape);
        assert_eq!(shapes[1], shape);
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
            let c = ShapeKind::Cube;
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
            let c = ShapeKind::Cube;
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
            let c = ShapeKind::Cylinder(Conic::default());
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
            let c = ShapeKind::Cylinder(Conic::default());
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
            let cyl = ShapeKind::Cylinder(Conic::new(1.0, 2.0, false));

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
            let cyl = ShapeKind::Cylinder(Conic::new(1.0, 2.0, true));

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
            let shape = ShapeKind::Cone(Conic::default());

            let r = Ray::new(origin, direction.normalize());
            let xs = shape.intersect(&r);

            assert_eq!(xs.len(), 2);
            assert_abs_diff_eq!(xs[0], t0, epsilon = Tuple::default_epsilon());
            assert_abs_diff_eq!(xs[1], t1, epsilon = Tuple::default_epsilon());
        }
    }

    #[test]
    fn intersect_cone_with_ray_parallel_to_one_half() {
        let shape = ShapeKind::Cone(Conic::default());

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
            let shape = ShapeKind::Cone(Conic::new(-0.5, 0.5, true));

            let r = Ray::new(origin, direction.normalize());
            let xs = shape.intersect(&r);

            assert_eq!(xs.len(), count);
        }
    }

    #[test]
    fn normal_on_translated_sphere() {
        let mut shape = Shape::new_sphere();
        shape.transform = Matrix::translation(0.0, 1.0, 0.0);

        let n = shape.normal_at(&Tuple::point(0.0, 1.70711, -FRAC_1_SQRT_2));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, FRAC_1_SQRT_2, -FRAC_1_SQRT_2));
    }

    #[test]
    fn normal_on_transformed_sphere() {
        let mut shape = Shape::new_sphere();
        shape.transform = Matrix::scaling(1.0, 0.5, 1.0) * Matrix::rotation_z(FRAC_PI_4);

        let n = shape.normal_at(&Tuple::point(
            0.0,
            2.0_f64.sqrt() / 2.0,
            -(2.0_f64.sqrt()) / 2.0,
        ));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 0.97014, -0.24254));
    }

    #[test]
    fn sphere_normal_x_axis() {
        let s = ShapeKind::Sphere;
        let n = s.normal_at(&Tuple::point(1.0, 0.0, 0.0));

        assert_abs_diff_eq!(n, Tuple::vector(1.0, 0.0, 0.0));
    }

    #[test]
    fn sphere_normal_y_axis() {
        let s = ShapeKind::Sphere;
        let n = s.normal_at(&Tuple::point(0.0, 1.0, 0.0));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 1.0, 0.0));
    }

    #[test]
    fn sphere_normal_z_axis() {
        let s = ShapeKind::Sphere;
        let n = s.normal_at(&Tuple::point(0.0, 0.0, 1.0));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 0.0, 1.0));
    }

    #[test]
    fn sphere_normal_non_axial_point() {
        let s = ShapeKind::Sphere;
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
        let s = ShapeKind::Sphere;
        let n = s.normal_at(&Tuple::point(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
        ));

        assert_abs_diff_eq!(n, n.normalize());
    }

    #[test]
    fn normal_of_plane_is_constant() {
        let plane = ShapeKind::Plane;
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
            let c = ShapeKind::Cube;
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
            let c = ShapeKind::Cylinder(Conic::default());
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
            let c = ShapeKind::Cylinder(Conic::new(1.0, 2.0, true));
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
            let c = ShapeKind::Cone(Conic::default());
            let n = c.normal_at(&point);

            assert_abs_diff_eq!(n, normal);
        }
    }
}
