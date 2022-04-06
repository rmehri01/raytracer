use approx::AbsDiffEq;

use crate::core::{matrix::Matrix, tuple::Tuple};

use super::{
    intersection::{Intersection, Intersections},
    ray::Ray,
};

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Object {
    Sphere(Sphere),
}

impl AbsDiffEq for Object {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-5
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        match (self, other) {
            (Object::Sphere(s1), Object::Sphere(s2)) => true,
        }
    }
}

/// A unit sphere with center at the origin.
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Sphere {
    transform: Matrix<4>,
}

impl Sphere {
    // TODO: should this be passed a transform
    pub fn new() -> Self {
        Self {
            transform: Matrix::identity(),
        }
    }

    pub fn intersect(&self, ray: &Ray) -> Intersections {
        let transformed_ray = ray.transform(&self.transform.inverse());
        let sphere_to_ray = transformed_ray.origin - Tuple::point(0.0, 0.0, 0.0);

        let a = transformed_ray.direction.dot(&transformed_ray.direction);
        let b = 2.0 * sphere_to_ray.dot(&transformed_ray.direction);
        let c = sphere_to_ray.dot(&sphere_to_ray) - 1.0;

        let discriminant = b * b - 4.0 * a * c;

        if discriminant < 0.0 {
            Intersections(Vec::new())
        } else {
            let t1 = (-b - discriminant.sqrt()) / (2.0 * a);
            let t2 = (-b + discriminant.sqrt()) / (2.0 * a);

            let i1 = Intersection::new(t1, Object::Sphere(*self));
            let i2 = Intersection::new(t2, Object::Sphere(*self));

            Intersections(vec![i1, i2])
        }
    }

    pub fn normal_at(&self, world_point: &Tuple) -> Tuple {
        let object_point = self.transform.inverse() * *world_point;
        let object_normal = object_point - Tuple::point(0.0, 0.0, 0.0);
        let mut world_normal = self.transform.inverse().transpose() * object_normal;
        world_normal.w = 0.0;

        world_normal.normalize()
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_4};

    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn ray_intersects_sphere() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let xs = s.intersect(&r);

        assert_eq!(xs.0.len(), 2);
        assert_abs_diff_eq!(xs[0].t, 4.0);
        assert_abs_diff_eq!(xs[1].t, 6.0);
    }

    #[test]
    fn ray_tangent_to_sphere() {
        let r = Ray::new(Tuple::point(0.0, 1.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let xs = s.intersect(&r);

        assert_eq!(xs.0.len(), 2);
        assert_abs_diff_eq!(xs[0].t, 5.0);
        assert_abs_diff_eq!(xs[1].t, 5.0);
    }

    #[test]
    fn ray_misses_sphere() {
        let r = Ray::new(Tuple::point(0.0, 2.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let intersects = s.intersect(&r);

        assert!(intersects.0.is_empty());
    }

    #[test]
    fn ray_originates_inside_sphere() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let xs = s.intersect(&r);

        assert_eq!(xs.0.len(), 2);
        assert_abs_diff_eq!(xs[0].t, -1.0);
        assert_abs_diff_eq!(xs[1].t, 1.0);
    }

    #[test]
    fn sphere_behind_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let xs = s.intersect(&r);

        assert_eq!(xs.0.len(), 2);
        assert_abs_diff_eq!(xs[0].t, -6.0);
        assert_abs_diff_eq!(xs[1].t, -4.0);
    }

    #[test]
    fn intersect_sets_object() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let xs = s.intersect(&r);

        assert_eq!(xs[0].object, Object::Sphere(s));
        assert_eq!(xs[1].object, Object::Sphere(s));
    }

    #[test]
    fn default_transformation() {
        let s = Sphere::new();

        assert_eq!(s.transform, Matrix::identity());
    }

    #[test]
    fn change_sphere_transformation() {
        let mut s = Sphere::new();
        s.transform = Matrix::translation(2.0, 3.0, 4.0);

        assert_eq!(s.transform, Matrix::translation(2.0, 3.0, 4.0));
    }

    #[test]
    fn intersect_scaled_sphere_with_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut s = Sphere::new();
        s.transform = Matrix::scaling(2.0, 2.0, 2.0);

        let xs = s.intersect(&r);

        assert_eq!(xs.0.len(), 2);
        assert_abs_diff_eq!(xs[0].t, 3.0);
        assert_abs_diff_eq!(xs[1].t, 7.0);
    }

    #[test]
    fn sphere_normal_x_axis() {
        let s = Sphere::new();
        let n = s.normal_at(&Tuple::point(1.0, 0.0, 0.0));

        assert_abs_diff_eq!(n, Tuple::vector(1.0, 0.0, 0.0));
    }

    #[test]
    fn sphere_normal_y_axis() {
        let s = Sphere::new();
        let n = s.normal_at(&Tuple::point(0.0, 1.0, 0.0));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 1.0, 0.0));
    }

    #[test]
    fn sphere_normal_z_axis() {
        let s = Sphere::new();
        let n = s.normal_at(&Tuple::point(0.0, 0.0, 1.0));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 0.0, 1.0));
    }

    #[test]
    fn sphere_normal_non_axial_point() {
        let s = Sphere::new();
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
        let s = Sphere::new();
        let n = s.normal_at(&Tuple::point(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
        ));

        assert_abs_diff_eq!(n, n.normalize());
    }

    #[test]
    fn normal_on_translated_sphere() {
        let mut s = Sphere::new();
        s.transform = Matrix::translation(0.0, 1.0, 0.0);

        let n = s.normal_at(&Tuple::point(0.0, 1.70711, -FRAC_1_SQRT_2));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, FRAC_1_SQRT_2, -FRAC_1_SQRT_2));
    }

    #[test]
    fn normal_on_transformed_sphere() {
        let mut s = Sphere::new();
        s.transform = Matrix::scaling(1.0, 0.5, 1.0) * Matrix::rotation_z(FRAC_PI_4);

        let n = s.normal_at(&Tuple::point(
            0.0,
            2.0_f64.sqrt() / 2.0,
            -(2.0_f64.sqrt()) / 2.0,
        ));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 0.97014, -0.24254));
    }
}