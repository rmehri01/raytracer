use crate::{core::tuple::Tuple, raytracer::ray::Ray};

/// A unit sphere with center at the origin.
#[derive(Debug, PartialEq, Clone, Copy, Default)]
pub struct Sphere {}

impl Sphere {
    pub fn intersect(&self, ray: &Ray) -> Vec<f64> {
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

    pub fn normal_at(&self, object_point: &Tuple) -> Tuple {
        *object_point - Tuple::point(0.0, 0.0, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn ray_intersects_sphere() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::default();

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], 4.0);
        assert_abs_diff_eq!(xs[1], 6.0);
    }

    #[test]
    fn ray_tangent_to_sphere() {
        let r = Ray::new(Tuple::point(0.0, 1.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::default();

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], 5.0);
        assert_abs_diff_eq!(xs[1], 5.0);
    }

    #[test]
    fn ray_misses_sphere() {
        let r = Ray::new(Tuple::point(0.0, 2.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::default();

        let intersects = s.intersect(&r);

        assert!(intersects.is_empty());
    }

    #[test]
    fn ray_originates_inside_sphere() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::default();

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], -1.0);
        assert_abs_diff_eq!(xs[1], 1.0);
    }

    #[test]
    fn sphere_behind_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::default();

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], -6.0);
        assert_abs_diff_eq!(xs[1], -4.0);
    }

    #[test]
    fn sphere_normal_x_axis() {
        let s = Sphere::default();
        let n = s.normal_at(&Tuple::point(1.0, 0.0, 0.0));

        assert_abs_diff_eq!(n, Tuple::vector(1.0, 0.0, 0.0));
    }

    #[test]
    fn sphere_normal_y_axis() {
        let s = Sphere::default();
        let n = s.normal_at(&Tuple::point(0.0, 1.0, 0.0));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 1.0, 0.0));
    }

    #[test]
    fn sphere_normal_z_axis() {
        let s = Sphere::default();
        let n = s.normal_at(&Tuple::point(0.0, 0.0, 1.0));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 0.0, 1.0));
    }

    #[test]
    fn sphere_normal_non_axial_point() {
        let s = Sphere::default();
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
        let s = Sphere::default();
        let n = s.normal_at(&Tuple::point(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
        ));

        assert_abs_diff_eq!(n, n.normalize());
    }
}
