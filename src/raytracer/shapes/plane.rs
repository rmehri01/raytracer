use approx::AbsDiffEq;

use crate::{core::tuple::Tuple, raytracer::ray::Ray};

/// A perfectly flat surface that extends infinitely in x and z.
#[derive(Debug, PartialEq, Clone, Copy, Default)]
pub struct Plane {}

impl Plane {
    pub fn intersect(&self, ray: &Ray) -> Vec<f64> {
        let epsilon = <Tuple as AbsDiffEq>::default_epsilon();
        if ray.direction.y.abs() < epsilon {
            Vec::new()
        } else {
            vec![-ray.origin.y / ray.direction.y]
        }
    }

    pub fn normal_at(&self, _p: &Tuple) -> Tuple {
        Tuple::vector(0.0, 1.0, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn normal_of_plane_is_constant() {
        let plane = Plane::default();
        let n1 = plane.normal_at(&Tuple::point(0.0, 0.0, 0.0));
        let n2 = plane.normal_at(&Tuple::point(10.0, 0.0, -10.0));
        let n3 = plane.normal_at(&Tuple::point(-5.0, 0.0, 150.0));

        assert_abs_diff_eq!(n1, Tuple::vector(0.0, 1.0, 0.0));
        assert_abs_diff_eq!(n2, Tuple::vector(0.0, 1.0, 0.0));
        assert_abs_diff_eq!(n3, Tuple::vector(0.0, 1.0, 0.0));
    }

    #[test]
    fn intersecting_a_ray_parallel_to_the_plane() {
        let plane = Plane::default();
        let r = Ray::new(Tuple::point(0.0, 10.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = plane.intersect(&r);

        assert!(xs.is_empty());
    }

    #[test]
    fn intersecting_a_coplanar_ray() {
        let plane = Plane::default();
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = plane.intersect(&r);

        assert!(xs.is_empty());
    }

    #[test]
    fn ray_intersecting_from_above() {
        let plane = Plane::default();
        let r = Ray::new(Tuple::point(0.0, 1.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));

        let xs = plane.intersect(&r);

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 1.0);
    }

    #[test]
    fn ray_intersecting_from_below() {
        let plane = Plane::default();
        let r = Ray::new(Tuple::point(0.0, -1.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));

        let xs = plane.intersect(&r);

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 1.0);
    }
}
