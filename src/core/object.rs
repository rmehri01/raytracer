use approx::AbsDiffEq;

use super::{intersection::Intersection, ray::Ray, tuple::Tuple};

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
pub struct Sphere {}

impl Sphere {
    pub fn new() -> Self {
        Self {}
    }

    fn intersects(&self, ray: &Ray) -> Option<(Intersection, Intersection)> {
        let sphere_to_ray = ray.origin - Tuple::point(0.0, 0.0, 0.0);

        let a = ray.direction.dot(&ray.direction);
        let b = 2.0 * sphere_to_ray.dot(&ray.direction);
        let c = sphere_to_ray.dot(&sphere_to_ray) - 1.0;

        let discriminant = b * b - 4.0 * a * c;

        if discriminant < 0.0 {
            None
        } else {
            let t1 = (-b - discriminant.sqrt()) / (2.0 * a);
            let t2 = (-b + discriminant.sqrt()) / (2.0 * a);

            let i1 = Intersection::new(t1, Object::Sphere(*self));
            let i2 = Intersection::new(t2, Object::Sphere(*self));
            Some((i1, i2))
        }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn ray_intersects_sphere() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let (i1, i2) = s.intersects(&r).expect("two valid intersections");

        assert_abs_diff_eq!(i1.t, 4.0);
        assert_abs_diff_eq!(i2.t, 6.0);
    }

    #[test]
    fn ray_tangent_to_sphere() {
        let r = Ray::new(Tuple::point(0.0, 1.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let (i1, i2) = s.intersects(&r).expect("two valid intersections");

        assert_abs_diff_eq!(i1.t, 5.0);
        assert_abs_diff_eq!(i2.t, 5.0);
    }

    #[test]
    fn ray_misses_sphere() {
        let r = Ray::new(Tuple::point(0.0, 2.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let intersects = s.intersects(&r);

        assert!(intersects.is_none());
    }

    #[test]
    fn ray_originates_inside_sphere() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let (i1, i2) = s.intersects(&r).expect("two valid intersections");

        assert_abs_diff_eq!(i1.t, -1.0);
        assert_abs_diff_eq!(i2.t, 1.0);
    }

    #[test]
    fn sphere_behind_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let (i1, i2) = s.intersects(&r).expect("two valid intersections");

        assert_abs_diff_eq!(i1.t, -6.0);
        assert_abs_diff_eq!(i2.t, -4.0);
    }

    #[test]
    fn intersect_sets_object() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = Sphere::new();

        let (i1, i2) = s.intersects(&r).expect("two valid intersections");

        assert_eq!(i1.object, Object::Sphere(s));
        assert_eq!(i2.object, Object::Sphere(s));
    }
}
