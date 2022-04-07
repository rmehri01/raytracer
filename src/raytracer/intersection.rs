use std::ops;

use approx::AbsDiffEq;

use super::sphere::Sphere;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Intersection {
    pub t: f64,
    pub object: Sphere,
}

impl Intersection {
    pub fn new(t: f64, object: Sphere) -> Self {
        Self { t, object }
    }
}

impl AbsDiffEq for Intersection {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-5
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.t.abs_diff_eq(&other.t, epsilon) && self.object.abs_diff_eq(&other.object, epsilon)
    }
}

// TODO: maybe sorted list instead of vector
pub struct Intersections(pub Vec<Intersection>);

impl Intersections {
    pub fn hit(&self) -> Option<Intersection> {
        self.0
            .iter()
            .filter(|i| i.t >= 0.0)
            .min_by(|i1, i2| i1.t.partial_cmp(&i2.t).unwrap())
            .cloned()
    }
}

impl ops::Index<usize> for Intersections {
    type Output = Intersection;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

#[cfg(test)]
mod tests {
    use approx::{assert_abs_diff_eq, assert_relative_eq};

    use crate::raytracer::sphere::Sphere;

    use super::*;

    #[test]
    fn test_intersection_new() {
        let s = Sphere::default();
        let intersection = Intersection::new(3.5, s);

        assert_relative_eq!(intersection.t, 3.5);
        assert_eq!(intersection.object, s);
    }

    #[test]
    fn hit_all_positive() {
        let s = Sphere::default();
        let i1 = Intersection::new(1.0, s);
        let i2 = Intersection::new(2.0, s);
        let xs = vec![i1, i2];

        let hit = Intersections(xs).hit().expect("valid hit");

        assert_abs_diff_eq!(hit, i1);
    }

    #[test]
    fn hit_some_negative() {
        let s = Sphere::default();
        let i1 = Intersection::new(-1.0, s);
        let i2 = Intersection::new(1.0, s);
        let xs = vec![i1, i2];

        let hit = Intersections(xs).hit().expect("valid hit");

        assert_abs_diff_eq!(hit, i2);
    }

    #[test]
    fn hit_all_negative() {
        let s = Sphere::default();
        let i1 = Intersection::new(-2.0, s);
        let i2 = Intersection::new(-1.0, s);
        let xs = vec![i1, i2];

        let hit = Intersections(xs).hit();

        assert_eq!(hit, None);
    }

    #[test]
    fn hit_lowest_nonnegative() {
        let s = Sphere::default();
        let i1 = Intersection::new(5.0, s);
        let i2 = Intersection::new(7.0, s);
        let i3 = Intersection::new(-3.0, s);
        let i4 = Intersection::new(2.0, s);
        let xs = vec![i1, i2, i3, i4];

        let hit = Intersections(xs).hit().expect("valid hit");

        assert_abs_diff_eq!(hit, i4);
    }
}
