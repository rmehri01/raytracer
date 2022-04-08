use std::collections::BTreeSet;

use approx::AbsDiffEq;
use ordered_float::OrderedFloat;

use crate::core::tuple::Tuple;

use super::sphere::Sphere;

// TODO: what to do about partial eq and eq
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Intersection {
    pub t: OrderedFloat<f64>,
    pub object: Sphere,
}

impl Intersection {
    pub fn new(t: f64, object: Sphere) -> Self {
        Self {
            t: OrderedFloat(t),
            object,
        }
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

impl Eq for Intersection {}

impl PartialOrd for Intersection {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.t.partial_cmp(&other.t)
    }
}

impl Ord for Intersection {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.t.cmp(&other.t)
    }
}

/// Encapsulates precomputed information for an intersection.
pub struct Computations {
    pub t: f64,
    pub object: Sphere,
    pub point: Tuple,
    pub eyev: Tuple,
    pub normalv: Tuple,
    pub inside: bool,
}

// TODO: need to keep duplicates?
pub struct Intersections(pub BTreeSet<Intersection>);

impl Intersections {
    pub fn hit(&self) -> Option<Intersection> {
        self.0.iter().find(|i| i.t >= OrderedFloat(0.0)).cloned()
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

        assert_relative_eq!(intersection.t.0, 3.5);
        assert_eq!(intersection.object, s);
    }

    #[test]
    fn hit_all_positive() {
        let s = Sphere::default();
        let i1 = Intersection::new(1.0, s);
        let i2 = Intersection::new(2.0, s);
        let xs = BTreeSet::from([i1, i2]);

        let hit = Intersections(xs).hit().expect("valid hit");

        assert_abs_diff_eq!(hit, i1);
    }

    #[test]
    fn hit_some_negative() {
        let s = Sphere::default();
        let i1 = Intersection::new(-1.0, s);
        let i2 = Intersection::new(1.0, s);
        let xs = BTreeSet::from([i1, i2]);

        let hit = Intersections(xs).hit().expect("valid hit");

        assert_abs_diff_eq!(hit, i2);
    }

    #[test]
    fn hit_all_negative() {
        let s = Sphere::default();
        let i1 = Intersection::new(-2.0, s);
        let i2 = Intersection::new(-1.0, s);
        let xs = BTreeSet::from([i1, i2]);

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
        let xs = BTreeSet::from([i1, i2, i3, i4]);

        let hit = Intersections(xs).hit().expect("valid hit");

        assert_abs_diff_eq!(hit, i4);
    }
}
