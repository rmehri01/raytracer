use std::collections::BTreeSet;

use approx::AbsDiffEq;
use ordered_float::OrderedFloat;

use super::object::Object;

// TODO: what to do about partial eq and eq
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Intersection {
    pub t: OrderedFloat<f64>,
    pub object: Object,
}

impl Intersection {
    pub fn new(t: f64, object: Object) -> Self {
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

// TODO: need to keep duplicates?
pub struct Intersections(pub BTreeSet<Intersection>);

impl Intersections {
    pub fn hit(&self) -> Option<Intersection> {
        self.0.iter().find(|i| i.t.0 >= 0.0).cloned()
    }
}

#[cfg(test)]
mod tests {
    use approx::{assert_abs_diff_eq, assert_relative_eq};

    use super::*;

    #[test]
    fn test_intersection_new() {
        let obj = Object::new_sphere();
        let intersection = Intersection::new(3.5, obj);

        assert_relative_eq!(intersection.t.0, 3.5);
        assert_eq!(intersection.object, obj);
    }

    #[test]
    fn hit_all_positive() {
        let obj = Object::new_sphere();
        let i1 = Intersection::new(1.0, obj);
        let i2 = Intersection::new(2.0, obj);
        let xs = BTreeSet::from([i1, i2]);

        let hit = Intersections(xs).hit().expect("valid hit");

        assert_abs_diff_eq!(hit, i1);
    }

    #[test]
    fn hit_some_negative() {
        let obj = Object::new_sphere();
        let i1 = Intersection::new(-1.0, obj);
        let i2 = Intersection::new(1.0, obj);
        let xs = BTreeSet::from([i1, i2]);

        let hit = Intersections(xs).hit().expect("valid hit");

        assert_abs_diff_eq!(hit, i2);
    }

    #[test]
    fn hit_all_negative() {
        let obj = Object::new_sphere();
        let i1 = Intersection::new(-2.0, obj);
        let i2 = Intersection::new(-1.0, obj);
        let xs = BTreeSet::from([i1, i2]);

        let hit = Intersections(xs).hit();

        assert_eq!(hit, None);
    }

    #[test]
    fn hit_lowest_nonnegative() {
        let obj = Object::new_sphere();
        let i1 = Intersection::new(5.0, obj);
        let i2 = Intersection::new(7.0, obj);
        let i3 = Intersection::new(-3.0, obj);
        let i4 = Intersection::new(2.0, obj);
        let xs = BTreeSet::from([i1, i2, i3, i4]);

        let hit = Intersections(xs).hit().expect("valid hit");

        assert_abs_diff_eq!(hit, i4);
    }
}
