use std::ops;

/// A tuple using left-hand coordinates where:
/// * `w = 0.0` is a vector
/// * `w = 1.0` is a point
#[derive(Debug, PartialEq)]
pub struct Tuple {
    x: f64,
    y: f64,
    z: f64,
    w: f64,
}

impl Tuple {
    #[must_use]
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Self {
        Self { x, y, z, w }
    }

    pub fn point(x: f64, y: f64, z: f64) -> Self {
        Self::new(x, y, z, 1.0)
    }

    pub fn vector(x: f64, y: f64, z: f64) -> Self {
        Self::new(x, y, z, 0.0)
    }
}

impl ops::Add for Tuple {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self::new(
            self.x + other.x,
            self.y + other.y,
            self.z + other.z,
            self.w + other.w,
        )
    }
}

impl ops::Neg for Tuple {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z, -self.w)
    }
}

impl ops::Sub for Tuple {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self + -other
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tuple_point() {
        let p = Tuple::new(4.3, -4.2, 3.1, 1.0);

        assert_eq!(p.x, 4.3);
        assert_eq!(p.y, -4.2);
        assert_eq!(p.z, 3.1);
        assert_eq!(p.w, 1.0);
    }

    #[test]
    fn test_tuple_vector() {
        let v = Tuple::new(4.3, -4.2, 3.1, 0.0);

        assert_eq!(v.x, 4.3);
        assert_eq!(v.y, -4.2);
        assert_eq!(v.z, 3.1);
        assert_eq!(v.w, 0.0);
    }

    #[test]
    fn test_tuple_point_equals_point() {
        let p1 = Tuple::point(4.3, -4.2, 3.1);
        let p2 = Tuple::new(4.3, -4.2, 3.1, 1.0);

        assert_eq!(p1, p2);
    }

    #[test]
    fn test_tuple_vector_equals_vector() {
        let v1 = Tuple::vector(4.3, -4.2, 3.1);
        let v2 = Tuple::new(4.3, -4.2, 3.1, 0.0);

        assert_eq!(v1, v2);
    }

    #[test]
    fn test_tuple_add() {
        let a1 = Tuple::new(3.0, -2.0, 5.0, 1.0);
        let a2 = Tuple::new(-2.0, 3.0, 1.0, 0.0);

        let result = a1 + a2;

        assert_eq!(result, Tuple::new(1.0, 1.0, 6.0, 1.0));
    }

    #[test]
    fn test_tuple_sub_two_points() {
        let p1 = Tuple::point(3.0, 2.0, 1.0);
        let p2 = Tuple::point(5.0, 6.0, 7.0);

        let result = p1 - p2;

        assert_eq!(result, Tuple::vector(-2.0, -4.0, -6.0));
    }

    #[test]
    fn test_tuple_sub_vector_from_point() {
        let p = Tuple::point(3.0, 2.0, 1.0);
        let v = Tuple::vector(5.0, 6.0, 7.0);

        let result = p - v;

        assert_eq!(result, Tuple::point(-2.0, -4.0, -6.0));
    }

    #[test]
    fn test_tuple_sub_two_vectors() {
        let v1 = Tuple::vector(3.0, 2.0, 1.0);
        let v2 = Tuple::vector(5.0, 6.0, 7.0);

        let result = v1 - v2;

        assert_eq!(result, Tuple::vector(-2.0, -4.0, -6.0));
    }

    #[test]
    fn test_sub_vector_zero_vector() {
        let zero = Tuple::vector(0.0, 0.0, 0.0);
        let v = Tuple::vector(1.0, -2.0, 3.0);

        let result = zero - v;

        assert_eq!(result, Tuple::vector(-1.0, 2.0, -3.0));
    }

    #[test]
    fn test_tuple_neg() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);

        let result = -a;

        assert_eq!(result, Tuple::new(-1.0, 2.0, -3.0, 4.0));
    }
}
