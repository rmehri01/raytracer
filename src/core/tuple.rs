use std::ops;

use approx::AbsDiffEq;

/// A tuple using left-hand coordinates where:
/// * `w = 0.0` is a vector
/// * `w = 1.0` is a point
#[derive(Debug, Clone, Copy)]
pub struct Tuple {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
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

    pub fn magnitude(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2) + self.w.powi(2)).sqrt()
    }

    pub fn normalize(&self) -> Self {
        let m = self.magnitude();

        Self::new(self.x / m, self.y / m, self.z / m, self.w / m)
    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w
    }

    pub fn cross(&self, other: &Self) -> Self {
        Self::vector(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    pub fn reflect(&self, normal: &Self) -> Self {
        *self - *normal * 2.0 * self.dot(normal)
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

impl ops::Mul<f64> for Tuple {
    type Output = Self;

    fn mul(self, other: f64) -> Self {
        Self::new(
            self.x * other,
            self.y * other,
            self.z * other,
            self.w * other,
        )
    }
}

impl ops::Div<f64> for Tuple {
    type Output = Self;

    fn div(self, other: f64) -> Self {
        self * (1.0 / other)
    }
}

impl ops::Index<usize> for Tuple {
    type Output = f64;

    fn index(&self, index: usize) -> &f64 {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            3 => &self.w,
            _ => panic!("Index out of bounds"),
        }
    }
}

impl ops::IndexMut<usize> for Tuple {
    fn index_mut(&mut self, index: usize) -> &mut f64 {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            3 => &mut self.w,
            _ => panic!("Index out of bounds"),
        }
    }
}

impl PartialEq for Tuple {
    fn eq(&self, other: &Self) -> bool {
        self.abs_diff_eq(other, Self::default_epsilon())
    }
}

impl AbsDiffEq for Tuple {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-4
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.x.abs_diff_eq(&other.x, epsilon)
            && self.y.abs_diff_eq(&other.y, epsilon)
            && self.z.abs_diff_eq(&other.z, epsilon)
            && self.w.abs_diff_eq(&other.w, epsilon)
    }
}

#[cfg(test)]
mod tests {
    use approx::{assert_abs_diff_eq, assert_relative_eq};

    use super::*;

    #[test]
    fn tuple_point() {
        let p = Tuple::new(4.3, -4.2, 3.1, 1.0);

        assert_abs_diff_eq!(p, Tuple::point(4.3, -4.2, 3.1));
    }

    #[test]
    fn tuple_vector() {
        let v = Tuple::new(4.3, -4.2, 3.1, 0.0);

        assert_abs_diff_eq!(v, Tuple::vector(4.3, -4.2, 3.1));
    }

    #[test]
    fn tuple_point_equals_point() {
        let p1 = Tuple::point(4.3, -4.2, 3.1);
        let p2 = Tuple::new(4.3, -4.2, 3.1, 1.0);

        assert_abs_diff_eq!(p1, p2);
    }

    #[test]
    fn tuple_vector_equals_vector() {
        let v1 = Tuple::vector(4.3, -4.2, 3.1);
        let v2 = Tuple::new(4.3, -4.2, 3.1, 0.0);

        assert_abs_diff_eq!(v1, v2);
    }

    #[test]
    fn tuple_add() {
        let a1 = Tuple::new(3.0, -2.0, 5.0, 1.0);
        let a2 = Tuple::new(-2.0, 3.0, 1.0, 0.0);

        assert_abs_diff_eq!(a1 + a2, Tuple::new(1.0, 1.0, 6.0, 1.0));
    }

    #[test]
    fn tuple_sub_two_points() {
        let p1 = Tuple::point(3.0, 2.0, 1.0);
        let p2 = Tuple::point(5.0, 6.0, 7.0);

        assert_abs_diff_eq!(p1 - p2, Tuple::vector(-2.0, -4.0, -6.0));
    }

    #[test]
    fn tuple_sub_vector_from_point() {
        let p = Tuple::point(3.0, 2.0, 1.0);
        let v = Tuple::vector(5.0, 6.0, 7.0);

        assert_abs_diff_eq!(p - v, Tuple::point(-2.0, -4.0, -6.0));
    }

    #[test]
    fn tuple_sub_two_vectors() {
        let v1 = Tuple::vector(3.0, 2.0, 1.0);
        let v2 = Tuple::vector(5.0, 6.0, 7.0);

        assert_abs_diff_eq!(v1 - v2, Tuple::vector(-2.0, -4.0, -6.0));
    }

    #[test]
    fn sub_vector_zero_vector() {
        let zero = Tuple::vector(0.0, 0.0, 0.0);
        let v = Tuple::vector(1.0, -2.0, 3.0);

        assert_abs_diff_eq!(zero - v, Tuple::vector(-1.0, 2.0, -3.0));
    }

    #[test]
    fn tuple_neg() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);

        assert_abs_diff_eq!(-a, Tuple::new(-1.0, 2.0, -3.0, 4.0));
    }

    #[test]
    fn tuple_mul_scalar() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);

        assert_abs_diff_eq!(a * 3.5, Tuple::new(3.5, -7.0, 10.5, -14.0));
    }

    #[test]
    fn tuple_mul_fraction() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);

        assert_abs_diff_eq!(a * 0.5, Tuple::new(0.5, -1.0, 1.5, -2.0));
    }

    #[test]
    fn tuple_div_scalar() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);

        assert_abs_diff_eq!(a / 2.0, Tuple::new(0.5, -1.0, 1.5, -2.0));
    }

    #[test]
    fn magnitude_x() {
        let v = Tuple::vector(1.0, 0.0, 0.0);

        assert_relative_eq!(v.magnitude(), 1.0);
    }

    #[test]
    fn magnitude_y() {
        let v = Tuple::vector(0.0, 1.0, 0.0);

        assert_relative_eq!(v.magnitude(), 1.0);
    }

    #[test]
    fn magnitude_z() {
        let v = Tuple::vector(0.0, 0.0, 1.0);

        assert_relative_eq!(v.magnitude(), 1.0);
    }

    #[test]
    fn magnitude_all() {
        let v = Tuple::vector(1.0, 2.0, 3.0);

        assert_relative_eq!(v.magnitude(), 14.0_f64.sqrt());
    }

    #[test]
    fn magnitude_negative() {
        let v = Tuple::vector(-1.0, -2.0, -3.0);

        assert_relative_eq!(v.magnitude(), 14.0_f64.sqrt());
    }

    #[test]
    fn normalize_x() {
        let v = Tuple::vector(4.0, 0.0, 0.0);

        assert_abs_diff_eq!(v.normalize(), Tuple::vector(1.0, 0.0, 0.0));
    }

    #[test]
    fn normalize_all() {
        let v = Tuple::vector(1.0, 2.0, 3.0);

        assert_abs_diff_eq!(
            v.normalize(),
            Tuple::vector(
                1.0 / 14.0_f64.sqrt(),
                2.0 / 14.0_f64.sqrt(),
                3.0 / 14.0_f64.sqrt()
            )
        );
    }

    #[test]
    fn magnitude_of_normalized() {
        let v = Tuple::vector(1.0, 2.0, 3.0);

        assert_relative_eq!(v.normalize().magnitude(), 1.0);
    }

    #[test]
    fn dot_product() {
        let v1 = Tuple::vector(1.0, 2.0, 3.0);
        let v2 = Tuple::vector(2.0, 3.0, 4.0);

        assert_relative_eq!(v1.dot(&v2), 20.0);
    }

    #[test]
    fn cross_product() {
        let v1 = Tuple::vector(1.0, 2.0, 3.0);
        let v2 = Tuple::vector(2.0, 3.0, 4.0);

        assert_abs_diff_eq!(v1.cross(&v2), Tuple::vector(-1.0, 2.0, -1.0));
        assert_abs_diff_eq!(v2.cross(&v1), Tuple::vector(1.0, -2.0, 1.0));
    }

    #[test]
    fn reflect_vector_45() {
        let v = Tuple::vector(1.0, -1.0, 0.0);
        let n = Tuple::vector(0.0, 1.0, 0.0);

        assert_abs_diff_eq!(v.reflect(&n), Tuple::vector(1.0, 1.0, 0.0));
    }

    #[test]
    fn reflect_vector_slanted_surface() {
        let v = Tuple::vector(0.0, -1.0, 0.0);
        let n = Tuple::vector(2.0_f64.sqrt() / 2.0, 2.0_f64.sqrt() / 2.0, 0.0);

        assert_abs_diff_eq!(v.reflect(&n), Tuple::vector(1.0, 0.0, 0.0));
    }
}
