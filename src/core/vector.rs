use std::ops;

use approx::AbsDiffEq;

use super::point::Point;

#[derive(Debug, Clone, Copy)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn magnitude(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }

    // TODO: mutate?
    pub fn normalize(&self) -> Self {
        let m = self.magnitude();

        Self::new(self.x / m, self.y / m, self.z / m)
    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(&self, other: &Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    pub fn reflect(&self, normal: &Self) -> Self {
        *self - *normal * 2.0 * self.dot(normal)
    }
}

impl ops::Add<Self> for Vector {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self::Output::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl ops::Add<Point> for Vector {
    type Output = Point;

    fn add(self, other: Point) -> Self::Output {
        Self::Output::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl ops::Sub for Vector {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl ops::Neg for Vector {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl ops::Mul<f64> for Vector {
    type Output = Self;

    fn mul(self, other: f64) -> Self {
        Self::new(self.x * other, self.y * other, self.z * other)
    }
}

impl ops::Div<f64> for Vector {
    type Output = Self;

    fn div(self, other: f64) -> Self {
        Self::new(self.x / other, self.y / other, self.z / other)
    }
}

impl PartialEq for Vector {
    fn eq(&self, other: &Self) -> bool {
        self.abs_diff_eq(other, Self::default_epsilon())
    }
}

impl AbsDiffEq for Vector {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-4
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.x.abs_diff_eq(&other.x, epsilon)
            && self.y.abs_diff_eq(&other.y, epsilon)
            && self.z.abs_diff_eq(&other.z, epsilon)
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn sub_two_vectors() {
        let v1 = Vector::new(3.0, 2.0, 1.0);
        let v2 = Vector::new(5.0, 6.0, 7.0);

        assert_abs_diff_eq!(v1 - v2, Vector::new(-2.0, -4.0, -6.0));
    }

    #[test]
    fn sub_vector_zero_vector() {
        let zero = Vector::new(0.0, 0.0, 0.0);
        let v = Vector::new(1.0, -2.0, 3.0);

        assert_abs_diff_eq!(zero - v, Vector::new(-1.0, 2.0, -3.0));
    }

    #[test]
    fn vector_neg() {
        let a = Vector::new(1.0, -2.0, 3.0);

        assert_abs_diff_eq!(-a, Vector::new(-1.0, 2.0, -3.0));
    }

    #[test]
    fn mul_scalar() {
        let a = Vector::new(1.0, -2.0, 3.0);

        assert_abs_diff_eq!(a * 3.5, Vector::new(3.5, -7.0, 10.5));
    }

    #[test]
    fn mul_fraction() {
        let a = Vector::new(1.0, -2.0, 3.0);

        assert_abs_diff_eq!(a * 0.5, Vector::new(0.5, -1.0, 1.5));
    }

    #[test]
    fn div_scalar() {
        let a = Vector::new(1.0, -2.0, 3.0);

        assert_abs_diff_eq!(a / 2.0, Vector::new(0.5, -1.0, 1.5));
    }

    #[test]
    fn magnitude_x() {
        let v = Vector::new(1.0, 0.0, 0.0);

        assert_abs_diff_eq!(v.magnitude(), 1.0);
    }

    #[test]
    fn magnitude_y() {
        let v = Vector::new(0.0, 1.0, 0.0);

        assert_abs_diff_eq!(v.magnitude(), 1.0);
    }

    #[test]
    fn magnitude_z() {
        let v = Vector::new(0.0, 0.0, 1.0);

        assert_abs_diff_eq!(v.magnitude(), 1.0);
    }

    #[test]
    fn magnitude_all() {
        let v = Vector::new(1.0, 2.0, 3.0);

        assert_abs_diff_eq!(v.magnitude(), 14.0_f64.sqrt());
    }

    #[test]
    fn magnitude_negative() {
        let v = Vector::new(-1.0, -2.0, -3.0);

        assert_abs_diff_eq!(v.magnitude(), 14.0_f64.sqrt());
    }

    #[test]
    fn normalize_x() {
        let v = Vector::new(4.0, 0.0, 0.0);

        assert_abs_diff_eq!(v.normalize(), Vector::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn normalize_all() {
        let v = Vector::new(1.0, 2.0, 3.0);

        assert_abs_diff_eq!(
            v.normalize(),
            Vector::new(
                1.0 / 14.0_f64.sqrt(),
                2.0 / 14.0_f64.sqrt(),
                3.0 / 14.0_f64.sqrt()
            )
        );
    }

    #[test]
    fn magnitude_of_normalized() {
        let v = Vector::new(1.0, 2.0, 3.0);

        assert_abs_diff_eq!(v.normalize().magnitude(), 1.0);
    }

    #[test]
    fn dot_product() {
        let v1 = Vector::new(1.0, 2.0, 3.0);
        let v2 = Vector::new(2.0, 3.0, 4.0);

        assert_abs_diff_eq!(v1.dot(&v2), 20.0);
    }

    #[test]
    fn cross_product() {
        let v1 = Vector::new(1.0, 2.0, 3.0);
        let v2 = Vector::new(2.0, 3.0, 4.0);

        assert_abs_diff_eq!(v1.cross(&v2), Vector::new(-1.0, 2.0, -1.0));
        assert_abs_diff_eq!(v2.cross(&v1), Vector::new(1.0, -2.0, 1.0));
    }

    #[test]
    fn reflect_vector_45() {
        let v = Vector::new(1.0, -1.0, 0.0);
        let n = Vector::new(0.0, 1.0, 0.0);

        assert_abs_diff_eq!(v.reflect(&n), Vector::new(1.0, 1.0, 0.0));
    }

    #[test]
    fn reflect_vector_slanted_surface() {
        let v = Vector::new(0.0, -1.0, 0.0);
        let n = Vector::new(2.0_f64.sqrt() / 2.0, 2.0_f64.sqrt() / 2.0, 0.0);

        assert_abs_diff_eq!(v.reflect(&n), Vector::new(1.0, 0.0, 0.0));
    }
}
