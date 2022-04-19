use std::ops;

use approx::AbsDiffEq;

use super::tuple::Tuple;

#[derive(Debug, Clone, Copy)]
pub struct Matrix<const N: usize>([[f64; N]; N]);

impl<const N: usize> Matrix<N> {
    pub fn identity() -> Self {
        let mut m = Self::zero();

        for i in 0..N {
            m[i][i] = 1.0;
        }

        m
    }

    fn zero() -> Self {
        Self([[0.0; N]; N])
    }

    pub fn transpose(&self) -> Self {
        let mut m = Self::zero();

        for i in 0..N {
            for j in 0..N {
                m[i][j] = self[j][i];
            }
        }

        m
    }
}

impl Matrix<2> {
    pub fn determinant(&self) -> f64 {
        self[0][0] * self[1][1] - self[0][1] * self[1][0]
    }
}

impl Matrix<3> {
    pub fn submatrix(&self, row: usize, col: usize) -> Matrix<2> {
        let mut m = Matrix::<2>::zero();

        let mut i = 0;
        for r in 0..3 {
            if r == row {
                continue;
            }

            let mut j = 0;
            for c in 0..3 {
                if c == col {
                    continue;
                }

                m[i][j] = self[r][c];
                j += 1;
            }

            i += 1;
        }

        m
    }

    pub fn minor(&self, row: usize, col: usize) -> f64 {
        self.submatrix(row, col).determinant()
    }

    pub fn cofactor(&self, row: usize, col: usize) -> f64 {
        let sign = if (row + col) % 2 == 0 { 1.0 } else { -1.0 };

        sign * self.minor(row, col)
    }

    pub fn determinant(&self) -> f64 {
        let mut det = 0.0;

        for column in 0..3 {
            det += self[0][column] * self.cofactor(0, column);
        }

        det
    }
}

impl Matrix<4> {
    pub fn submatrix(&self, row: usize, col: usize) -> Matrix<3> {
        let mut m = Matrix::<3>::zero();

        let mut i = 0;
        for r in 0..4 {
            if r == row {
                continue;
            }

            let mut j = 0;
            for c in 0..4 {
                if c == col {
                    continue;
                }

                m[i][j] = self[r][c];
                j += 1;
            }

            i += 1;
        }

        m
    }

    pub fn minor(&self, row: usize, col: usize) -> f64 {
        self.submatrix(row, col).determinant()
    }

    pub fn cofactor(&self, row: usize, col: usize) -> f64 {
        let sign = if (row + col) % 2 == 0 { 1.0 } else { -1.0 };

        sign * self.minor(row, col)
    }

    pub fn determinant(&self) -> f64 {
        let mut det = 0.0;

        for column in 0..4 {
            det += self[0][column] * self.cofactor(0, column);
        }

        det
    }

    pub fn inverse(&self) -> Self {
        let determinant = self.determinant();
        assert!(!(determinant == 0.0), "matrix is not invertible");

        let mut m = Self::zero();

        for row in 0..4 {
            for col in 0..4 {
                let cofactor = self.cofactor(row, col);
                m[col][row] = cofactor / determinant;
            }
        }

        m
    }

    pub fn translation(x: f64, y: f64, z: f64) -> Self {
        Self([
            [1.0, 0.0, 0.0, x],
            [0.0, 1.0, 0.0, y],
            [0.0, 0.0, 1.0, z],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }

    pub fn scaling(x: f64, y: f64, z: f64) -> Self {
        Self([
            [x, 0.0, 0.0, 0.0],
            [0.0, y, 0.0, 0.0],
            [0.0, 0.0, z, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }

    pub fn rotation_x(radians: f64) -> Self {
        let cos_r = radians.cos();
        let sin_r = radians.sin();

        Self([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, cos_r, -sin_r, 0.0],
            [0.0, sin_r, cos_r, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }

    pub fn rotation_y(radians: f64) -> Self {
        let cos_r = radians.cos();
        let sin_r = radians.sin();

        Self([
            [cos_r, 0.0, sin_r, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [-sin_r, 0.0, cos_r, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }

    pub fn rotation_z(radians: f64) -> Self {
        let cos_r = radians.cos();
        let sin_r = radians.sin();

        Self([
            [cos_r, -sin_r, 0.0, 0.0],
            [sin_r, cos_r, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }

    pub fn shearing(xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64) -> Self {
        Self([
            [1.0, xy, xz, 0.0],
            [yx, 1.0, yz, 0.0],
            [zx, zy, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }

    pub fn view_transform(from: Tuple, to: Tuple, up: Tuple) -> Self {
        let forward = (to - from).normalize();
        let upn = up.normalize();
        let left = forward.cross(&upn);
        let true_up = left.cross(&forward);

        let orientation = Self([
            [left.x, left.y, left.z, 0.0],
            [true_up.x, true_up.y, true_up.z, 0.0],
            [-forward.x, -forward.y, -forward.z, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);

        orientation * Self::translation(-from.x, -from.y, -from.z)
    }
}

impl<const N: usize> ops::Index<usize> for Matrix<N> {
    type Output = [f64; N];

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<const N: usize> ops::IndexMut<usize> for Matrix<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<const N: usize> PartialEq for Matrix<N> {
    fn eq(&self, other: &Self) -> bool {
        self.abs_diff_eq(other, Self::default_epsilon())
    }
}

impl<const N: usize> AbsDiffEq for Matrix<N> {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-5
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0
            .iter()
            .zip(other.0.iter())
            .all(|(a, b)| a.abs_diff_eq(b, epsilon))
    }
}

impl<const N: usize> ops::Mul for Matrix<N> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        let mut result = Self([[0.0; N]; N]);

        for i in 0..N {
            for k in 0..N {
                for j in 0..N {
                    result[i][j] += self[i][k] * other[k][j];
                }
            }
        }

        result
    }
}

impl ops::Mul<Tuple> for Matrix<4> {
    type Output = Tuple;

    fn mul(self, other: Tuple) -> Self::Output {
        let mut result = Tuple::new(0.0, 0.0, 0.0, 0.0);

        for i in 0..4 {
            for j in 0..4 {
                result[i] += self[i][j] * other[j];
            }
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

    use approx::{assert_abs_diff_eq, assert_abs_diff_ne, assert_relative_eq};

    use super::*;

    #[test]
    fn construct_4x4_matrix() {
        let m = Matrix([
            [1.0, 2.0, 3.0, 4.0],
            [5.5, 6.5, 7.5, 8.5],
            [9.0, 10.0, 11.0, 12.0],
            [13.5, 14.5, 15.5, 16.5],
        ]);

        assert_relative_eq!(m[0][0], 1.0);
        assert_relative_eq!(m[0][3], 4.0);
        assert_relative_eq!(m[1][0], 5.5);
        assert_relative_eq!(m[1][2], 7.5);
        assert_relative_eq!(m[2][2], 11.0);
        assert_relative_eq!(m[3][0], 13.5);
        assert_relative_eq!(m[3][2], 15.5);
    }

    #[test]
    fn construct_3x3_matrix() {
        let m = Matrix([[-3.0, 5.0, 0.0], [1.0, -2.0, -7.0], [0.0, 1.0, 1.0]]);

        assert_relative_eq!(m[0][0], -3.0);
        assert_relative_eq!(m[1][1], -2.0);
        assert_relative_eq!(m[2][2], 1.0);
    }

    #[test]
    fn construct_2x2_matrix() {
        let m = Matrix([[-3.0, 5.0], [1.0, -2.0]]);

        assert_relative_eq!(m[0][0], -3.0);
        assert_relative_eq!(m[0][1], 5.0);
        assert_relative_eq!(m[1][0], 1.0);
        assert_relative_eq!(m[1][1], -2.0);
    }

    #[test]
    fn equality_with_identical() {
        let m1 = Matrix([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);
        let m2 = Matrix([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);

        assert_abs_diff_eq!(m1, m2);
    }

    #[test]
    fn equality_with_different() {
        let m1 = Matrix([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);
        let m2 = Matrix([
            [2.0, 3.0, 4.0, 5.0],
            [6.0, 7.0, 8.0, 9.0],
            [8.0, 7.0, 6.0, 5.0],
            [4.0, 3.0, 2.0, 1.0],
        ]);

        assert_abs_diff_ne!(m1, m2);
    }

    #[test]
    fn multiply_matrices() {
        let a = Matrix([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);

        let b = Matrix([
            [-2.0, 1.0, 2.0, 3.0],
            [3.0, 2.0, 1.0, -1.0],
            [4.0, 3.0, 6.0, 5.0],
            [1.0, 2.0, 7.0, 8.0],
        ]);

        assert_abs_diff_eq!(
            a * b,
            Matrix([
                [20.0, 22.0, 50.0, 48.0],
                [44.0, 54.0, 114.0, 108.0],
                [40.0, 58.0, 110.0, 102.0],
                [16.0, 26.0, 46.0, 42.0],
            ])
        );
    }

    #[test]
    fn multiply_matrix_and_tuple() {
        let a = Matrix([
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 4.0, 2.0],
            [8.0, 6.0, 4.0, 1.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);
        let t = Tuple::new(1.0, 2.0, 3.0, 1.0);

        assert_abs_diff_eq!(a * t, Tuple::new(18.0, 24.0, 33.0, 1.0));
    }

    #[test]
    fn multiply_matrix_and_identity() {
        let a = Matrix([
            [0.0, 1.0, 2.0, 4.0],
            [1.0, 2.0, 4.0, 8.0],
            [2.0, 4.0, 8.0, 16.0],
            [4.0, 8.0, 16.0, 32.0],
        ]);

        assert_abs_diff_eq!(a * Matrix::identity(), a);
    }

    #[test]
    fn transpose_matrix() {
        let a = Matrix([
            [0.0, 9.0, 3.0, 0.0],
            [9.0, 8.0, 0.0, 8.0],
            [1.0, 8.0, 5.0, 3.0],
            [0.0, 0.0, 5.0, 8.0],
        ]);

        assert_abs_diff_eq!(
            a.transpose(),
            Matrix([
                [0.0, 9.0, 1.0, 0.0],
                [9.0, 8.0, 8.0, 0.0],
                [3.0, 0.0, 5.0, 5.0],
                [0.0, 8.0, 3.0, 8.0],
            ])
        );
    }

    #[test]
    fn determinant_2x2() {
        let a = Matrix([[1.0, 5.0], [-3.0, 2.0]]);

        assert_relative_eq!(a.determinant(), 17.0);
    }

    #[test]
    fn submatrix_3x3() {
        let a = Matrix([[1.0, 5.0, 0.0], [-3.0, 2.0, 7.0], [0.0, 6.0, -3.0]]);

        assert_abs_diff_eq!(a.submatrix(0, 2), Matrix([[-3.0, 2.0], [0.0, 6.0]]));
    }

    #[test]
    fn submatrix_4x4() {
        let a = Matrix([
            [-6.0, 1.0, 1.0, 6.0],
            [-8.0, 5.0, 8.0, 6.0],
            [-1.0, 0.0, 8.0, 2.0],
            [-7.0, 1.0, -1.0, 1.0],
        ]);

        assert_abs_diff_eq!(
            a.submatrix(2, 1),
            Matrix([[-6.0, 1.0, 6.0], [-8.0, 8.0, 6.0], [-7.0, -1.0, 1.0]])
        );
    }

    #[test]
    fn minor_3x3() {
        let a = Matrix([[3.0, 5.0, 0.0], [2.0, -1.0, -7.0], [6.0, -1.0, 5.0]]);

        let b = a.submatrix(1, 0);

        assert_relative_eq!(b.determinant(), 25.0);
        assert_relative_eq!(a.minor(1, 0), 25.0);
    }

    #[test]
    fn cofactor_3x3() {
        let a = Matrix([[3.0, 5.0, 0.0], [2.0, -1.0, -7.0], [6.0, -1.0, 5.0]]);

        assert_relative_eq!(a.minor(0, 0), -12.0);
        assert_relative_eq!(a.cofactor(0, 0), -12.0);
        assert_relative_eq!(a.minor(1, 0), 25.0);
        assert_relative_eq!(a.cofactor(1, 0), -25.0);
    }

    #[test]
    fn determinant_3x3() {
        let a = Matrix([[1.0, 2.0, 6.0], [-5.0, 8.0, -4.0], [2.0, 6.0, 4.0]]);

        assert_relative_eq!(a.cofactor(0, 0), 56.0);
        assert_relative_eq!(a.cofactor(0, 1), 12.0);
        assert_relative_eq!(a.cofactor(0, 2), -46.0);
        assert_relative_eq!(a.determinant(), -196.0);
    }

    #[test]
    fn determinant_4x4() {
        let a = Matrix([
            [-2.0, -8.0, 3.0, 5.0],
            [-3.0, 1.0, 7.0, 3.0],
            [1.0, 2.0, -9.0, 6.0],
            [-6.0, 7.0, 7.0, -9.0],
        ]);

        assert_relative_eq!(a.cofactor(0, 0), 690.0);
        assert_relative_eq!(a.cofactor(0, 1), 447.0);
        assert_relative_eq!(a.cofactor(0, 2), 210.0);
        assert_relative_eq!(a.cofactor(0, 3), 51.0);
        assert_relative_eq!(a.determinant(), -4071.0);
    }

    #[test]
    fn check_invertible() {
        let a = Matrix([
            [6.0, 4.0, 4.0, 4.0],
            [5.0, 5.0, 7.0, 6.0],
            [4.0, -9.0, 3.0, -7.0],
            [9.0, 1.0, 7.0, -6.0],
        ]);

        assert_relative_eq!(a.determinant(), -2120.0);
    }

    #[test]
    fn check_noninvertible() {
        let a = Matrix([
            [-4.0, 2.0, -2.0, -3.0],
            [9.0, 6.0, 2.0, 6.0],
            [0.0, -5.0, 1.0, -5.0],
            [0.0, 0.0, 0.0, 0.0],
        ]);

        assert_relative_eq!(a.determinant(), 0.0);
    }

    #[test]
    fn inverse_4x4() {
        let a = Matrix([
            [-5.0, 2.0, 6.0, -8.0],
            [1.0, -5.0, 1.0, 8.0],
            [7.0, 7.0, -6.0, -7.0],
            [1.0, -3.0, 7.0, 4.0],
        ]);

        let inv = a.inverse();

        assert_relative_eq!(a.determinant(), 532.0);
        assert_relative_eq!(a.cofactor(2, 3), -160.0);
        assert_relative_eq!(inv[3][2], -160.0 / 532.0);
        assert_relative_eq!(a.cofactor(3, 2), 105.0);
        assert_relative_eq!(inv[2][3], 105.0 / 532.0);
        assert_abs_diff_eq!(
            inv,
            Matrix([
                [0.21805, 0.45113, 0.24060, -0.04511],
                [-0.80827, -1.45677, -0.44361, 0.52068],
                [-0.07895, -0.22368, -0.05263, 0.19737],
                [-0.52256, -0.81391, -0.30075, 0.30639],
            ])
        );
    }

    #[test]
    fn inverse_4x4_2() {
        let a = Matrix([
            [8.0, -5.0, 9.0, 2.0],
            [7.0, 5.0, 6.0, 1.0],
            [-6.0, 0.0, 9.0, 6.0],
            [-3.0, 0.0, -9.0, -4.0],
        ]);

        assert_abs_diff_eq!(
            a.inverse(),
            Matrix([
                [-0.15385, -0.15385, -0.28205, -0.53846],
                [-0.07692, 0.12308, 0.02564, 0.03077],
                [0.35897, 0.35897, 0.43590, 0.92308],
                [-0.69231, -0.69231, -0.76923, -1.92308],
            ])
        );
    }

    #[test]
    fn inverse_4x4_3() {
        let a = Matrix([
            [9.0, 3.0, 0.0, 9.0],
            [-5.0, -2.0, -6.0, -3.0],
            [-4.0, 9.0, 6.0, 4.0],
            [-7.0, 6.0, 6.0, 2.0],
        ]);

        assert_abs_diff_eq!(
            a.inverse(),
            Matrix([
                [-0.04074, -0.07778, 0.14444, -0.22222],
                [-0.07778, 0.03333, 0.36667, -0.33333],
                [-0.02901, -0.14630, -0.10926, 0.12963],
                [0.17778, 0.06667, -0.26667, 0.33333],
            ])
        );
    }

    #[test]
    fn multiply_product_inverse() {
        let a = Matrix([
            [3.0, -9.0, 7.0, 3.0],
            [3.0, -8.0, 2.0, -9.0],
            [-4.0, 4.0, 4.0, 1.0],
            [-6.0, 5.0, -1.0, 1.0],
        ]);

        let b = Matrix([
            [8.0, 2.0, 2.0, 2.0],
            [3.0, -1.0, 7.0, 0.0],
            [7.0, 0.0, 5.0, 4.0],
            [6.0, -2.0, 0.0, 5.0],
        ]);

        assert_abs_diff_eq!(a * b * b.inverse(), a);
    }

    #[test]
    fn invert_identity() {
        let indentity: Matrix<4> = Matrix::identity();

        assert_abs_diff_eq!(indentity.inverse(), Matrix::identity());
    }

    #[test]
    fn multiply_inverse() {
        let a = Matrix([
            [3.0, -9.0, 7.0, 3.0],
            [3.0, -8.0, 2.0, -9.0],
            [-4.0, 4.0, 4.0, 1.0],
            [-6.0, 5.0, -1.0, 1.0],
        ]);

        assert_abs_diff_eq!(a * a.inverse(), Matrix::identity());
    }

    #[test]
    fn inverse_transpose_transpose_inverse() {
        let a = Matrix([
            [3.0, -9.0, 7.0, 3.0],
            [3.0, -8.0, 2.0, -9.0],
            [-4.0, 4.0, 4.0, 1.0],
            [-6.0, 5.0, -1.0, 1.0],
        ]);

        assert_abs_diff_eq!(a.inverse().transpose(), a.transpose().inverse());
    }

    #[test]
    fn multiply_almost_identity_by_tuple() {
        let mut a: Matrix<4> = Matrix::identity();
        a[0][0] = -2.0;

        let b = a * Tuple::point(3.0, 2.0, 1.0);

        assert_relative_eq!(b.x, -6.0);
        assert_relative_eq!(b.y, 2.0);
        assert_relative_eq!(b.z, 1.0);
    }

    #[test]
    fn multiply_by_translation() {
        let a = Matrix::translation(5.0, -3.0, 2.0);
        let p = Tuple::point(-3.0, 4.0, 5.0);

        assert_abs_diff_eq!(a * p, Tuple::point(2.0, 1.0, 7.0));
    }

    #[test]
    fn multiply_by_inverse_translation() {
        let a = Matrix::translation(5.0, -3.0, 2.0);
        let p = Tuple::point(-3.0, 4.0, 5.0);

        assert_abs_diff_eq!(a.inverse() * p, Tuple::point(-8.0, 7.0, 3.0));
    }

    #[test]
    fn translation_doesnt_affect_vectors() {
        let a = Matrix::translation(5.0, -3.0, 2.0);
        let v = Tuple::vector(-3.0, 4.0, 5.0);

        assert_abs_diff_eq!(a * v, v);
    }

    #[test]
    fn scaling_matrix_applied_to_point() {
        let a = Matrix::scaling(2.0, 3.0, 4.0);
        let p = Tuple::point(-4.0, 6.0, 8.0);

        assert_abs_diff_eq!(a * p, Tuple::point(-8.0, 18.0, 32.0));
    }

    #[test]
    fn scaling_matrix_applied_to_vector() {
        let a = Matrix::scaling(2.0, 3.0, 4.0);
        let v = Tuple::vector(-4.0, 6.0, 8.0);

        assert_abs_diff_eq!(a * v, Tuple::vector(-8.0, 18.0, 32.0));
    }

    #[test]
    fn multiply_inverse_of_scaling_matrix() {
        let a = Matrix::scaling(2.0, 3.0, 4.0);
        let p = Tuple::point(-4.0, 6.0, 8.0);

        assert_abs_diff_eq!(a.inverse() * p, Tuple::point(-2.0, 2.0, 2.0));
    }

    #[test]
    fn reflection_is_scaling_by_negative_value() {
        let a = Matrix::scaling(-1.0, 1.0, 1.0);
        let p = Tuple::point(2.0, 3.0, 4.0);

        assert_abs_diff_eq!(a * p, Tuple::point(-2.0, 3.0, 4.0));
    }

    #[test]
    fn rotation_around_x_axis() {
        let p = Tuple::point(0.0, 1.0, 0.0);
        let half_quarter = Matrix::rotation_x(FRAC_PI_4);
        let full_quarter = Matrix::rotation_x(FRAC_PI_2);

        assert_abs_diff_eq!(
            half_quarter * p,
            Tuple::point(0.0, 2.0_f64.sqrt() / 2.0, 2.0_f64.sqrt() / 2.0)
        );
        assert_abs_diff_eq!(full_quarter * p, Tuple::point(0.0, 0.0, 1.0));
    }

    #[test]
    fn rotation_around_y_axis() {
        let p = Tuple::point(0.0, 0.0, 1.0);
        let half_quarter = Matrix::rotation_y(FRAC_PI_4);
        let full_quarter = Matrix::rotation_y(FRAC_PI_2);

        assert_abs_diff_eq!(
            half_quarter * p,
            Tuple::point(2.0_f64.sqrt() / 2.0, 0.0, 2.0_f64.sqrt() / 2.0)
        );
        assert_abs_diff_eq!(full_quarter * p, Tuple::point(1.0, 0.0, 0.0));
    }

    #[test]
    fn rotation_around_z_axis() {
        let p = Tuple::point(0.0, 1.0, 0.0);
        let half_quarter = Matrix::rotation_z(FRAC_PI_4);
        let full_quarter = Matrix::rotation_z(FRAC_PI_2);

        assert_abs_diff_eq!(
            half_quarter * p,
            Tuple::point(-(2.0_f64.sqrt()) / 2.0, 2.0_f64.sqrt() / 2.0, 0.0)
        );
        assert_abs_diff_eq!(full_quarter * p, Tuple::point(-1.0, 0.0, 0.0));
    }

    #[test]
    fn shearing_moves_x_in_proportion_to_y() {
        let a = Matrix::shearing(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let p = Tuple::point(2.0, 3.0, 4.0);

        assert_abs_diff_eq!(a * p, Tuple::point(5.0, 3.0, 4.0));
    }

    #[test]
    fn shearing_moves_x_in_proportion_to_z() {
        let a = Matrix::shearing(0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
        let p = Tuple::point(2.0, 3.0, 4.0);

        assert_abs_diff_eq!(a * p, Tuple::point(6.0, 3.0, 4.0));
    }

    #[test]
    fn shearing_moves_y_in_proportion_to_x() {
        let a = Matrix::shearing(0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
        let p = Tuple::point(2.0, 3.0, 4.0);

        assert_abs_diff_eq!(a * p, Tuple::point(2.0, 5.0, 4.0));
    }

    #[test]
    fn shearing_moves_y_in_proportion_to_z() {
        let a = Matrix::shearing(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        let p = Tuple::point(2.0, 3.0, 4.0);

        assert_abs_diff_eq!(a * p, Tuple::point(2.0, 7.0, 4.0));
    }

    #[test]
    fn shearing_moves_z_in_proportion_to_x() {
        let a = Matrix::shearing(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        let p = Tuple::point(2.0, 3.0, 4.0);

        assert_abs_diff_eq!(a * p, Tuple::point(2.0, 3.0, 6.0));
    }

    #[test]
    fn shearing_moves_z_in_proportion_to_y() {
        let a = Matrix::shearing(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let p = Tuple::point(2.0, 3.0, 4.0);

        assert_abs_diff_eq!(a * p, Tuple::point(2.0, 3.0, 7.0));
    }

    #[test]
    fn transformations_applied_in_sequence() {
        let p = Tuple::point(1.0, 0.0, 1.0);
        let a = Matrix::rotation_x(FRAC_PI_2);
        let b = Matrix::scaling(5.0, 5.0, 5.0);
        let c = Matrix::translation(10.0, 5.0, 7.0);

        let p2 = a * p;
        assert_abs_diff_eq!(p2, Tuple::point(1.0, -1.0, 0.0));

        let p3 = b * p2;
        assert_abs_diff_eq!(p3, Tuple::point(5.0, -5.0, 0.0));

        let p4 = c * p3;
        assert_abs_diff_eq!(p4, Tuple::point(15.0, 0.0, 7.0));
    }

    #[test]
    fn chained_transformations() {
        let p = Tuple::point(1.0, 0.0, 1.0);
        let a = Matrix::rotation_x(FRAC_PI_2);
        let b = Matrix::scaling(5.0, 5.0, 5.0);
        let c = Matrix::translation(10.0, 5.0, 7.0);

        let t = c * b * a;

        assert_abs_diff_eq!(t * p, Tuple::point(15.0, 0.0, 7.0));
    }

    #[test]
    fn view_transformation_default_orientation() {
        let from = Tuple::point(0.0, 0.0, 0.0);
        let to = Tuple::point(0.0, 0.0, -1.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);

        let t = Matrix::view_transform(from, to, up);

        assert_abs_diff_eq!(t, Matrix::identity());
    }

    #[test]
    fn view_transformation_positive_z_direction() {
        let from = Tuple::point(0.0, 0.0, 0.0);
        let to = Tuple::point(0.0, 0.0, 1.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);

        let t = Matrix::view_transform(from, to, up);

        assert_abs_diff_eq!(t, Matrix::scaling(-1.0, 1.0, -1.0));
    }

    #[test]
    fn view_transformation_moves_world() {
        let from = Tuple::point(0.0, 0.0, 8.0);
        let to = Tuple::point(0.0, 0.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);

        let t = Matrix::view_transform(from, to, up);

        assert_abs_diff_eq!(t, Matrix::translation(0.0, 0.0, -8.0));
    }

    #[test]
    fn arbitrary_view_transformation() {
        let from = Tuple::point(1.0, 3.0, 2.0);
        let to = Tuple::point(4.0, -2.0, 8.0);
        let up = Tuple::vector(1.0, 1.0, 0.0);

        let t = Matrix::view_transform(from, to, up);

        assert_abs_diff_eq!(
            t,
            Matrix([
                [-0.50709, 0.50709, 0.67612, -2.36643],
                [0.76772, 0.60609, 0.12122, -2.82843],
                [-0.35857, 0.59761, -0.71714, 0.00000],
                [0.00000, 0.00000, 0.00000, 1.00000]
            ])
        );
    }
}
