use std::ops;

use approx::AbsDiffEq;

use super::tuple::Tuple;

#[derive(Debug, PartialEq, Clone, Copy)]
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

impl<const N: usize> AbsDiffEq for Matrix<N> {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-6
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

        let result_t = a * t;

        assert_relative_eq!(result_t.x, 18.0);
        assert_relative_eq!(result_t.y, 24.0);
        assert_relative_eq!(result_t.z, 33.0);
        assert_relative_eq!(result_t.w, 1.0);
    }

    #[test]
    fn multiply_matrix_and_identity() {
        let a = Matrix([
            [0.0, 1.0, 2.0, 4.0],
            [1.0, 2.0, 4.0, 8.0],
            [2.0, 4.0, 8.0, 16.0],
            [4.0, 8.0, 16.0, 32.0],
        ]);

        let result = a * Matrix::identity();

        assert_abs_diff_eq!(result, a);
    }

    #[test]
    fn transpose_matrix() {
        let a = Matrix([
            [0.0, 9.0, 3.0, 0.0],
            [9.0, 8.0, 0.0, 8.0],
            [1.0, 8.0, 5.0, 3.0],
            [0.0, 0.0, 5.0, 8.0],
        ]);

        let result = a.transpose();

        assert_abs_diff_eq!(
            result,
            Matrix([
                [0.0, 9.0, 1.0, 0.0],
                [9.0, 8.0, 8.0, 0.0],
                [3.0, 0.0, 5.0, 5.0],
                [0.0, 8.0, 3.0, 8.0],
            ])
        );
    }
}
