use std::ops;

pub struct Matrix<const N: usize>([[f64; N]; N]);

impl<const N: usize> ops::Index<usize> for Matrix<N> {
    type Output = [f64; N];

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

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
}
