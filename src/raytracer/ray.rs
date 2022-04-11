use crate::core::{matrix::Matrix, tuple::Tuple};

/// A ray in 3D space with a origin position and direction vector.
#[derive(Debug)]
pub struct Ray {
    pub origin: Tuple,
    pub direction: Tuple,
}

impl Ray {
    pub fn new(origin: Tuple, direction: Tuple) -> Self {
        Self { origin, direction }
    }

    pub fn position(&self, t: f64) -> Tuple {
        self.origin + self.direction * t
    }

    pub fn transform(&self, m: &Matrix<4>) -> Self {
        Self {
            origin: *m * self.origin,
            direction: *m * self.direction,
        }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn ray_new() {
        let origin = Tuple::point(1.0, 2.0, 3.0);
        let direction = Tuple::vector(4.0, 5.0, 6.0);
        let r = Ray::new(origin, direction);

        assert_abs_diff_eq!(r.origin, origin);
        assert_abs_diff_eq!(r.direction, direction);
    }

    #[test]
    fn point_from_distance() {
        let r = Ray::new(Tuple::point(2.0, 3.0, 4.0), Tuple::vector(1.0, 0.0, 0.0));

        assert_abs_diff_eq!(r.position(0.0), Tuple::point(2.0, 3.0, 4.0));
        assert_abs_diff_eq!(r.position(1.0), Tuple::point(3.0, 3.0, 4.0));
        assert_abs_diff_eq!(r.position(-1.0), Tuple::point(1.0, 3.0, 4.0));
        assert_abs_diff_eq!(r.position(2.5), Tuple::point(4.5, 3.0, 4.0));
    }

    #[test]
    fn ray_translation() {
        let r = Ray::new(Tuple::point(1.0, 2.0, 3.0), Tuple::vector(0.0, 1.0, 0.0));
        let m = Matrix::translation(3.0, 4.0, 5.0);

        let r2 = r.transform(&m);

        assert_abs_diff_eq!(r2.origin, Tuple::point(4.0, 6.0, 8.0));
        assert_abs_diff_eq!(r2.direction, Tuple::vector(0.0, 1.0, 0.0));
    }

    #[test]
    fn ray_scaling() {
        let r = Ray::new(Tuple::point(1.0, 2.0, 3.0), Tuple::vector(0.0, 1.0, 0.0));
        let m = Matrix::scaling(2.0, 3.0, 4.0);

        let r2 = r.transform(&m);

        assert_abs_diff_eq!(r2.origin, Tuple::point(2.0, 6.0, 12.0));
        assert_abs_diff_eq!(r2.direction, Tuple::vector(0.0, 3.0, 0.0));
    }
}
