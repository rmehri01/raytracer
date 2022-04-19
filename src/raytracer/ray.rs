use crate::core::{matrix::Matrix, point::Point, vector::Vector};

/// A ray in 3D space with a origin position and direction vector.
#[derive(Debug)]
pub struct Ray {
    pub origin: Point,
    pub direction: Vector,
}

impl Ray {
    pub fn new(origin: Point, direction: Vector) -> Self {
        Self { origin, direction }
    }

    pub fn position(&self, t: f64) -> Point {
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
        let origin = Point::new(1.0, 2.0, 3.0);
        let direction = Vector::new(4.0, 5.0, 6.0);
        let r = Ray::new(origin, direction);

        assert_abs_diff_eq!(r.origin, origin);
        assert_abs_diff_eq!(r.direction, direction);
    }

    #[test]
    fn point_from_distance() {
        let r = Ray::new(Point::new(2.0, 3.0, 4.0), Vector::new(1.0, 0.0, 0.0));

        assert_abs_diff_eq!(r.position(0.0), Point::new(2.0, 3.0, 4.0));
        assert_abs_diff_eq!(r.position(1.0), Point::new(3.0, 3.0, 4.0));
        assert_abs_diff_eq!(r.position(-1.0), Point::new(1.0, 3.0, 4.0));
        assert_abs_diff_eq!(r.position(2.5), Point::new(4.5, 3.0, 4.0));
    }

    #[test]
    fn ray_translation() {
        let r = Ray::new(Point::new(1.0, 2.0, 3.0), Vector::new(0.0, 1.0, 0.0));
        let m = Matrix::translation(3.0, 4.0, 5.0);

        let r2 = r.transform(&m);

        assert_abs_diff_eq!(r2.origin, Point::new(4.0, 6.0, 8.0));
        assert_abs_diff_eq!(r2.direction, Vector::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn ray_scaling() {
        let r = Ray::new(Point::new(1.0, 2.0, 3.0), Vector::new(0.0, 1.0, 0.0));
        let m = Matrix::scaling(2.0, 3.0, 4.0);

        let r2 = r.transform(&m);

        assert_abs_diff_eq!(r2.origin, Point::new(2.0, 6.0, 12.0));
        assert_abs_diff_eq!(r2.direction, Vector::new(0.0, 3.0, 0.0));
    }
}
