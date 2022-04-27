use crate::core::{matrix::Transformation, point::Point};

use super::{ray::Ray, shapes::check_axis};

/// An axis-aligned bounding box that can be used to quickly determine if a ray
/// might intersect with anything in the box.
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Bounds {
    minimum: Point,
    maximum: Point,
}

impl Bounds {
    pub fn new(minimum: Point, maximum: Point) -> Self {
        Self { minimum, maximum }
    }

    /// Extends the bounding box to include all points in the transformed box.
    pub fn transform(&mut self, transform: &Transformation) {
        let corners = [
            self.minimum,
            Point::new(self.minimum.x, self.minimum.y, self.maximum.z),
            Point::new(self.minimum.x, self.maximum.y, self.minimum.z),
            Point::new(self.minimum.x, self.maximum.y, self.maximum.z),
            Point::new(self.maximum.x, self.minimum.y, self.minimum.z),
            Point::new(self.maximum.x, self.minimum.y, self.maximum.z),
            Point::new(self.maximum.x, self.maximum.y, self.minimum.z),
            self.maximum,
        ];

        for corner in corners {
            let transformed_point = *transform * corner;
            self.add_point(&transformed_point);
        }
    }

    /// Combines this bounding box with another bounding box.
    pub fn union(&mut self, other: &Self) {
        self.add_point(&other.minimum);
        self.add_point(&other.maximum);
    }

    /// Extends the bounding box to include the given point.
    pub fn add_point(&mut self, point: &Point) {
        self.minimum.x = self.minimum.x.min(point.x);
        self.minimum.y = self.minimum.y.min(point.y);
        self.minimum.z = self.minimum.z.min(point.z);
        self.maximum.x = self.maximum.x.max(point.x);
        self.maximum.y = self.maximum.y.max(point.y);
        self.maximum.z = self.maximum.z.max(point.z);
    }

    pub fn intersects(&self, ray: &Ray) -> bool {
        let (x_t_min, x_t_max) = check_axis(
            ray.origin.x,
            ray.direction.x,
            self.minimum.x,
            self.maximum.x,
        );
        let (y_t_min, y_t_max) = check_axis(
            ray.origin.y,
            ray.direction.y,
            self.minimum.y,
            self.maximum.y,
        );
        let (z_t_min, z_t_max) = check_axis(
            ray.origin.z,
            ray.direction.z,
            self.minimum.z,
            self.maximum.z,
        );

        let t_min = x_t_min.max(y_t_min).max(z_t_min);
        let t_max = x_t_max.min(y_t_max).min(z_t_max);

        t_min <= t_max
    }
}

impl Default for Bounds {
    fn default() -> Self {
        Self {
            minimum: Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
            maximum: Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::FRAC_PI_4;

    use crate::core::{matrix::Matrix, vector::Vector};

    use super::*;

    #[test]
    fn add_point() {
        let mut bounds = Bounds::default();

        bounds.add_point(&Point::new(-3.0, 2.0, 4.0));
        bounds.add_point(&Point::new(-5.0, 4.0, 4.0));

        assert_eq!(bounds.minimum, Point::new(-5.0, 2.0, 4.0));
        assert_eq!(bounds.maximum, Point::new(-3.0, 4.0, 4.0));
    }

    #[test]
    fn union() {
        let mut bounds = Bounds::default();

        bounds.add_point(&Point::new(-3.0, 2.0, 4.0));
        bounds.add_point(&Point::new(-5.0, 4.0, 4.0));

        let mut other = Bounds::default();

        other.add_point(&Point::new(-10.0, 2.0, -1.0));
        other.add_point(&Point::new(-4.0, 4.0, 10.0));

        bounds.union(&other);

        assert_eq!(bounds.minimum, Point::new(-10.0, 2.0, -1.0));
        assert_eq!(bounds.maximum, Point::new(-3.0, 4.0, 10.0));
    }

    #[test]
    fn intersects() {
        let bounds = Bounds::new(Point::new(-3.0, 2.0, 4.0), Point::new(5.0, 4.0, 6.0));
        let ray = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 1.0, 1.0));

        assert!(bounds.intersects(&ray));
    }

    #[test]
    fn does_not_intersect() {
        let bounds = Bounds::new(Point::new(-3.0, 2.0, 4.0), Point::new(5.0, 4.0, 6.0));
        let ray = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));

        assert!(!bounds.intersects(&ray));
    }

    #[test]
    fn transform_translation() {
        let mut bounds = Bounds::new(Point::new(-1.0, -1.0, -1.0), Point::new(1.0, 1.0, 1.0));

        bounds.transform(&Matrix::translation(1.0, 0.0, 0.0));

        assert_eq!(bounds.minimum, Point::new(-1.0, -1.0, -1.0));
        assert_eq!(bounds.maximum, Point::new(2.0, 1.0, 1.0));
    }

    #[test]
    fn transform_rotation() {
        let mut bounds = Bounds::new(Point::new(-1.0, -1.0, -1.0), Point::new(1.0, 1.0, 1.0));

        bounds.transform(&Matrix::rotation_x(FRAC_PI_4));

        assert_eq!(
            bounds.minimum,
            Point::new(-1.0, -(2.0_f64.sqrt()), -(2.0_f64.sqrt()))
        );
        assert_eq!(
            bounds.maximum,
            Point::new(1.0, 2.0_f64.sqrt(), 2.0_f64.sqrt())
        );
    }
}
