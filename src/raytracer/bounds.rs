use crate::core::{matrix::Transformation, point::Point};

use super::ray::Ray;

/// An axis-aligned bounding box that can be used to quickly determine if a ray
/// might intersect with anything in the box.
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
        let mut bounds = Self::default();
        let corners = [
            self.minimum,
            Point::new(self.minimum.x, self.minimum.y, self.maximum.z),
            Point::new(self.minimum.x, self.maximum.y, self.minimum.z),
            Point::new(self.maximum.x, self.minimum.y, self.minimum.z),
            Point::new(self.minimum.x, self.maximum.y, self.maximum.z),
            Point::new(self.maximum.x, self.minimum.y, self.maximum.z),
            Point::new(self.maximum.x, self.maximum.y, self.minimum.z),
            self.maximum,
        ];

        for corner in &corners {
            let transformed_point = *transform * *corner;
            bounds.add_point(&transformed_point);
        }
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

    // TODO: duplicated
    pub fn intersects(&self, ray: &Ray) -> bool {
        let (x_t_min, x_t_max) = Self::check_axis(
            ray.origin.x,
            ray.direction.x,
            self.minimum.x,
            self.maximum.x,
        );
        let (y_t_min, y_t_max) = Self::check_axis(
            ray.origin.y,
            ray.direction.y,
            self.minimum.y,
            self.maximum.y,
        );
        let (z_t_min, z_t_max) = Self::check_axis(
            ray.origin.z,
            ray.direction.z,
            self.minimum.z,
            self.maximum.z,
        );

        let t_min = x_t_min.max(y_t_min).max(z_t_min);
        let t_max = x_t_max.min(y_t_max).min(z_t_max);

        t_min <= t_max
    }

    fn check_axis(origin: f64, direction: f64, min: f64, max: f64) -> (f64, f64) {
        let t_min_numerator = min - origin;
        let t_max_numerator = max - origin;

        let t_min = t_min_numerator / direction;
        let t_max = t_max_numerator / direction;

        if t_min > t_max {
            (t_max, t_min)
        } else {
            (t_min, t_max)
        }
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
    use super::*;

    // TODO: tests
}
