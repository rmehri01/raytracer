use crate::core::{matrix::Matrix, tuple::Tuple};

use super::ray::Ray;

/// An axis-aligned bounding box that can be used to quickly determine if a ray
/// might intersect with anything in the box.
pub struct Bounds {
    minimum: Tuple,
    maximum: Tuple,
}

impl Bounds {
    pub fn new(minimum: Tuple, maximum: Tuple) -> Self {
        Self { minimum, maximum }
    }

    // TODO: this and transform not sure if better way
    pub fn add_point(&mut self, point: &Tuple) {
        self.minimum.x = self.minimum.x.min(point.x);
        self.minimum.y = self.minimum.y.min(point.y);
        self.minimum.z = self.minimum.z.min(point.z);
        self.maximum.x = self.maximum.x.max(point.x);
        self.maximum.y = self.maximum.y.max(point.y);
        self.maximum.z = self.maximum.z.max(point.z);
    }

    pub fn transform(self, transform: &Matrix<4>) -> Self {
        let mut bounds = Self::default();
        let corners = [
            self.minimum,
            Tuple::point(self.minimum.x, self.minimum.y, self.maximum.z),
            Tuple::point(self.minimum.x, self.maximum.y, self.minimum.z),
            Tuple::point(self.maximum.x, self.minimum.y, self.minimum.z),
            Tuple::point(self.minimum.x, self.maximum.y, self.maximum.z),
            Tuple::point(self.maximum.x, self.minimum.y, self.maximum.z),
            Tuple::point(self.maximum.x, self.maximum.y, self.minimum.z),
            self.maximum,
        ];

        for corner in corners.iter() {
            let transformed_point = *transform * *corner;
            bounds.add_point(&transformed_point);
        }

        self
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
        let tmin_numerator = min - origin;
        let tmax_numerator = max - origin;

        let tmin = tmin_numerator / direction;
        let tmax = tmax_numerator / direction;

        if tmin > tmax {
            (tmax, tmin)
        } else {
            (tmin, tmax)
        }
    }
}

impl Default for Bounds {
    fn default() -> Self {
        Self {
            minimum: Tuple::point(f64::INFINITY, f64::INFINITY, f64::INFINITY),
            maximum: Tuple::point(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // TODO: tests
}
