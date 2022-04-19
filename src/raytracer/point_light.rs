use crate::{core::point::Point, graphics::color::Color};

/// A light source existing in a single point in space.
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct PointLight {
    pub position: Point,
    pub intensity: Color,
}

impl PointLight {
    pub fn new(position: Point, intensity: Color) -> Self {
        Self {
            position,
            intensity,
        }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn test_point_light_position() {
        let light = PointLight::new(Point::new(0.0, 0.0, 0.0), Color::WHITE);

        assert_abs_diff_eq!(light.position, Point::new(0.0, 0.0, 0.0));
        assert_abs_diff_eq!(light.intensity, Color::WHITE);
    }
}
