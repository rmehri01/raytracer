use crate::{core::tuple::Tuple, graphics::color::Color};

/// A light source existing in a single point in space.
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct PointLight {
    pub position: Tuple,
    pub intensity: Color,
}

impl PointLight {
    pub fn new(position: Tuple, intensity: Color) -> Self {
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
        let light = PointLight::new(Tuple::point(0.0, 0.0, 0.0), Color::WHITE);

        assert_abs_diff_eq!(light.position, Tuple::point(0.0, 0.0, 0.0));
        assert_abs_diff_eq!(light.intensity, Color::WHITE);
    }
}
