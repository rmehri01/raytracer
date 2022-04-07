use approx::AbsDiffEq;

use crate::{core::tuple::Tuple, graphics::color::Color};

/// A light source existing in a single point in space.
#[derive(Debug, PartialEq)]
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

impl AbsDiffEq for PointLight {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-5
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.position.abs_diff_eq(&other.position, epsilon)
            && self.intensity.abs_diff_eq(&other.intensity, epsilon)
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn test_point_light_position() {
        let light = PointLight::new(Tuple::point(0.0, 0.0, 0.0), Color::new(1.0, 1.0, 1.0));

        assert_abs_diff_eq!(light.position, Tuple::point(0.0, 0.0, 0.0));
        assert_abs_diff_eq!(light.intensity, Color::new(1.0, 1.0, 1.0));
    }
}
