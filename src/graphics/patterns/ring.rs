use approx::AbsDiffEq;

use crate::{core::tuple::Tuple, graphics::color::Color};

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Ring {
    a: Color,
    b: Color,
}

impl Ring {
    pub fn new(a: Color, b: Color) -> Self {
        Self { a, b }
    }

    pub fn pattern_at(&self, point: &Tuple) -> Color {
        let distance = (point.x * point.x + point.z * point.z).sqrt();
        if distance.floor() as i32 % 2 == 0 {
            self.a
        } else {
            self.b
        }
    }
}

impl AbsDiffEq for Ring {
    type Epsilon = <Color as AbsDiffEq>::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        Color::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.a.abs_diff_eq(&other.a, epsilon) && self.b.abs_diff_eq(&other.b, epsilon)
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn ring_extends_in_x_and_z() {
        let ring = Ring::new(Color::WHITE, Color::BLACK);

        assert_abs_diff_eq!(ring.pattern_at(&Tuple::point(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(ring.pattern_at(&Tuple::point(1.0, 0.0, 0.0)), Color::BLACK);
        assert_abs_diff_eq!(ring.pattern_at(&Tuple::point(0.0, 0.0, 1.0)), Color::BLACK);
        assert_abs_diff_eq!(
            ring.pattern_at(&Tuple::point(0.708, 0.0, 0.708)),
            Color::BLACK
        );
    }
}