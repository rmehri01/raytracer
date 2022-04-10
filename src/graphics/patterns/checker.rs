use approx::AbsDiffEq;

use crate::{core::tuple::Tuple, graphics::color::Color};

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Checker {
    a: Color,
    b: Color,
}

impl Checker {
    pub fn new(a: Color, b: Color) -> Self {
        Self { a, b }
    }

    pub fn pattern_at(&self, point: &Tuple) -> Color {
        let sum_floors = point.x.floor() + point.y.floor() + point.z.floor();

        if sum_floors as i32 % 2 == 0 {
            self.a
        } else {
            self.b
        }
    }
}

impl AbsDiffEq for Checker {
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
    fn checkers_repeat_in_x() {
        let pattern = Checker::new(Color::WHITE, Color::BLACK);

        assert_abs_diff_eq!(
            pattern.pattern_at(&Tuple::point(0.0, 0.0, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            pattern.pattern_at(&Tuple::point(0.99, 0.0, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            pattern.pattern_at(&Tuple::point(1.0, 0.0, 0.0)),
            Color::BLACK
        );
    }

    #[test]
    fn checkers_repeat_in_y() {
        let pattern = Checker::new(Color::WHITE, Color::BLACK);

        assert_abs_diff_eq!(
            pattern.pattern_at(&Tuple::point(0.0, 0.0, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            pattern.pattern_at(&Tuple::point(0.0, 0.99, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            pattern.pattern_at(&Tuple::point(0.0, 1.01, 0.0)),
            Color::BLACK
        );
    }

    #[test]
    fn checkers_repeat_in_z() {
        let pattern = Checker::new(Color::WHITE, Color::BLACK);

        assert_abs_diff_eq!(
            pattern.pattern_at(&Tuple::point(0.0, 0.0, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            pattern.pattern_at(&Tuple::point(0.0, 0.0, 0.99)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            pattern.pattern_at(&Tuple::point(0.0, 0.0, 1.01)),
            Color::BLACK
        );
    }
}
