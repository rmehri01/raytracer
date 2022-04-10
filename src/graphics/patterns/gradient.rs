use approx::AbsDiffEq;

use crate::{core::tuple::Tuple, graphics::color::Color};

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Gradient {
    start: Color,
    end: Color,
}

impl Gradient {
    pub fn new(start: Color, end: Color) -> Self {
        Self { start, end }
    }

    pub fn pattern_at(&self, point: &Tuple) -> Color {
        let distance = self.end - self.start;
        let fraction = point.x - point.x.floor();

        self.start + distance * fraction
    }
}

impl AbsDiffEq for Gradient {
    type Epsilon = <Color as AbsDiffEq>::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        Color::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.start.abs_diff_eq(&other.start, epsilon) && self.end.abs_diff_eq(&other.end, epsilon)
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn linearly_interpolates_between_colors() {
        let gradient = Gradient::new(Color::WHITE, Color::BLACK);

        assert_abs_diff_eq!(
            gradient.pattern_at(&Tuple::point(0.0, 0.0, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            gradient.pattern_at(&Tuple::point(0.25, 0.25, 0.25)),
            Color::new(0.75, 0.75, 0.75)
        );
        assert_abs_diff_eq!(
            gradient.pattern_at(&Tuple::point(0.5, 0.5, 0.5)),
            Color::new(0.5, 0.5, 0.5)
        );
        assert_abs_diff_eq!(
            gradient.pattern_at(&Tuple::point(0.75, 0.75, 0.75)),
            Color::new(0.25, 0.25, 0.25)
        );
    }
}
