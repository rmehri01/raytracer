use approx::AbsDiffEq;

use crate::core::tuple::Tuple;

use super::color::Color;

/// A striped pattern that alternates between two colors.
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Stripe {
    a: Color,
    b: Color,
}

impl Stripe {
    pub fn new(a: Color, b: Color) -> Self {
        Self { a, b }
    }

    pub fn stripe_at(&self, point: &Tuple) -> Color {
        if point.x.floor() as i32 % 2 == 0 {
            self.a
        } else {
            self.b
        }
    }
}

impl AbsDiffEq for Stripe {
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
    fn stripe_new() {
        let stripe = Stripe::new(Color::WHITE, Color::BLACK);

        assert_abs_diff_eq!(stripe.a, Color::WHITE);
        assert_abs_diff_eq!(stripe.b, Color::BLACK);
    }

    #[test]
    fn stripe_constant_in_y() {
        let stripe = Stripe::new(Color::WHITE, Color::BLACK);

        assert_abs_diff_eq!(stripe.stripe_at(&Tuple::point(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.stripe_at(&Tuple::point(0.0, 1.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.stripe_at(&Tuple::point(0.0, 2.0, 0.0)), Color::WHITE);
    }

    #[test]
    fn stripe_constant_in_z() {
        let stripe = Stripe::new(Color::WHITE, Color::BLACK);

        assert_abs_diff_eq!(stripe.stripe_at(&Tuple::point(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.stripe_at(&Tuple::point(0.0, 0.0, 1.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.stripe_at(&Tuple::point(0.0, 0.0, 2.0)), Color::WHITE);
    }

    #[test]
    fn stripe_alternates_in_x() {
        let stripe = Stripe::new(Color::WHITE, Color::BLACK);

        assert_abs_diff_eq!(stripe.stripe_at(&Tuple::point(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.stripe_at(&Tuple::point(0.9, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.stripe_at(&Tuple::point(1.0, 0.0, 0.0)), Color::BLACK);
        assert_abs_diff_eq!(
            stripe.stripe_at(&Tuple::point(-0.1, 0.0, 0.0)),
            Color::BLACK
        );
        assert_abs_diff_eq!(
            stripe.stripe_at(&Tuple::point(-1.0, 0.0, 0.0)),
            Color::BLACK
        );
        assert_abs_diff_eq!(
            stripe.stripe_at(&Tuple::point(-1.1, 0.0, 0.0)),
            Color::WHITE
        );
    }
}
