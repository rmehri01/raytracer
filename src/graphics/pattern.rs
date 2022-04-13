use crate::{
    core::{matrix::Matrix, tuple::Tuple},
    raytracer::shape::Shape,
};

use super::color::Color;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Pattern {
    pub transform: Matrix<4>,
    kind: PatternKind,
}

impl Pattern {
    fn new(transform: Matrix<4>, kind: PatternKind) -> Self {
        Self { transform, kind }
    }

    pub fn new_stripe(a: Color, b: Color) -> Self {
        Self::new(Matrix::identity(), PatternKind::Stripe { a, b })
    }

    pub fn new_gradient(start: Color, end: Color) -> Self {
        Self::new(Matrix::identity(), PatternKind::Gradient { start, end })
    }

    pub fn new_ring(a: Color, b: Color) -> Self {
        Self::new(Matrix::identity(), PatternKind::Ring { a, b })
    }

    pub fn new_checker(a: Color, b: Color) -> Self {
        Self::new(Matrix::identity(), PatternKind::Checker { a, b })
    }

    #[cfg(test)]
    pub fn new_test() -> Self {
        Self::new(Matrix::identity(), PatternKind::Test)
    }

    pub fn pattern_at_shape(&self, shape: &Shape, world_point: &Tuple) -> Color {
        let object_point = shape.transform.inverse() * *world_point;
        let pattern_point = self.transform.inverse() * object_point;

        self.kind.pattern_at(&pattern_point)
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
enum PatternKind {
    Stripe {
        a: Color,
        b: Color,
    },
    Gradient {
        start: Color,
        end: Color,
    },
    Ring {
        a: Color,
        b: Color,
    },
    Checker {
        a: Color,
        b: Color,
    },
    #[cfg(test)]
    Test,
}

impl PatternKind {
    fn pattern_at(&self, pattern_point: &Tuple) -> Color {
        match self {
            Self::Stripe { a, b } => Self::stripe_at(*a, *b, pattern_point),
            Self::Gradient { start, end } => Self::gradient_at(*start, *end, pattern_point),
            Self::Ring { a, b } => Self::ring_at(*a, *b, pattern_point),
            Self::Checker { a, b } => Self::checker_at(*a, *b, pattern_point),
            #[cfg(test)]
            Self::Test => Color::new(pattern_point.x, pattern_point.y, pattern_point.z),
        }
    }

    fn stripe_at(a: Color, b: Color, point: &Tuple) -> Color {
        if point.x.floor() as i32 % 2 == 0 {
            a
        } else {
            b
        }
    }

    fn gradient_at(start: Color, end: Color, point: &Tuple) -> Color {
        let distance = end - start;
        let fraction = point.x - point.x.floor();

        start + distance * fraction
    }

    fn ring_at(a: Color, b: Color, point: &Tuple) -> Color {
        let distance = (point.x * point.x + point.z * point.z).sqrt();

        if distance.floor() as i32 % 2 == 0 {
            a
        } else {
            b
        }
    }

    fn checker_at(a: Color, b: Color, point: &Tuple) -> Color {
        let sum_floors = point.x.floor() + point.y.floor() + point.z.floor();

        if sum_floors as i32 % 2 == 0 {
            a
        } else {
            b
        }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn pattern_with_shape_transformation() {
        let mut shape = Shape::new_sphere();
        shape.transform = Matrix::scaling(2.0, 2.0, 2.0);

        let pattern = Pattern::new_test();

        assert_abs_diff_eq!(
            pattern.pattern_at_shape(&shape, &Tuple::point(2.0, 3.0, 4.0)),
            Color::new(1.0, 1.5, 2.0)
        );
    }

    #[test]
    fn pattern_with_pattern_transformation() {
        let shape = Shape::new_sphere();
        let pattern = Pattern::new(Matrix::scaling(2.0, 2.0, 2.0), PatternKind::Test);

        assert_abs_diff_eq!(
            pattern.pattern_at_shape(&shape, &Tuple::point(2.0, 3.0, 4.0)),
            Color::new(1.0, 1.5, 2.0)
        );
    }

    #[test]
    fn pattern_with_both_transformations() {
        let mut shape = Shape::new_sphere();
        shape.transform = Matrix::scaling(2.0, 2.0, 2.0);

        let pattern = Pattern::new(Matrix::translation(0.5, 1.0, 1.5), PatternKind::Test);

        assert_abs_diff_eq!(
            pattern.pattern_at_shape(&shape, &Tuple::point(2.5, 3.0, 3.5)),
            Color::new(0.75, 0.5, 0.25)
        );
    }

    #[test]
    fn stripe_constant_in_y() {
        let stripe = PatternKind::Stripe {
            a: Color::WHITE,
            b: Color::BLACK,
        };

        assert_abs_diff_eq!(
            stripe.pattern_at(&Tuple::point(0.0, 0.0, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            stripe.pattern_at(&Tuple::point(0.0, 1.0, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            stripe.pattern_at(&Tuple::point(0.0, 2.0, 0.0)),
            Color::WHITE
        );
    }

    #[test]
    fn stripe_constant_in_z() {
        let stripe = PatternKind::Stripe {
            a: Color::WHITE,
            b: Color::BLACK,
        };

        assert_abs_diff_eq!(
            stripe.pattern_at(&Tuple::point(0.0, 0.0, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            stripe.pattern_at(&Tuple::point(0.0, 0.0, 1.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            stripe.pattern_at(&Tuple::point(0.0, 0.0, 2.0)),
            Color::WHITE
        );
    }

    #[test]
    fn stripe_alternates_in_x() {
        let stripe = PatternKind::Stripe {
            a: Color::WHITE,
            b: Color::BLACK,
        };

        assert_abs_diff_eq!(
            stripe.pattern_at(&Tuple::point(0.0, 0.0, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            stripe.pattern_at(&Tuple::point(0.9, 0.0, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            stripe.pattern_at(&Tuple::point(1.0, 0.0, 0.0)),
            Color::BLACK
        );
        assert_abs_diff_eq!(
            stripe.pattern_at(&Tuple::point(-0.1, 0.0, 0.0)),
            Color::BLACK
        );
        assert_abs_diff_eq!(
            stripe.pattern_at(&Tuple::point(-1.0, 0.0, 0.0)),
            Color::BLACK
        );
        assert_abs_diff_eq!(
            stripe.pattern_at(&Tuple::point(-1.1, 0.0, 0.0)),
            Color::WHITE
        );
    }

    #[test]
    fn linearly_interpolates_between_colors() {
        let gradient = PatternKind::Gradient {
            start: Color::WHITE,
            end: Color::BLACK,
        };

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

    #[test]
    fn ring_extends_in_x_and_z() {
        let ring = PatternKind::Ring {
            a: Color::WHITE,
            b: Color::BLACK,
        };

        assert_abs_diff_eq!(ring.pattern_at(&Tuple::point(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(ring.pattern_at(&Tuple::point(1.0, 0.0, 0.0)), Color::BLACK);
        assert_abs_diff_eq!(ring.pattern_at(&Tuple::point(0.0, 0.0, 1.0)), Color::BLACK);
        assert_abs_diff_eq!(
            ring.pattern_at(&Tuple::point(0.708, 0.0, 0.708)),
            Color::BLACK
        );
    }

    #[test]
    fn checkers_repeat_in_x() {
        let pattern = PatternKind::Checker {
            a: Color::WHITE,
            b: Color::BLACK,
        };

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
        let pattern = PatternKind::Checker {
            a: Color::WHITE,
            b: Color::BLACK,
        };

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
        let pattern = PatternKind::Checker {
            a: Color::WHITE,
            b: Color::BLACK,
        };

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
