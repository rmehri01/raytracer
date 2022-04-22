use crate::core::{
    matrix::{Matrix, Transformation},
    point::Point,
};

use super::color::Color;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Pattern {
    transform: Transformation,
    transform_inversed: Transformation,
    kind: PatternKind,
}

impl Pattern {
    fn new(kind: PatternKind) -> Self {
        Self {
            transform: Matrix::identity(),
            transform_inversed: Matrix::identity(),
            kind,
        }
    }

    pub fn new_stripe(a: Color, b: Color) -> Self {
        Self::new(PatternKind::Stripe { a, b })
    }

    pub fn new_gradient(start: Color, end: Color) -> Self {
        Self::new(PatternKind::Gradient { start, end })
    }

    pub fn new_ring(a: Color, b: Color) -> Self {
        Self::new(PatternKind::Ring { a, b })
    }

    pub fn new_checker(a: Color, b: Color) -> Self {
        Self::new(PatternKind::Checker { a, b })
    }

    #[cfg(test)]
    pub fn new_test() -> Self {
        Self::new(PatternKind::Test)
    }

    pub fn with_transform(mut self, transform: Transformation) -> Self {
        self.transform = transform;
        self.transform_inversed = transform.inverse();
        self
    }

    pub fn pattern_at_object_point(&self, object_point: &Point) -> Color {
        let pattern_point = self.transform_inversed * *object_point;

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
    fn pattern_at(&self, pattern_point: &Point) -> Color {
        match self {
            Self::Stripe { a, b } => Self::stripe_at(pattern_point, *a, *b),
            Self::Gradient { start, end } => Self::gradient_at(pattern_point, *start, *end),
            Self::Ring { a, b } => Self::ring_at(pattern_point, *a, *b),
            Self::Checker { a, b } => Self::checker_at(pattern_point, *a, *b),
            #[cfg(test)]
            Self::Test => Color::new(pattern_point.x, pattern_point.y, pattern_point.z),
        }
    }

    fn stripe_at(point: &Point, a: Color, b: Color) -> Color {
        if point.x.floor() as i32 % 2 == 0 {
            a
        } else {
            b
        }
    }

    fn gradient_at(point: &Point, start: Color, end: Color) -> Color {
        let distance = end - start;
        let fraction = point.x - point.x.floor();

        start + distance * fraction
    }

    fn ring_at(point: &Point, a: Color, b: Color) -> Color {
        let distance = point.x.hypot(point.z);

        if distance.floor() as i32 % 2 == 0 {
            a
        } else {
            b
        }
    }

    fn checker_at(point: &Point, a: Color, b: Color) -> Color {
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

    use crate::raytracer::shape::Shape;

    use super::*;

    #[test]
    fn pattern_with_shape_transformation() {
        let shape = Shape::new_sphere().with_transform(Matrix::scaling(2.0, 2.0, 2.0));
        let pattern = Pattern::new_test();

        let object_point = shape.world_to_object(&Point::new(2.0, 3.0, 4.0), &im_rc::Vector::new());

        assert_abs_diff_eq!(
            pattern.pattern_at_object_point(&object_point),
            Color::new(1.0, 1.5, 2.0)
        );
    }

    #[test]
    fn pattern_with_pattern_transformation() {
        let shape = Shape::new_sphere();
        let pattern =
            Pattern::new(PatternKind::Test).with_transform(Matrix::scaling(2.0, 2.0, 2.0));

        let object_point = shape.world_to_object(&Point::new(2.0, 3.0, 4.0), &im_rc::Vector::new());

        assert_abs_diff_eq!(
            pattern.pattern_at_object_point(&object_point),
            Color::new(1.0, 1.5, 2.0)
        );
    }

    #[test]
    fn pattern_with_both_transformations() {
        let shape = Shape::new_sphere().with_transform(Matrix::scaling(2.0, 2.0, 2.0));
        let pattern =
            Pattern::new(PatternKind::Test).with_transform(Matrix::translation(0.5, 1.0, 1.5));

        let object_point = shape.world_to_object(&Point::new(2.5, 3.0, 3.5), &im_rc::Vector::new());

        assert_abs_diff_eq!(
            pattern.pattern_at_object_point(&object_point),
            Color::new(0.75, 0.5, 0.25)
        );
    }

    #[test]
    fn stripe_constant_in_y() {
        let stripe = PatternKind::Stripe {
            a: Color::WHITE,
            b: Color::BLACK,
        };

        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.0, 1.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.0, 2.0, 0.0)), Color::WHITE);
    }

    #[test]
    fn stripe_constant_in_z() {
        let stripe = PatternKind::Stripe {
            a: Color::WHITE,
            b: Color::BLACK,
        };

        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.0, 0.0, 1.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.0, 0.0, 2.0)), Color::WHITE);
    }

    #[test]
    fn stripe_alternates_in_x() {
        let stripe = PatternKind::Stripe {
            a: Color::WHITE,
            b: Color::BLACK,
        };

        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.9, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(1.0, 0.0, 0.0)), Color::BLACK);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(-0.1, 0.0, 0.0)), Color::BLACK);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(-1.0, 0.0, 0.0)), Color::BLACK);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(-1.1, 0.0, 0.0)), Color::WHITE);
    }

    #[test]
    fn linearly_interpolates_between_colors() {
        let gradient = PatternKind::Gradient {
            start: Color::WHITE,
            end: Color::BLACK,
        };

        assert_abs_diff_eq!(
            gradient.pattern_at(&Point::new(0.0, 0.0, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            gradient.pattern_at(&Point::new(0.25, 0.25, 0.25)),
            Color::new(0.75, 0.75, 0.75)
        );
        assert_abs_diff_eq!(
            gradient.pattern_at(&Point::new(0.5, 0.5, 0.5)),
            Color::new(0.5, 0.5, 0.5)
        );
        assert_abs_diff_eq!(
            gradient.pattern_at(&Point::new(0.75, 0.75, 0.75)),
            Color::new(0.25, 0.25, 0.25)
        );
    }

    #[test]
    fn ring_extends_in_x_and_z() {
        let ring = PatternKind::Ring {
            a: Color::WHITE,
            b: Color::BLACK,
        };

        assert_abs_diff_eq!(ring.pattern_at(&Point::new(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(ring.pattern_at(&Point::new(1.0, 0.0, 0.0)), Color::BLACK);
        assert_abs_diff_eq!(ring.pattern_at(&Point::new(0.0, 0.0, 1.0)), Color::BLACK);
        assert_abs_diff_eq!(
            ring.pattern_at(&Point::new(0.708, 0.0, 0.708)),
            Color::BLACK
        );
    }

    #[test]
    fn checkers_repeat_in_x() {
        let pattern = PatternKind::Checker {
            a: Color::WHITE,
            b: Color::BLACK,
        };

        assert_abs_diff_eq!(pattern.pattern_at(&Point::new(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(
            pattern.pattern_at(&Point::new(0.99, 0.0, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(pattern.pattern_at(&Point::new(1.0, 0.0, 0.0)), Color::BLACK);
    }

    #[test]
    fn checkers_repeat_in_y() {
        let pattern = PatternKind::Checker {
            a: Color::WHITE,
            b: Color::BLACK,
        };

        assert_abs_diff_eq!(pattern.pattern_at(&Point::new(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(
            pattern.pattern_at(&Point::new(0.0, 0.99, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            pattern.pattern_at(&Point::new(0.0, 1.01, 0.0)),
            Color::BLACK
        );
    }

    #[test]
    fn checkers_repeat_in_z() {
        let pattern = PatternKind::Checker {
            a: Color::WHITE,
            b: Color::BLACK,
        };

        assert_abs_diff_eq!(pattern.pattern_at(&Point::new(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(
            pattern.pattern_at(&Point::new(0.0, 0.0, 0.99)),
            Color::WHITE
        );
        assert_abs_diff_eq!(
            pattern.pattern_at(&Point::new(0.0, 0.0, 1.01)),
            Color::BLACK
        );
    }
}
