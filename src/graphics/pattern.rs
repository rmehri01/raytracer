use approx::AbsDiffEq;

use crate::{
    core::{matrix::Matrix, tuple::Tuple},
    raytracer::object::Object,
};

use super::{color::Color, patterns::stripe::Stripe};

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Pattern {
    transform: Matrix<4>,
    kind: PatternKind,
}

impl Pattern {
    fn new(transform: Matrix<4>, kind: PatternKind) -> Self {
        Self { transform, kind }
    }

    pub fn new_stripe(a: Color, b: Color) -> Self {
        Self::new(Matrix::identity(), PatternKind::Stripe(Stripe::new(a, b)))
    }

    pub fn pattern_at_object(&self, object: &Object, world_point: &Tuple) -> Color {
        let object_point = object.transform.inverse() * *world_point;
        let pattern_point = self.transform.inverse() * object_point;

        self.kind.pattern_at(&pattern_point)
    }
}

impl AbsDiffEq for Pattern {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-5
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.transform.abs_diff_eq(&other.transform, epsilon)
            && self.kind.abs_diff_eq(&other.kind, epsilon)
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
enum PatternKind {
    Stripe(Stripe),
}

impl PatternKind {
    fn pattern_at(&self, pattern_point: &Tuple) -> Color {
        match self {
            PatternKind::Stripe(stripe) => stripe.pattern_at(pattern_point),
        }
    }
}

impl AbsDiffEq for PatternKind {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-5
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        match (self, other) {
            (PatternKind::Stripe(a), PatternKind::Stripe(b)) => a.abs_diff_eq(b, epsilon),
        }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn stripes_with_object_transformation() {
        let mut object = Object::new_sphere();
        object.transform = Matrix::scaling(2.0, 2.0, 2.0);

        let pattern = Pattern::new_stripe(Color::WHITE, Color::BLACK);

        assert_abs_diff_eq!(
            pattern.pattern_at_object(&object, &Tuple::point(1.5, 0.0, 0.0)),
            Color::WHITE
        );
    }

    #[test]
    fn stripes_with_pattern_transformation() {
        let object = Object::new_sphere();

        let mut pattern = Pattern::new_stripe(Color::WHITE, Color::BLACK);
        pattern.transform = Matrix::scaling(2.0, 2.0, 2.0);

        assert_abs_diff_eq!(
            pattern.pattern_at_object(&object, &Tuple::point(1.5, 0.0, 0.0)),
            Color::WHITE
        );
    }

    #[test]
    fn stripes_with_both_transformations() {
        let mut object = Object::new_sphere();
        object.transform = Matrix::scaling(2.0, 2.0, 2.0);

        let mut pattern = Pattern::new_stripe(Color::WHITE, Color::BLACK);
        pattern.transform = Matrix::translation(0.5, 0.0, 0.0);

        assert_abs_diff_eq!(
            pattern.pattern_at_object(&object, &Tuple::point(2.5, 0.0, 0.0)),
            Color::WHITE
        );
    }
}
