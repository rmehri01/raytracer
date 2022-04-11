use crate::{
    core::{matrix::Matrix, tuple::Tuple},
    raytracer::object::Object,
};

use super::{
    color::Color,
    patterns::{checker::Checker, gradient::Gradient, ring::Ring, stripe::Stripe},
};

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
        Self::new(Matrix::identity(), PatternKind::Stripe(Stripe::new(a, b)))
    }

    pub fn new_gradient(start: Color, end: Color) -> Self {
        Self::new(
            Matrix::identity(),
            PatternKind::Gradient(Gradient::new(start, end)),
        )
    }

    pub fn new_ring(a: Color, b: Color) -> Self {
        Self::new(Matrix::identity(), PatternKind::Ring(Ring::new(a, b)))
    }

    pub fn new_checker(a: Color, b: Color) -> Self {
        Self::new(Matrix::identity(), PatternKind::Checker(Checker::new(a, b)))
    }

    pub fn pattern_at_object(&self, object: &Object, world_point: &Tuple) -> Color {
        let object_point = object.transform.inverse() * *world_point;
        let pattern_point = self.transform.inverse() * object_point;

        self.kind.pattern_at(&pattern_point)
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
enum PatternKind {
    Stripe(Stripe),
    Gradient(Gradient),
    Ring(Ring),
    Checker(Checker),
}

impl PatternKind {
    fn pattern_at(&self, pattern_point: &Tuple) -> Color {
        match self {
            PatternKind::Stripe(stripe) => stripe.pattern_at(pattern_point),
            PatternKind::Gradient(gradient) => gradient.pattern_at(pattern_point),
            PatternKind::Ring(ring) => ring.pattern_at(pattern_point),
            PatternKind::Checker(checker) => checker.pattern_at(pattern_point),
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
