use noise::{NoiseFn, SuperSimplex};

use crate::{
    core::{
        matrix::{Matrix, Transformation},
        point::Point,
    },
    raytracer::shapes::Primitive,
};

use super::color::Color;

#[derive(Debug, PartialEq, Clone)]
pub struct Pattern {
    transform: Transformation,
    transform_inversed: Transformation,
    kind: Kind,
}

impl Pattern {
    fn new(kind: Kind) -> Self {
        Self {
            transform: Matrix::identity(),
            transform_inversed: Matrix::identity(),
            kind,
        }
    }

    pub fn new_solid(color: Color) -> Self {
        Self::new(Kind::Solid(color))
    }

    pub fn new_blend(a: Self, b: Self) -> Self {
        Self::new_mixture(Mixture::Blend, a, b)
    }

    pub fn new_perturb(pat: Self) -> Self {
        Self::new(Kind::Perturb(Box::new(pat)))
    }

    pub fn new_stripe(a: Self, b: Self) -> Self {
        Self::new_mixture(Mixture::Stripe, a, b)
    }

    pub fn new_gradient(start: Self, end: Self) -> Self {
        Self::new_mixture(Mixture::Gradient, start, end)
    }

    pub fn new_radial_gradient(start: Self, end: Self) -> Self {
        Self::new_mixture(Mixture::RadialGradient, start, end)
    }

    pub fn new_ring(a: Self, b: Self) -> Self {
        Self::new_mixture(Mixture::Ring, a, b)
    }

    pub fn new_checker(a: Self, b: Self) -> Self {
        Self::new_mixture(Mixture::Checker, a, b)
    }

    fn new_mixture(mixture: Mixture, a: Self, b: Self) -> Self {
        Self::new(Kind::Mixture(mixture, Box::new(a), Box::new(b)))
    }

    #[cfg(test)]
    pub fn new_test() -> Self {
        Self::new(Kind::Test)
    }

    pub fn with_transform(mut self, transform: Transformation) -> Self {
        self.transform = transform;
        self.transform_inversed = transform.inverse();
        self
    }

    pub fn pattern_at_shape(
        &self,
        shape: &Primitive,
        world_point: &Point,
        trail: &im_rc::Vector<Transformation>,
    ) -> Color {
        let object_point = shape.world_to_object(world_point, trail);

        self.pattern_at(&object_point)
    }

    pub fn pattern_at(&self, object_point: &Point) -> Color {
        let pattern_point = self.transform_inversed * *object_point;

        self.kind.pattern_at(&pattern_point)
    }
}

#[derive(Debug, PartialEq, Clone)]
enum Kind {
    Solid(Color),
    Perturb(Box<Pattern>),
    Mixture(Mixture, Box<Pattern>, Box<Pattern>),
    #[cfg(test)]
    Test,
}

impl Kind {
    fn pattern_at(&self, pattern_point: &Point) -> Color {
        match self {
            Self::Solid(color) => *color,
            Self::Perturb(pat) => {
                const PERTURB_AMOUNT: f64 = 0.4;

                let simplex = SuperSimplex::new();
                let perturbed_point = Point::new(
                    pattern_point.x
                        + simplex.get([pattern_point.x, pattern_point.y, pattern_point.z])
                            * PERTURB_AMOUNT,
                    pattern_point.y
                        + simplex.get([pattern_point.x, pattern_point.y, pattern_point.z + 1.0])
                            * PERTURB_AMOUNT,
                    pattern_point.z
                        + simplex.get([pattern_point.x, pattern_point.y, pattern_point.z + 2.0])
                            * PERTURB_AMOUNT,
                );

                pat.pattern_at(&perturbed_point)
            }
            Self::Mixture(mixture, a, b) => mixture.pattern_at(
                pattern_point,
                a.pattern_at(pattern_point),
                b.pattern_at(pattern_point),
            ),
            #[cfg(test)]
            Self::Test => Color::new(pattern_point.x, pattern_point.y, pattern_point.z),
        }
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Mixture {
    Blend,
    Stripe,
    Gradient,
    RadialGradient,
    Ring,
    Checker,
}

impl Mixture {
    fn pattern_at(self, pattern_point: &Point, a: Color, b: Color) -> Color {
        match self {
            Self::Blend => (a + b) * 0.5,
            Self::Stripe => Self::stripe_at(pattern_point, a, b),
            Self::Gradient => Self::gradient_at(pattern_point, a, b),
            Self::RadialGradient => Self::radial_gradient_at(pattern_point, a, b),
            Self::Ring => Self::ring_at(pattern_point, a, b),
            Self::Checker => Self::checker_at(pattern_point, a, b),
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

    fn radial_gradient_at(point: &Point, start: Color, end: Color) -> Color {
        let distance = (*point - Point::new(0.0, 0.0, 0.0)).magnitude();
        let fraction = distance - distance.floor();

        start + (end - start) * fraction
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

    use crate::raytracer::shapes::{Primitive, SetProperties};

    use super::*;

    #[test]
    fn pattern_with_shape_transformation() {
        let shape = Primitive::new_sphere().with_transform(Matrix::scaling(2.0, 2.0, 2.0));
        let pattern = Pattern::new_test();

        assert_abs_diff_eq!(
            pattern.pattern_at_shape(&shape, &Point::new(2.0, 3.0, 4.0), &im_rc::Vector::new()),
            Color::new(1.0, 1.5, 2.0)
        );
    }

    #[test]
    fn pattern_with_pattern_transformation() {
        let shape = Primitive::new_sphere();
        let pattern = Pattern::new(Kind::Test).with_transform(Matrix::scaling(2.0, 2.0, 2.0));

        assert_abs_diff_eq!(
            pattern.pattern_at_shape(&shape, &Point::new(2.0, 3.0, 4.0), &im_rc::Vector::new()),
            Color::new(1.0, 1.5, 2.0)
        );
    }

    #[test]
    fn pattern_with_both_transformations() {
        let shape = Primitive::new_sphere().with_transform(Matrix::scaling(2.0, 2.0, 2.0));
        let pattern = Pattern::new(Kind::Test).with_transform(Matrix::translation(0.5, 1.0, 1.5));

        assert_abs_diff_eq!(
            pattern.pattern_at_shape(&shape, &Point::new(2.5, 3.0, 3.5), &im_rc::Vector::new()),
            Color::new(0.75, 0.5, 0.25)
        );
    }

    #[test]
    fn solid_pattern() {
        let color = Color::new(1.0, 0.2, 0.1);
        let pattern = Kind::Solid(color);

        assert_abs_diff_eq!(pattern.pattern_at(&Point::new(0.0, 0.0, 0.0)), color);
        assert_abs_diff_eq!(pattern.pattern_at(&Point::new(1.0, -1.0, 0.5)), color);
    }

    #[test]
    fn stripe_constant_in_y() {
        let stripe = Pattern::new_stripe(
            Pattern::new_solid(Color::WHITE),
            Pattern::new_solid(Color::BLACK),
        );

        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.0, 1.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.0, 2.0, 0.0)), Color::WHITE);
    }

    #[test]
    fn stripe_constant_in_z() {
        let stripe = Pattern::new_stripe(
            Pattern::new_solid(Color::WHITE),
            Pattern::new_solid(Color::BLACK),
        );

        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.0, 0.0, 1.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.0, 0.0, 2.0)), Color::WHITE);
    }

    #[test]
    fn stripe_alternates_in_x() {
        let stripe = Pattern::new_stripe(
            Pattern::new_solid(Color::WHITE),
            Pattern::new_solid(Color::BLACK),
        );

        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(0.9, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(1.0, 0.0, 0.0)), Color::BLACK);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(-0.1, 0.0, 0.0)), Color::BLACK);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(-1.0, 0.0, 0.0)), Color::BLACK);
        assert_abs_diff_eq!(stripe.pattern_at(&Point::new(-1.1, 0.0, 0.0)), Color::WHITE);
    }

    #[test]
    fn linearly_interpolates_between_colors() {
        let gradient = Pattern::new_gradient(
            Pattern::new_solid(Color::WHITE),
            Pattern::new_solid(Color::BLACK),
        );

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
        let ring = Pattern::new_ring(
            Pattern::new_solid(Color::WHITE),
            Pattern::new_solid(Color::BLACK),
        );

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
        let pattern = Pattern::new_checker(
            Pattern::new_solid(Color::WHITE),
            Pattern::new_solid(Color::BLACK),
        );

        assert_abs_diff_eq!(pattern.pattern_at(&Point::new(0.0, 0.0, 0.0)), Color::WHITE);
        assert_abs_diff_eq!(
            pattern.pattern_at(&Point::new(0.99, 0.0, 0.0)),
            Color::WHITE
        );
        assert_abs_diff_eq!(pattern.pattern_at(&Point::new(1.0, 0.0, 0.0)), Color::BLACK);
    }

    #[test]
    fn checkers_repeat_in_y() {
        let pattern = Pattern::new_checker(
            Pattern::new_solid(Color::WHITE),
            Pattern::new_solid(Color::BLACK),
        );

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
        let pattern = Pattern::new_checker(
            Pattern::new_solid(Color::WHITE),
            Pattern::new_solid(Color::BLACK),
        );

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
