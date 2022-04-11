use std::ops;

use approx::AbsDiffEq;

#[derive(Debug, Clone, Copy)]
pub struct Color {
    pub red: f64,
    pub green: f64,
    pub blue: f64,
}

impl Color {
    pub fn new(red: f64, green: f64, blue: f64) -> Self {
        Self { red, green, blue }
    }

    pub const BLACK: Self = Self {
        red: 0.0,
        green: 0.0,
        blue: 0.0,
    };

    pub const WHITE: Self = Self {
        red: 1.0,
        green: 1.0,
        blue: 1.0,
    };
}

impl ops::Add for Color {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            red: self.red + other.red,
            green: self.green + other.green,
            blue: self.blue + other.blue,
        }
    }
}

impl ops::Sub for Color {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            red: self.red - other.red,
            green: self.green - other.green,
            blue: self.blue - other.blue,
        }
    }
}

impl ops::Mul<f64> for Color {
    type Output = Self;

    fn mul(self, other: f64) -> Self {
        Self {
            red: self.red * other,
            green: self.green * other,
            blue: self.blue * other,
        }
    }
}

impl ops::Mul for Color {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self {
            red: self.red * other.red,
            green: self.green * other.green,
            blue: self.blue * other.blue,
        }
    }
}

impl PartialEq for Color {
    fn eq(&self, other: &Self) -> bool {
        self.abs_diff_eq(other, Self::default_epsilon())
    }
}

impl AbsDiffEq for Color {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-4
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.red.abs_diff_eq(&other.red, epsilon)
            && self.green.abs_diff_eq(&other.green, epsilon)
            && self.blue.abs_diff_eq(&other.blue, epsilon)
    }
}

#[cfg(test)]
mod tests {
    use approx::{assert_abs_diff_eq, assert_relative_eq};

    use super::*;

    #[test]
    fn color_new() {
        let color = Color::new(-0.5, 0.4, 1.7);

        assert_relative_eq!(color.red, -0.5);
        assert_relative_eq!(color.green, 0.4);
        assert_relative_eq!(color.blue, 1.7);
    }

    #[test]
    fn color_add() {
        let color1 = Color::new(0.9, 0.6, 0.75);
        let color2 = Color::new(0.7, 0.1, 0.25);

        assert_abs_diff_eq!(color1 + color2, Color::new(1.6, 0.7, 1.0));
    }

    #[test]
    fn color_sub() {
        let color1 = Color::new(0.9, 0.6, 0.75);
        let color2 = Color::new(0.7, 0.1, 0.25);

        assert_abs_diff_eq!(color1 - color2, Color::new(0.2, 0.5, 0.5));
    }

    #[test]
    fn color_mul_scalar() {
        let color1 = Color::new(0.2, 0.3, 0.4);

        assert_abs_diff_eq!(color1 * 2.0, Color::new(0.4, 0.6, 0.8));
    }

    #[test]
    fn color_mul_color() {
        let color1 = Color::new(1.0, 0.2, 0.4);
        let color2 = Color::new(0.9, 1.0, 0.1);

        assert_abs_diff_eq!(color1 * color2, Color::new(0.9, 0.2, 0.04));
    }
}
