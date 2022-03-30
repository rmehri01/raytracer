use std::ops;

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

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn test_color_new() {
        let color = Color::new(-0.5, 0.4, 1.7);
        assert_relative_eq!(color.red, -0.5);
        assert_relative_eq!(color.green, 0.4);
        assert_relative_eq!(color.blue, 1.7);
    }

    #[test]
    fn test_color_add() {
        let color1 = Color::new(0.9, 0.6, 0.75);
        let color2 = Color::new(0.7, 0.1, 0.25);

        let color3 = color1 + color2;

        assert_relative_eq!(color3.red, 1.6);
        assert_relative_eq!(color3.green, 0.7);
        assert_relative_eq!(color3.blue, 1.0);
    }

    #[test]
    fn test_color_sub() {
        let color1 = Color::new(0.9, 0.6, 0.75);
        let color2 = Color::new(0.7, 0.1, 0.25);

        let color3 = color1 - color2;

        assert_relative_eq!(color3.red, 0.2);
        assert_relative_eq!(color3.green, 0.5);
        assert_relative_eq!(color3.blue, 0.5);
    }

    #[test]
    fn test_color_mul_scalar() {
        let color1 = Color::new(0.2, 0.3, 0.4);

        let color2 = color1 * 2.0;

        assert_relative_eq!(color2.red, 0.4);
        assert_relative_eq!(color2.green, 0.6);
        assert_relative_eq!(color2.blue, 0.8);
    }

    #[test]
    fn test_color_mul_color() {
        let color1 = Color::new(1.0, 0.2, 0.4);
        let color2 = Color::new(0.9, 1.0, 0.1);

        let color3 = color1 * color2;

        assert_relative_eq!(color3.red, 0.9);
        assert_relative_eq!(color3.green, 0.2);
        assert_relative_eq!(color3.blue, 0.04);
    }
}
