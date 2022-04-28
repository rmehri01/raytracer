use approx::AbsDiffEq;

use crate::graphics::{color::Color, pattern::Pattern};

#[derive(Debug, Clone)]
pub struct Material {
    pub color: Color,
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64,
    pub reflective: f64,
    pub transparency: f64,
    pub refractive_index: f64,
    pub pattern: Option<Pattern>,
}

impl Default for Material {
    fn default() -> Self {
        Self {
            color: Color::WHITE,
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.0,
            reflective: 0.0,
            transparency: 0.0,
            refractive_index: 1.0,
            pattern: None,
        }
    }
}

impl PartialEq for Material {
    fn eq(&self, other: &Self) -> bool {
        self.abs_diff_eq(other, Self::default_epsilon())
    }
}

impl AbsDiffEq for Material {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-5
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.color == other.color
            && self.ambient.abs_diff_eq(&other.ambient, epsilon)
            && self.diffuse.abs_diff_eq(&other.diffuse, epsilon)
            && self.specular.abs_diff_eq(&other.specular, epsilon)
            && self.shininess.abs_diff_eq(&other.shininess, epsilon)
            && self.reflective.abs_diff_eq(&other.reflective, epsilon)
            && self.transparency.abs_diff_eq(&other.transparency, epsilon)
            && self
                .refractive_index
                .abs_diff_eq(&other.refractive_index, epsilon)
            && self.pattern == other.pattern
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn default_material() {
        let material = Material::default();

        assert_abs_diff_eq!(
            material,
            Material {
                color: Color::WHITE,
                ambient: 0.1,
                diffuse: 0.9,
                specular: 0.9,
                shininess: 200.0,
                reflective: 0.0,
                transparency: 0.0,
                refractive_index: 1.0,
                pattern: None,
            }
        );
    }
}
