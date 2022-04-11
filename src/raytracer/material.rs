use approx::AbsDiffEq;

use crate::{
    core::tuple::Tuple,
    graphics::{color::Color, pattern::Pattern},
};

use super::{object::Object, point_light::PointLight};

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Material {
    pub color: Color,
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64,
    pub reflective: f64,
    pub pattern: Option<Pattern>,
}

impl Material {
    /// Returns the color of the material using the Phuong Reflection Model.
    pub fn lighting(
        &self,
        object: &Object,
        light: &PointLight,
        position: &Tuple,
        eye_v: &Tuple,
        normal_v: &Tuple,
        in_shadow: bool,
    ) -> Color {
        let color = self.pattern.map_or(self.color, |pattern| {
            pattern.pattern_at_object(object, position)
        });

        let effective_color = color * light.intensity;
        let light_v = (light.position - *position).normalize();
        let ambient = effective_color * self.ambient;

        let diffuse;
        let specular;

        let light_dot_normal = light_v.dot(normal_v);
        if light_dot_normal < 0.0 || in_shadow {
            diffuse = Color::BLACK;
            specular = Color::BLACK;
        } else {
            diffuse = effective_color * self.diffuse * light_dot_normal;

            let reflect_v = (-light_v).reflect(normal_v);
            let reflect_dot_eye = reflect_v.dot(eye_v);

            specular = if reflect_dot_eye <= 0.0 {
                Color::BLACK
            } else {
                let factor = reflect_dot_eye.powf(self.shininess);
                light.intensity * self.specular * factor
            };
        }

        ambient + diffuse + specular
    }
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
            pattern: None,
        }
    }
}

impl AbsDiffEq for Material {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-5
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.color.abs_diff_eq(&other.color, epsilon)
            && self.ambient.abs_diff_eq(&other.ambient, epsilon)
            && self.diffuse.abs_diff_eq(&other.diffuse, epsilon)
            && self.specular.abs_diff_eq(&other.specular, epsilon)
            && self.shininess.abs_diff_eq(&other.shininess, epsilon)
            && self
                .pattern
                .and_then(|p1| other.pattern.map(|p2| p1.abs_diff_eq(&p2, epsilon)))
                .unwrap_or(true)
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
                pattern: None,
            }
        );
    }

    #[test]
    fn lighting_with_eye_between_light_and_surface() {
        let material = Material::default();
        let position = Tuple::point(0.0, 0.0, 0.0);

        let eye_v = Tuple::vector(0.0, 0.0, -1.0);
        let normal_v = Tuple::vector(0.0, 0.0, -1.0);
        let light = PointLight::new(Tuple::point(0.0, 0.0, -10.0), Color::WHITE);
        let object = Object::new_sphere();

        let result = material.lighting(&object, &light, &position, &eye_v, &normal_v, false);

        assert_abs_diff_eq!(result, Color::new(1.9, 1.9, 1.9));
    }

    #[test]
    fn lighting_with_eye_between_light_and_surface_eye_offset_45_degrees() {
        let material = Material::default();
        let position = Tuple::point(0.0, 0.0, 0.0);

        let eye_v = Tuple::vector(0.0, 2_f64.sqrt() / 2.0, -(2_f64.sqrt()) / 2.0);
        let normal_v = Tuple::vector(0.0, 0.0, -1.0);
        let light = PointLight::new(Tuple::point(0.0, 0.0, -10.0), Color::WHITE);
        let object = Object::new_sphere();

        let result = material.lighting(&object, &light, &position, &eye_v, &normal_v, false);

        assert_abs_diff_eq!(result, Color::WHITE);
    }

    #[test]
    fn lighting_with_eye_opposite_surface_light_offset_45_degrees() {
        let material = Material::default();
        let position = Tuple::point(0.0, 0.0, 0.0);

        let eye_v = Tuple::vector(0.0, 0.0, -1.0);
        let normal_v = Tuple::vector(0.0, 0.0, -1.0);
        let light = PointLight::new(Tuple::point(0.0, 10.0, -10.0), Color::WHITE);
        let object = Object::new_sphere();

        let result = material.lighting(&object, &light, &position, &eye_v, &normal_v, false);

        assert_abs_diff_eq!(result, Color::new(0.7364, 0.7364, 0.7364));
    }

    #[test]
    fn lighting_with_eye_in_path_of_reflection_vector() {
        let material = Material::default();
        let position = Tuple::point(0.0, 0.0, 0.0);

        let eye_v = Tuple::vector(0.0, -(2_f64.sqrt()) / 2.0, -(2_f64.sqrt()) / 2.0);
        let normal_v = Tuple::vector(0.0, 0.0, -1.0);
        let light = PointLight::new(Tuple::point(0.0, 10.0, -10.0), Color::WHITE);
        let object = Object::new_sphere();

        let result = material.lighting(&object, &light, &position, &eye_v, &normal_v, false);

        assert_abs_diff_eq!(result, Color::new(1.6364, 1.6364, 1.6364));
    }

    #[test]
    fn lighting_with_light_behind_surface() {
        let material = Material::default();
        let position = Tuple::point(0.0, 0.0, 0.0);

        let eye_v = Tuple::vector(0.0, 0.0, -1.0);
        let normal_v = Tuple::vector(0.0, 0.0, -1.0);
        let light = PointLight::new(Tuple::point(0.0, 0.0, 10.0), Color::WHITE);
        let object = Object::new_sphere();

        let result = material.lighting(&object, &light, &position, &eye_v, &normal_v, false);

        assert_abs_diff_eq!(result, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn lighting_with_surface_in_shadow() {
        let material = Material::default();
        let position = Tuple::point(0.0, 0.0, 0.0);

        let eye_v = Tuple::vector(0.0, 0.0, -1.0);
        let normal_v = Tuple::vector(0.0, 0.0, -1.0);
        let light = PointLight::new(Tuple::point(0.0, 0.0, -10.0), Color::WHITE);
        let object = Object::new_sphere();

        assert_abs_diff_eq!(
            material.lighting(&object, &light, &position, &eye_v, &normal_v, true),
            Color::new(0.1, 0.1, 0.1)
        );
    }

    #[test]
    fn lighting_with_pattern_applied() {
        let material = Material {
            ambient: 1.0,
            diffuse: 0.0,
            specular: 0.0,
            pattern: Some(Pattern::new_stripe(Color::WHITE, Color::BLACK)),
            ..Material::default()
        };

        let eye_v = Tuple::vector(0.0, 0.0, -1.0);
        let normal_v = Tuple::vector(0.0, 0.0, -1.0);
        let light = PointLight::new(Tuple::point(0.0, 0.0, -10.0), Color::WHITE);
        let object = Object::new_sphere();

        let c1 = material.lighting(
            &object,
            &light,
            &Tuple::point(0.9, 0.0, 0.0),
            &eye_v,
            &normal_v,
            false,
        );
        let c2 = material.lighting(
            &object,
            &light,
            &Tuple::point(1.1, 0.0, 0.0),
            &eye_v,
            &normal_v,
            false,
        );

        assert_abs_diff_eq!(c1, Color::WHITE);
        assert_abs_diff_eq!(c2, Color::BLACK);
    }
}
