use enum_dispatch::enum_dispatch;

use crate::{
    core::matrix::{Matrix, Transformation},
    raytracer::{intersection::Intersections, material::Material, ray::Ray},
};

use super::{Compound, Single};

#[enum_dispatch]
#[derive(Debug, PartialEq, Clone)]
pub enum Shape {
    Single,
    Compound,
}

impl Shape {
    pub fn includes(&self, other: &Single) -> bool {
        match self {
            Shape::Single(s) => s == other,
            Shape::Compound(c) => c.includes(other),
        }
    }
}

#[enum_dispatch(Shape)]
pub trait Intersect: HasProperties {
    fn intersect(&self, ray: &Ray, trail: &im_rc::Vector<Transformation>) -> Intersections {
        let local_ray = ray.transform(&self.properties().transform_inversed);

        self.local_intersect(&local_ray, trail)
    }

    fn local_intersect(&self, ray: &Ray, trail: &im_rc::Vector<Transformation>) -> Intersections;
}

#[derive(Debug, PartialEq, Clone)]
pub struct Properties {
    pub transform: Transformation,
    pub transform_inversed: Transformation,
    pub material: Material,
}

impl Default for Properties {
    fn default() -> Self {
        Self {
            transform: Matrix::identity(),
            transform_inversed: Matrix::identity(),
            material: Material::default(),
        }
    }
}

#[enum_dispatch(Shape)]
pub trait HasProperties: Sized {
    fn properties(&self) -> &Properties;
    fn properties_mut(&mut self) -> &mut Properties;
}

pub trait SetProperties: HasProperties {
    fn with_transform(mut self, transform: Transformation) -> Self {
        self.properties_mut().transform = transform;
        self.properties_mut().transform_inversed = transform.inverse();
        self
    }

    fn with_material(mut self, material: Material) -> Self {
        self.properties_mut().material = material;
        self
    }
}

impl<T: HasProperties> SetProperties for T {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_transformation() {
        let shape = Single::new_sphere();

        assert_eq!(shape.properties.transform, Matrix::identity());
    }

    #[test]
    fn change_transformation() {
        let shape = Single::new_sphere().with_transform(Matrix::translation(2.0, 3.0, 4.0));

        assert_eq!(
            shape.properties.transform,
            Matrix::translation(2.0, 3.0, 4.0)
        );
    }

    #[test]
    fn default_material() {
        let shape = Single::new_sphere();

        assert_eq!(shape.properties.material, Material::default());
    }

    #[test]
    fn assign_material() {
        let material = Material {
            ambient: 1.0,
            ..Material::default()
        };

        let shape = Single::new_sphere().with_material(material);

        assert_eq!(shape.properties.material, material);
    }
}
