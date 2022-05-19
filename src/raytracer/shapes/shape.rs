use enum_dispatch::enum_dispatch;

use crate::{
    core::matrix::{Matrix, Transformation},
    raytracer::{bounds::Bounds, intersection::Intersections, material::Material, ray::Ray},
};

use super::{Compound, Primitive};

#[enum_dispatch]
#[derive(Debug, PartialEq, Clone)]
pub enum Shape {
    Primitive,
    Compound,
}

impl Shape {
    pub fn includes(&self, other: &Primitive) -> bool {
        match self {
            Shape::Primitive(p) => p == other,
            Shape::Compound(c) => c.includes(other),
        }
    }
}

#[enum_dispatch(Shape)]
pub trait Intersect: HasProperties {
    fn intersect(&self, ray: &Ray, trail: &im_rc::Vector<Transformation>) -> Intersections {
        let local_ray = ray.transform(&self.properties().inverse_transform);

        self.local_intersect(&local_ray, trail)
    }

    fn local_intersect(&self, ray: &Ray, trail: &im_rc::Vector<Transformation>) -> Intersections;
}

#[derive(Debug, PartialEq, Clone)]
pub struct Properties {
    pub transform: Transformation,
    pub(crate) inverse_transform: Transformation,
    pub material: Material,
    pub(crate) bounds: Bounds,
}

impl Default for Properties {
    fn default() -> Self {
        Self {
            transform: Matrix::identity(),
            inverse_transform: Matrix::identity(),
            material: Material::default(),
            bounds: Bounds::default(),
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
        self.properties_mut().inverse_transform = transform.inverse();
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
        let shape = Primitive::new_sphere();

        assert_eq!(shape.properties().transform, Matrix::identity());
    }

    #[test]
    fn change_transformation() {
        let shape = Primitive::new_sphere().with_transform(Matrix::translation(2.0, 3.0, 4.0));

        assert_eq!(
            shape.properties().transform,
            Matrix::translation(2.0, 3.0, 4.0)
        );
    }

    #[test]
    fn default_material() {
        let shape = Primitive::new_sphere();

        assert_eq!(shape.properties().material, Material::default());
    }

    #[test]
    fn assign_material() {
        let material = Material {
            ambient: 1.0,
            ..Material::default()
        };

        let shape = Primitive::new_sphere().with_material(material.clone());

        assert_eq!(shape.properties().material, material);
    }
}
