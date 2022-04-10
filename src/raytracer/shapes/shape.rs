use crate::{core::tuple::Tuple, raytracer::ray::Ray};

use super::{plane::Plane, sphere::Sphere};

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Shape {
    Sphere(Sphere),
    Plane(Plane),
}

// TODO: maybe trait, some way to get around this?
impl Shape {
    pub fn intersect(&self, ray: &Ray) -> Vec<f64> {
        match self {
            Shape::Sphere(s) => s.intersect(ray),
            Shape::Plane(p) => p.intersect(ray),
        }
    }

    pub fn normal_at(&self, point: &Tuple) -> Tuple {
        match self {
            Shape::Sphere(s) => s.normal_at(point),
            Shape::Plane(p) => p.normal_at(point),
        }
    }
}
