use crate::{core::tuple::Tuple, raytracer::ray::Ray};

use super::sphere::Sphere;

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Shape {
    Sphere(Sphere),
}

impl Shape {
    pub fn intersect(&self, ray: &Ray) -> Vec<f64> {
        match self {
            Shape::Sphere(s) => s.intersect(ray),
        }
    }

    pub fn normal_at(&self, point: &Tuple) -> Tuple {
        match self {
            Shape::Sphere(s) => s.normal_at(point),
        }
    }
}
