use std::collections::BTreeSet;

use crate::{
    core::{matrix::Matrix, tuple::Tuple},
    graphics::color::Color,
};

use super::{
    intersection::Intersections, material::Material, point_light::PointLight, ray::Ray,
    sphere::Sphere,
};

#[derive(Debug, PartialEq)]
pub struct World {
    objects: Vec<Sphere>,
    light: Option<PointLight>,
}

impl World {
    pub fn empty() -> Self {
        Self {
            objects: Vec::new(),
            light: None,
        }
    }

    pub fn intersect(&self, ray: &Ray) -> Intersections {
        let intersects = self
            .objects
            .iter()
            .flat_map(|object| object.intersect(ray).0)
            .collect::<BTreeSet<_>>();

        Intersections(intersects)
    }
}

impl Default for World {
    fn default() -> Self {
        let light = Some(PointLight::new(
            Tuple::point(-10.0, 10.0, -10.0),
            Color::new(1.0, 1.0, 1.0),
        ));

        let s1 = Sphere {
            material: Material {
                color: Color::new(0.8, 1.0, 0.6),
                diffuse: 0.7,
                specular: 0.2,
                ..Material::default()
            },
            ..Sphere::default()
        };
        let s2 = Sphere {
            transform: Matrix::scaling(0.5, 0.5, 0.5),
            ..Sphere::default()
        };

        Self {
            objects: vec![s1, s2],
            light,
        }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn empty_world() {
        let world = World::empty();

        assert_eq!(world.objects.len(), 0);
        assert_eq!(world.light, None);
    }

    #[test]
    fn default_world() {
        let world = World::default();

        assert_eq!(world.objects.len(), 2);
        assert_abs_diff_eq!(
            world.light.expect("light exists"),
            PointLight::new(Tuple::point(-10.0, 10.0, -10.0), Color::new(1.0, 1.0, 1.0))
        );
    }

    #[test]
    fn intersect_world_with_ray() {
        let world = World::default();

        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = world
            .intersect(&r)
            .0
            .iter()
            .map(|x| x.t.0)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 4);
        assert_abs_diff_eq!(xs[0], 4.0);
        assert_abs_diff_eq!(xs[1], 4.5);
        assert_abs_diff_eq!(xs[2], 5.5);
        assert_abs_diff_eq!(xs[3], 6.0);
    }
}
