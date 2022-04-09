use std::collections::BTreeSet;

use crate::{
    core::{matrix::Matrix, tuple::Tuple},
    graphics::color::Color,
};

use super::{
    intersection::{Computations, Intersection, Intersections},
    material::Material,
    point_light::PointLight,
    ray::Ray,
    sphere::Sphere,
};

#[derive(Debug, PartialEq)]
pub struct World {
    pub objects: Vec<Sphere>,
    // TODO: does this need to be optional
    pub light: Option<PointLight>,
}

impl World {
    pub fn empty() -> Self {
        Self {
            objects: Vec::new(),
            light: None,
        }
    }

    pub fn color_at(&self, ray: &Ray) -> Color {
        let intersections = self.intersect(ray);
        let hit = intersections.hit();

        match hit {
            Some(hit) => {
                let comps = Self::prepare_computations(&hit, ray);
                self.shade_hit(comps)
            }
            None => Color::new(0.0, 0.0, 0.0),
        }
    }

    fn intersect(&self, ray: &Ray) -> Intersections {
        let intersects = self
            .objects
            .iter()
            .flat_map(|object| object.intersect(ray).0)
            .collect::<BTreeSet<_>>();

        Intersections(intersects)
    }

    fn shade_hit(&self, comps: Computations) -> Color {
        comps.object.material.lighting(
            &self.light.expect("world should have light"),
            &comps.point,
            &comps.eyev,
            &comps.normalv,
        )
    }

    fn prepare_computations(intersection: &Intersection, ray: &Ray) -> Computations {
        let point = ray.position(intersection.t.0);
        let eyev = -ray.direction;
        let normalv = intersection.object.normal_at(&point);
        let inside = normalv.dot(&eyev) < 0.0;

        Computations {
            t: intersection.t.0,
            object: intersection.object,
            point,
            eyev,
            normalv: if inside { -normalv } else { normalv },
            inside,
        }
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

    #[test]
    fn shade_intersection() {
        let world = World::default();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = &world.objects[0];
        let i = Intersection::new(4.0, *shape);
        let comps = World::prepare_computations(&i, &r);

        let color = world.shade_hit(comps);

        assert_abs_diff_eq!(color, Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn shade_intersection_from_inside() {
        let world = World {
            light: Some(PointLight::new(
                Tuple::point(0.0, 0.25, 0.0),
                Color::new(1.0, 1.0, 1.0),
            )),
            ..World::default()
        };

        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = &world.objects[1];
        let i = Intersection::new(0.5, *shape);
        let comps = World::prepare_computations(&i, &r);

        let color = world.shade_hit(comps);

        assert_abs_diff_eq!(color, Color::new(0.90498, 0.90498, 0.90498));
    }

    #[test]
    fn precomputing_state_of_intersection() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::default();
        let i = Intersection::new(4.0, shape);

        let comps = World::prepare_computations(&i, &r);

        assert_abs_diff_eq!(comps.t, i.t);
        assert_abs_diff_eq!(comps.object, i.object);
        assert_abs_diff_eq!(comps.point, Tuple::point(0.0, 0.0, -1.0));
        assert_abs_diff_eq!(comps.eyev, Tuple::vector(0.0, 0.0, -1.0));
        assert_abs_diff_eq!(comps.normalv, Tuple::vector(0.0, 0.0, -1.0));
    }

    #[test]
    fn hit_when_intersection_outside() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::default();
        let i = Intersection::new(4.0, shape);

        let comps = World::prepare_computations(&i, &r);

        assert!(!comps.inside);
    }

    #[test]
    fn hit_when_intersection_inside() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Sphere::default();
        let i = Intersection::new(1.0, shape);

        let comps = World::prepare_computations(&i, &r);

        assert_abs_diff_eq!(comps.point, Tuple::point(0.0, 0.0, 1.0));
        assert_abs_diff_eq!(comps.eyev, Tuple::vector(0.0, 0.0, -1.0));
        assert!(comps.inside);
        assert_abs_diff_eq!(comps.normalv, Tuple::vector(0.0, 0.0, -1.0));
    }

    #[test]
    fn color_ray_misses() {
        let world = World::default();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 1.0, 0.0));

        let color = world.color_at(&r);

        assert_abs_diff_eq!(color, Color::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn color_ray_hits() {
        let world = World::default();
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));

        let color = world.color_at(&r);

        assert_abs_diff_eq!(color, Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn color_intersection_behind_ray() {
        let mut world = World::default();
        world.objects[0].material.ambient = 1.0;
        world.objects[1].material.ambient = 1.0;

        let r = Ray::new(Tuple::point(0.0, 0.0, 0.75), Tuple::vector(0.0, 0.0, -1.0));

        let color = world.color_at(&r);

        assert_abs_diff_eq!(color, world.objects[1].material.color);
    }
}
