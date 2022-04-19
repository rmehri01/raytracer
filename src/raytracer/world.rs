use crate::{
    core::{matrix::Matrix, point::Point},
    graphics::color::Color,
};

use super::{
    intersection::{Computations, Intersections},
    material::Material,
    point_light::PointLight,
    ray::Ray,
    shape::Shape,
};

#[derive(Debug, PartialEq)]
pub struct World {
    pub shapes: Vec<Shape>,
    // TODO: does this need to be optional
    pub light: Option<PointLight>,
}

impl World {
    pub fn empty() -> Self {
        Self {
            shapes: Vec::new(),
            light: None,
        }
    }

    pub fn color_at(&self, ray: &Ray, remaining_recursions: u8) -> Color {
        let intersections = self.intersect(ray);
        let hit = intersections.hit();

        match hit {
            Some(hit) => {
                let comps = hit.prepare_computations(ray, &intersections);
                self.shade_hit(&comps, remaining_recursions)
            }
            None => Color::BLACK,
        }
    }

    fn intersect(&self, ray: &Ray) -> Intersections {
        let intersects = self
            .shapes
            .iter()
            .flat_map(|shape| shape.intersect(ray, im::Vector::new()).0)
            .collect();

        Intersections(intersects)
    }

    fn shade_hit(&self, comps: &Computations, remaining_recursions: u8) -> Color {
        let is_shadowed = self.is_shadowed(&comps.over_point);
        let surface = comps.shape.material.lighting(
            &comps.shape.world_to_object(&comps.over_point, &comps.trail),
            &comps.over_point,
            &self.light.expect("world should have light"),
            &comps.eye_v,
            &comps.normal_v,
            is_shadowed,
        );

        let reflected = self.reflected_color(comps, remaining_recursions);
        let refracted = self.refracted_color(comps, remaining_recursions);

        let material = comps.shape.material;
        if material.reflective > 0.0 && material.transparency > 0.0 {
            let reflectance = comps.schlick();
            surface + reflected * reflectance + refracted * (1.0 - reflectance)
        } else {
            surface + reflected + refracted
        }
    }

    fn is_shadowed(&self, point: &Point) -> bool {
        let v = self.light.expect("world should have light").position - *point;
        let distance = v.magnitude();
        let direction = v.normalize();

        let r = Ray::new(*point, direction);
        let intersections = self.intersect(&r);

        intersections.hit().map_or(false, |hit| hit.t < distance)
    }

    fn reflected_color(&self, comps: &Computations, remaining_recursions: u8) -> Color {
        if comps.shape.material.reflective == 0.0 || remaining_recursions == 0 {
            Color::BLACK
        } else {
            let reflect_ray = Ray::new(comps.over_point, comps.reflect_v);
            let color = self.color_at(&reflect_ray, remaining_recursions - 1);

            color * comps.shape.material.reflective
        }
    }

    fn refracted_color(&self, comps: &Computations, remaining_recursions: u8) -> Color {
        let n_ratio = comps.n1 / comps.n2;
        let cos_i = comps.eye_v.dot(&comps.normal_v);
        let sin2_t = n_ratio.powi(2) * (1.0 - cos_i.powi(2));

        if comps.shape.material.transparency == 0.0 || remaining_recursions == 0 || sin2_t > 1.0 {
            Color::BLACK
        } else {
            let cos_t = (1.0 - sin2_t).sqrt();
            let direction = comps.normal_v * (n_ratio * cos_i - cos_t) - comps.eye_v * n_ratio;
            let refract_ray = Ray::new(comps.under_point, direction);

            self.color_at(&refract_ray, remaining_recursions - 1)
                * comps.shape.material.transparency
        }
    }
}

impl Default for World {
    fn default() -> Self {
        let light = Some(PointLight::new(
            Point::new(-10.0, 10.0, -10.0),
            Color::WHITE,
        ));

        let s1 = Shape::new_sphere().with_material(Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ..Material::default()
        });

        let s2 = Shape::new_sphere().with_transform(Matrix::scaling(0.5, 0.5, 0.5));

        Self {
            shapes: vec![s1, s2],
            light,
        }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use crate::{
        core::vector::Vector, graphics::pattern::Pattern, raytracer::intersection::Intersection,
    };

    use super::*;

    #[test]
    fn empty_world() {
        let world = World::empty();

        assert_eq!(world.shapes.len(), 0);
        assert_eq!(world.light, None);
    }

    #[test]
    fn default_world() {
        let world = World::default();

        assert_eq!(world.shapes.len(), 2);
        assert_eq!(
            world.light.expect("light exists"),
            PointLight::new(Point::new(-10.0, 10.0, -10.0), Color::WHITE)
        );
    }

    #[test]
    fn intersect_world_with_ray() {
        let world = World::default();

        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));

        let xs = world
            .intersect(&r)
            .0
            .iter()
            .map(|x| x.t)
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
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let shape = &world.shapes[0];
        let i = Intersection::new(4.0, shape, im::Vector::new());
        let comps = i.prepare_computations(&r, &Intersections::new([i.clone()]));

        let color = world.shade_hit(&comps, 5);

        assert_abs_diff_eq!(color, Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn shade_intersection_from_inside() {
        let world = World {
            light: Some(PointLight::new(Point::new(0.0, 0.25, 0.0), Color::WHITE)),
            ..World::default()
        };

        let r = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));
        let shape = &world.shapes[1];
        let i = Intersection::new(0.5, shape, im::Vector::new());
        let comps = i.prepare_computations(&r, &Intersections::new([i.clone()]));

        let color = world.shade_hit(&comps, 5);

        assert_abs_diff_eq!(color, Color::new(0.90498, 0.90498, 0.90498));
    }

    #[test]
    fn color_ray_misses() {
        let world = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 1.0, 0.0));

        let color = world.color_at(&r, 5);

        assert_abs_diff_eq!(color, Color::BLACK);
    }

    #[test]
    fn color_ray_hits() {
        let world = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));

        let color = world.color_at(&r, 5);

        assert_abs_diff_eq!(color, Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn color_intersection_behind_ray() {
        let mut world = World::default();
        world.shapes[0].material.ambient = 1.0;
        world.shapes[1].material.ambient = 1.0;

        let r = Ray::new(Point::new(0.0, 0.0, 0.75), Vector::new(0.0, 0.0, -1.0));

        let color = world.color_at(&r, 5);

        assert_abs_diff_eq!(color, world.shapes[1].material.color);
    }

    #[test]
    fn no_shadow_when_no_shape_between_light_and_point() {
        let w = World::default();
        let p = Point::new(0.0, 10.0, 0.0);

        assert!(!w.is_shadowed(&p));
    }

    #[test]
    fn shadow_when_shape_between_light_and_point() {
        let w = World::default();
        let p = Point::new(10.0, -10.0, 10.0);

        assert!(w.is_shadowed(&p));
    }

    #[test]
    fn no_shadow_when_shape_behind_light() {
        let w = World::default();
        let p = Point::new(-20.0, 20.0, -20.0);

        assert!(!w.is_shadowed(&p));
    }

    #[test]
    fn no_shadow_when_shape_behind_point() {
        let w = World::default();
        let p = Point::new(-2.0, 2.0, -2.0);

        assert!(!w.is_shadowed(&p));
    }

    #[test]
    fn shade_hit_given_intersection_in_shadow() {
        let s1 = Shape::new_sphere();
        let s2 = Shape::new_sphere().with_transform(Matrix::translation(0.0, 0.0, 10.0));

        let w = World {
            light: Some(PointLight::new(Point::new(0.0, 0.0, -10.0), Color::WHITE)),
            shapes: vec![s1, s2.clone()],
        };

        let r = Ray::new(Point::new(0.0, 0.0, 5.0), Vector::new(0.0, 0.0, 1.0));
        let i = Intersection::new(4.0, &s2, im::Vector::new());
        let comps = i.prepare_computations(&r, &Intersections::new([i.clone()]));

        assert_abs_diff_eq!(w.shade_hit(&comps, 5), Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn reflected_color_nonreflective_material() {
        let mut w = World::default();
        w.shapes[1].material.ambient = 1.0;

        let r = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));
        let i = Intersection::new(1.0, &w.shapes[1], im::Vector::new());

        let comps = i.prepare_computations(&r, &Intersections::new([i.clone()]));

        assert_abs_diff_eq!(w.reflected_color(&comps, 5), Color::BLACK);
    }

    #[test]
    fn reflected_color_reflective_material() {
        let shape = Shape::new_plane()
            .with_transform(Matrix::translation(0.0, -1.0, 0.0))
            .with_material(Material {
                reflective: 0.5,
                ..Material::default()
            });

        let mut w = World::default();
        w.shapes.push(shape.clone());

        let r = Ray::new(
            Point::new(0.0, 0.0, -3.0),
            Vector::new(0.0, -(2.0_f64.sqrt()) / 2.0, 2.0_f64.sqrt() / 2.0),
        );
        let i = Intersection::new(2.0_f64.sqrt(), &shape, im::Vector::new());

        let comps = i.prepare_computations(&r, &Intersections::new([i.clone()]));

        assert_abs_diff_eq!(
            w.reflected_color(&comps, 5),
            Color::new(0.19032, 0.2379, 0.14274)
        );
    }

    #[test]
    fn shade_hit_with_reflective_material() {
        let shape = Shape::new_plane()
            .with_transform(Matrix::translation(0.0, -1.0, 0.0))
            .with_material(Material {
                reflective: 0.5,
                ..Material::default()
            });

        let mut w = World::default();
        w.shapes.push(shape.clone());

        let r = Ray::new(
            Point::new(0.0, 0.0, -3.0),
            Vector::new(0.0, -(2.0_f64.sqrt()) / 2.0, 2.0_f64.sqrt() / 2.0),
        );
        let i = Intersection::new(2.0_f64.sqrt(), &shape, im::Vector::new());

        let comps = i.prepare_computations(&r, &Intersections::new([i.clone()]));

        assert_abs_diff_eq!(
            w.shade_hit(&comps, 5),
            Color::new(0.87677, 0.92436, 0.82918)
        );
    }

    #[test]
    fn mutually_reflective_surfaces() {
        let lower = Shape::new_plane()
            .with_transform(Matrix::translation(0.0, -1.0, 0.0))
            .with_material(Material {
                reflective: 1.0,
                ..Material::default()
            });

        let upper = Shape::new_plane()
            .with_transform(Matrix::translation(0.0, 1.0, 0.0))
            .with_material(Material {
                reflective: 1.0,
                ..Material::default()
            });

        let mut w = World::default();
        w.shapes.push(lower);
        w.shapes.push(upper);

        let r = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0));

        w.color_at(&r, 5);
    }

    #[test]
    fn reflected_color_at_maximum_recursion_depth() {
        let shape = Shape::new_sphere()
            .with_transform(Matrix::translation(0.0, -1.0, 0.0))
            .with_material(Material {
                reflective: 0.5,
                ..Material::default()
            });

        let mut w = World::default();
        w.shapes.push(shape.clone());

        let r = Ray::new(
            Point::new(0.0, 0.0, -3.0),
            Vector::new(0.0, -(2.0_f64.sqrt()) / 2.0, 2.0_f64.sqrt() / 2.0),
        );
        let i = Intersection::new(2.0_f64.sqrt(), &shape, im::Vector::new());

        let comps = i.prepare_computations(&r, &Intersections::new([i.clone()]));

        assert_abs_diff_eq!(w.reflected_color(&comps, 0), Color::BLACK);
    }

    #[test]
    fn refracted_color_of_opaque_surface() {
        let w = World::default();
        let shape = &w.shapes[0];

        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let xs = Intersections::new([
            Intersection::new(4.0, shape, im::Vector::new()),
            Intersection::new(6.0, shape, im::Vector::new()),
        ]);

        let comps = xs.0.iter().next().unwrap().prepare_computations(&r, &xs);

        assert_eq!(w.refracted_color(&comps, 5), Color::BLACK);
    }

    #[test]
    fn refracted_color_at_max_recursion_depth() {
        let mut w = World::default();
        w.shapes[0].material.transparency = 1.0;
        w.shapes[0].material.refractive_index = 1.5;

        let shape = &w.shapes[0];

        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let xs = Intersections::new([
            Intersection::new(4.0, shape, im::Vector::new()),
            Intersection::new(6.0, shape, im::Vector::new()),
        ]);

        let comps = xs.0.iter().next().unwrap().prepare_computations(&r, &xs);

        assert_eq!(w.refracted_color(&comps, 0), Color::BLACK);
    }

    #[test]
    fn refracted_color_under_total_internal_refraction() {
        let mut w = World::default();
        w.shapes[0].material.transparency = 1.0;
        w.shapes[0].material.refractive_index = 1.5;

        let shape = &w.shapes[0];

        let r = Ray::new(
            Point::new(0.0, 0.0, 2.0_f64.sqrt() / 2.0),
            Vector::new(0.0, 1.0, 0.0),
        );
        let xs = Intersections::new([
            Intersection::new(-(2.0_f64.sqrt()) / 2.0, shape, im::Vector::new()),
            Intersection::new(2.0_f64.sqrt() / 2.0, shape, im::Vector::new()),
        ]);

        let comps = xs.0.iter().nth(1).unwrap().prepare_computations(&r, &xs);

        assert_eq!(w.refracted_color(&comps, 5), Color::BLACK);
    }

    #[test]
    fn refracted_color_with_refracted_ray() {
        let mut w = World::default();
        w.shapes[0].material.ambient = 1.0;
        w.shapes[0].material.pattern = Some(Pattern::new_test());
        w.shapes[1].material.transparency = 1.0;
        w.shapes[1].material.refractive_index = 1.5;

        let r = Ray::new(Point::new(0.0, 0.0, 0.1), Vector::new(0.0, 1.0, 0.0));
        let xs = Intersections::new([
            Intersection::new(-0.9899, &w.shapes[0], im::Vector::new()),
            Intersection::new(-0.4899, &w.shapes[1], im::Vector::new()),
            Intersection::new(0.4899, &w.shapes[1], im::Vector::new()),
            Intersection::new(0.9899, &w.shapes[0], im::Vector::new()),
        ]);

        let comps = xs.0.iter().nth(2).unwrap().prepare_computations(&r, &xs);

        assert_eq!(
            w.refracted_color(&comps, 5),
            Color::new(0.0, 0.99888, 0.04725)
        );
    }

    #[test]
    fn shade_hit_with_transparent_material() {
        let floor = Shape::new_plane()
            .with_transform(Matrix::translation(0.0, -1.0, 0.0))
            .with_material(Material {
                transparency: 0.5,
                refractive_index: 1.5,
                ..Material::default()
            });

        let ball = Shape::new_sphere()
            .with_transform(Matrix::translation(0.0, -3.5, -0.5))
            .with_material(Material {
                color: Color::new(1.0, 0.0, 0.0),
                ambient: 0.5,
                ..Material::default()
            });

        let mut w = World::default();
        w.shapes.push(floor.clone());
        w.shapes.push(ball);

        let r = Ray::new(
            Point::new(0.0, 0.0, -3.0),
            Vector::new(0.0, -(2.0_f64.sqrt()) / 2.0, 2.0_f64.sqrt() / 2.0),
        );
        let xs = Intersections::new([Intersection::new(2.0_f64.sqrt(), &floor, im::Vector::new())]);

        let comps = xs.0.iter().next().unwrap().prepare_computations(&r, &xs);

        assert_eq!(
            w.shade_hit(&comps, 5),
            Color::new(0.93642, 0.68642, 0.68642)
        );
    }

    #[test]
    fn shade_hit_with_reflective_transparent_material() {
        let floor = Shape::new_plane()
            .with_transform(Matrix::translation(0.0, -1.0, 0.0))
            .with_material(Material {
                reflective: 0.5,
                transparency: 0.5,
                refractive_index: 1.5,
                ..Material::default()
            });

        let ball = Shape::new_sphere()
            .with_transform(Matrix::translation(0.0, -3.5, -0.5))
            .with_material(Material {
                color: Color::new(1.0, 0.0, 0.0),
                ambient: 0.5,
                ..Material::default()
            });

        let mut w = World::default();
        w.shapes.push(floor.clone());
        w.shapes.push(ball);

        let r = Ray::new(
            Point::new(0.0, 0.0, -3.0),
            Vector::new(0.0, -(2.0_f64.sqrt()) / 2.0, 2.0_f64.sqrt() / 2.0),
        );
        let xs = Intersections::new([Intersection::new(2.0_f64.sqrt(), &floor, im::Vector::new())]);

        let comps = xs.0.iter().next().unwrap().prepare_computations(&r, &xs);

        assert_eq!(
            w.shade_hit(&comps, 5),
            Color::new(0.93391, 0.69643, 0.69243)
        );
    }
}
