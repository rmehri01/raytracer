#![allow(dead_code)]

use crate::{
    core::tuple::Tuple,
    graphics::{canvas::Canvas, color::Color},
    raytracer::{object::Sphere, ray::Ray},
};

pub fn plot_sphere(path: &str) {
    let side_len = 500;
    let mut canvas = Canvas::new(side_len, side_len);

    let ray_origin = Tuple::point(0.0, 0.0, -5.0);
    let wall_z = 10.0;
    let wall_size = 7.0;

    let pixel_size = wall_size / side_len as f64;
    let half = wall_size / 2.0;

    let color = Color::new(1.0, 0.0, 0.0);
    let sphere = Sphere::new();

    for y in 0..side_len {
        let world_y = half - pixel_size * y as f64;

        for x in 0..side_len {
            let world_x = -half + pixel_size * x as f64;
            let position = Tuple::point(world_x, world_y, wall_z);

            let r = Ray::new(ray_origin, (position - ray_origin).normalize());
            let xs = sphere.intersect(&r);

            if xs.hit().is_some() {
                canvas.write_pixel(x, y, color);
            }
        }
    }

    canvas.write_ppm(path);
}