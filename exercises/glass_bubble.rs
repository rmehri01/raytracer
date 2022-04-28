use std::f64::consts::FRAC_PI_3;

use raytracer::{
    core::{matrix::Matrix, point::Point, vector::Vector},
    graphics::{color::Color, pattern::Pattern},
    raytracer::{
        camera::Camera,
        material::Material,
        point_light::PointLight,
        shapes::{Primitive, SetProperties},
        world::World,
    },
};

fn main() {
    render_glass_bubble("images/glass_bubble.ppm");
}

fn render_glass_bubble(path: &str) {
    let floor_pattern = Pattern::new_checker(
        Pattern::new_solid(Color::WHITE),
        Pattern::new_solid(Color::BLACK),
    );
    let floor = Primitive::new_plane()
        .with_material(Material {
            specular: 0.0,
            pattern: Some(floor_pattern),
            ..Material::default()
        })
        .with_transform(Matrix::translation(0.0, -10.0, 0.0))
        .as_shape();

    let glass = Primitive::new_sphere()
        .with_material(Material {
            diffuse: 0.1,
            shininess: 300.0,
            reflective: 1.0,
            transparency: 1.0,
            refractive_index: 1.52,
            ..Material::default()
        })
        .as_shape();

    let air = Primitive::new_sphere()
        .with_material(Material {
            diffuse: 0.1,
            shininess: 300.0,
            reflective: 1.0,
            transparency: 1.0,
            refractive_index: 1.0,
            ..Material::default()
        })
        .with_transform(Matrix::scaling(0.5, 0.5, 0.5))
        .as_shape();

    let world = World::new(
        vec![floor, glass, air],
        vec![PointLight::new(
            Point::new(20.0, 10.0, 0.0),
            Color::new(0.6, 0.6, 0.6),
        )],
    );

    let camera = Camera::new(2048, 1080, FRAC_PI_3).with_transform(Matrix::view_transform(
        Point::new(0.0, 3.5, 0.0),
        Point::new(0.0, 0.0, 0.0),
        Vector::new(0.0, 0.0, 1.0),
    ));

    camera.render(&world).write_ppm(path);
}
