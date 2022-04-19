use std::f64::consts::FRAC_PI_3;

use raytracer::{
    core::{matrix::Matrix, tuple::Tuple},
    graphics::{color::Color, pattern::Pattern},
    io::obj,
    raytracer::{
        camera::Camera, material::Material, point_light::PointLight, shape::Shape, world::World,
    },
};

fn main() {
    render_teapot("images/teapot.ppm");
}

fn render_teapot(path: &str) {
    let floor_material = Material {
        color: Color::new(1.0, 0.9, 0.9),
        specular: 0.0,
        pattern: Some(Pattern::new_checker(Color::WHITE, Color::BLACK)),
        ..Material::default()
    };

    let floor = Shape::new_plane().with_material(floor_material);
    let teapot = obj::parse_file("exercises/resources/teapot.obj").unwrap();

    let world = World {
        light: Some(PointLight::new(
            Tuple::point(-10.0, 10.0, -10.0),
            Color::WHITE,
        )),
        shapes: vec![floor, teapot],
    };

    let camera = Camera::new(500, 500, FRAC_PI_3).with_transform(Matrix::view_transform(
        Tuple::point(0.0, 1.5, -5.0),
        Tuple::point(0.0, 1.0, 0.0),
        Tuple::vector(0.0, 1.0, 0.0),
    ));

    camera.render(world).write_ppm(path);
}
