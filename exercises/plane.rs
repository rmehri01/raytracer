use std::f64::consts::FRAC_PI_3;

use raytracer::{
    core::{matrix::Matrix, tuple::Tuple},
    graphics::color::Color,
    raytracer::{
        camera::Camera, material::Material, point_light::PointLight, shape::Shape, world::World,
    },
};

fn render_plane_scene(path: &str) {
    let floor_material = Material {
        color: Color::new(1.0, 0.9, 0.9),
        specular: 0.0,
        ..Material::default()
    };

    let mut floor = Shape::new_plane();
    floor.material = floor_material;

    let mut middle = Shape::new_sphere();
    middle.transform = Matrix::translation(-0.5, 1.0, 0.5);
    middle.material = Material {
        color: Color::new(0.1, 1.0, 0.5),
        diffuse: 0.7,
        specular: 0.3,
        ..Material::default()
    };

    let mut right = Shape::new_sphere();
    right.transform = Matrix::translation(1.5, 0.5, -0.5) * Matrix::scaling(0.5, 0.5, 0.5);
    right.material = Material {
        color: Color::new(0.5, 1.0, 0.1),
        diffuse: 0.7,
        specular: 0.3,
        ..Material::default()
    };

    let mut left = Shape::new_sphere();
    left.transform = Matrix::translation(-1.5, 0.33, -0.75) * Matrix::scaling(0.33, 0.33, 0.33);
    left.material = Material {
        color: Color::new(1.0, 0.8, 0.1),
        diffuse: 0.7,
        specular: 0.3,
        ..Material::default()
    };

    let world = World {
        light: Some(PointLight::new(
            Tuple::point(-10.0, 10.0, -10.0),
            Color::WHITE,
        )),
        shapes: vec![floor, middle, right, left],
    };

    let mut camera = Camera::new(2000, 1000, FRAC_PI_3);
    camera.transform = Matrix::view_transform(
        Tuple::point(0.0, 1.5, -5.0),
        Tuple::point(0.0, 1.0, 0.0),
        Tuple::vector(0.0, 1.0, 0.0),
    );

    camera.render(world).write_ppm(path);
}

fn main() {
    render_plane_scene("images/plane.ppm");
}
