use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_4};

use raytracer::{
    core::{matrix::Matrix, tuple::Tuple},
    graphics::color::Color,
    raytracer::{
        camera::Camera, material::Material, object::Object, point_light::PointLight, world::World,
    },
};

fn render_scene(path: &str) {
    let floor_material = Material {
        color: Color::new(1.0, 0.9, 0.9),
        specular: 0.0,
        ..Material::default()
    };

    let mut floor = Object::new_sphere();
    floor.transform = Matrix::scaling(10.0, 0.01, 10.0);
    floor.material = floor_material;

    let mut left_wall = Object::new_sphere();
    left_wall.transform = Matrix::translation(0.0, 0.0, 5.0)
        * Matrix::rotation_y(-FRAC_PI_4)
        * Matrix::rotation_x(FRAC_PI_2)
        * Matrix::scaling(10.0, 0.01, 10.0);
    left_wall.material = floor_material;

    let mut right_wall = Object::new_sphere();
    right_wall.transform = Matrix::translation(0.0, 0.0, 5.0)
        * Matrix::rotation_y(FRAC_PI_4)
        * Matrix::rotation_x(FRAC_PI_2)
        * Matrix::scaling(10.0, 0.01, 10.0);
    right_wall.material = floor_material;

    let mut middle = Object::new_sphere();
    middle.transform = Matrix::translation(-0.5, 1.0, 0.5);
    middle.material = Material {
        color: Color::new(0.1, 1.0, 0.5),
        diffuse: 0.7,
        specular: 0.3,
        ..Material::default()
    };

    let mut right = Object::new_sphere();
    right.transform = Matrix::translation(1.5, 0.5, -0.5) * Matrix::scaling(0.5, 0.5, 0.5);
    right.material = Material {
        color: Color::new(0.5, 1.0, 0.1),
        diffuse: 0.7,
        specular: 0.3,
        ..Material::default()
    };

    let mut left = Object::new_sphere();
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
        objects: vec![floor, left_wall, right_wall, middle, right, left],
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
    render_scene("images/scene.ppm");
}
