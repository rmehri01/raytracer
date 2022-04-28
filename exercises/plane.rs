use std::f64::consts::FRAC_PI_3;

use raytracer::{
    core::{matrix::Matrix, point::Point, vector::Vector},
    graphics::color::Color,
    raytracer::{
        camera::Camera,
        material::Material,
        point_light::PointLight,
        shapes::{Primitive, SetProperties},
        world::World,
    },
};

fn main() {
    render_plane_scene("images/ppm/plane.ppm");
}

fn render_plane_scene(path: &str) {
    let floor_material = Material {
        color: Color::new(1.0, 0.9, 0.9),
        specular: 0.0,
        ..Material::default()
    };

    let floor = Primitive::new_plane()
        .with_material(floor_material)
        .as_shape();

    let middle = Primitive::new_sphere()
        .with_transform(Matrix::translation(-0.5, 1.0, 0.5))
        .with_material(Material {
            color: Color::new(0.1, 1.0, 0.5),
            diffuse: 0.7,
            specular: 0.3,
            ..Material::default()
        })
        .as_shape();

    let right = Primitive::new_sphere()
        .with_transform(Matrix::translation(1.5, 0.5, -0.5) * Matrix::scaling(0.5, 0.5, 0.5))
        .with_material(Material {
            color: Color::new(0.5, 1.0, 0.1),
            diffuse: 0.7,
            specular: 0.3,
            ..Material::default()
        })
        .as_shape();

    let left = Primitive::new_sphere()
        .with_transform(Matrix::translation(-1.5, 0.33, -0.75) * Matrix::scaling(0.33, 0.33, 0.33))
        .with_material(Material {
            color: Color::new(1.0, 0.8, 0.1),
            diffuse: 0.7,
            specular: 0.3,
            ..Material::default()
        })
        .as_shape();

    let world = World::new(
        vec![floor, middle, right, left],
        vec![PointLight::new(
            Point::new(-10.0, 10.0, -10.0),
            Color::WHITE,
        )],
    );

    let camera = Camera::new(2048, 1080, FRAC_PI_3).with_transform(Matrix::view_transform(
        Point::new(0.0, 1.5, -5.0),
        Point::new(0.0, 1.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
    ));

    camera.render(&world).write_ppm(path);
}
