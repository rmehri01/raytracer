use std::f64::consts::FRAC_PI_3;

use raytracer::{
    core::{matrix::Matrix, point::Point, vector::Vector},
    graphics::{color::Color, pattern::Pattern},
    raytracer::{
        camera::Camera,
        material::Material,
        point_light::PointLight,
        shapes::{Compound, Operation, Primitive, SetProperties},
        world::World,
    },
};

fn main() {
    render_csg("images/csg.ppm");
}

fn render_csg(path: &str) {
    let floor_pattern = Pattern::new_checker(
        Pattern::new_solid(Color::WHITE),
        Pattern::new_solid(Color::BLACK),
    );
    let floor_material = Material {
        color: Color::new(1.0, 0.9, 0.9),
        specular: 0.0,
        pattern: Some(floor_pattern),
        ..Material::default()
    };

    let floor = Primitive::new_plane()
        .with_material(floor_material)
        .as_shape();
    let cube = Primitive::new_cube()
        .with_material(Material {
            color: Color::new(0.0, 0.0, 1.0),
            ..Material::default()
        })
        .as_shape();
    let sphere = Primitive::new_sphere()
        .with_material(Material {
            color: Color::new(1.0, 0.0, 0.0),
            ..Material::default()
        })
        .with_transform(Matrix::translation(0.5, 0.5, -0.5))
        .as_shape();
    let csg = Compound::new_csg(Operation::Difference, cube, sphere).as_shape();

    let world = World::new(
        vec![floor, csg],
        vec![PointLight::new(
            Point::new(-10.0, 10.0, -10.0),
            Color::WHITE,
        )],
    );

    let camera = Camera::new(2000, 1000, FRAC_PI_3).with_transform(Matrix::view_transform(
        Point::new(0.0, 1.5, -5.0),
        Point::new(0.0, 1.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
    ));

    camera.render(&world).write_ppm(path);
}
