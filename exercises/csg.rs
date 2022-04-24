use std::f64::consts::FRAC_PI_3;

use raytracer::{
    core::{matrix::Matrix, point::Point, vector::Vector},
    graphics::{color::Color, pattern::Pattern},
    raytracer::{
        camera::Camera,
        material::Material,
        point_light::PointLight,
        shapes::{Compound, Operation, SetProperties, Single},
        world::World,
    },
};

fn main() {
    render_csg("images/csg.ppm");
}

fn render_csg(path: &str) {
    let floor_material = Material {
        color: Color::new(1.0, 0.9, 0.9),
        specular: 0.0,
        pattern: Some(Pattern::new_checker(Color::WHITE, Color::BLACK)),
        ..Material::default()
    };

    let floor = Single::new_plane().with_material(floor_material).as_shape();
    let cube = Single::new_cube()
        .with_material(Material {
            color: Color::new(0.0, 0.0, 1.0),
            ..Material::default()
        })
        .as_shape();
    let sphere = Single::new_sphere()
        .with_material(Material {
            color: Color::new(1.0, 0.0, 0.0),
            ..Material::default()
        })
        .with_transform(Matrix::translation(0.5, 0.5, -0.5))
        .as_shape();
    let csg = Compound::new_csg(Operation::Difference, cube, sphere).as_shape();

    let world = World {
        light: Some(PointLight::new(
            Point::new(-10.0, 10.0, -10.0),
            Color::WHITE,
        )),
        shapes: vec![floor, csg],
    };

    let camera = Camera::new(2000, 1000, FRAC_PI_3).with_transform(Matrix::view_transform(
        Point::new(0.0, 1.5, -5.0),
        Point::new(0.0, 1.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
    ));

    camera.render(&world).write_ppm(path);
}
