use std::f64::consts::FRAC_PI_3;

use raytracer::{
    core::{matrix::Matrix, tuple::Tuple},
    graphics::{color::Color, pattern::Pattern},
    raytracer::{
        camera::Camera,
        material::Material,
        point_light::PointLight,
        shape::{Operation, Shape},
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

    let floor = Shape::new_plane().with_material(floor_material);
    let cube = Shape::new_cube().with_material(Material {
        color: Color::new(0.0, 0.0, 1.0),
        ..Material::default()
    });
    let sphere = Shape::new_sphere()
        .with_material(Material {
            color: Color::new(1.0, 0.0, 0.0),
            ..Material::default()
        })
        .with_transform(Matrix::translation(0.5, 0.5, -0.5));
    let csg = Shape::new_csg(Operation::Difference, cube, sphere);

    let world = World {
        light: Some(PointLight::new(
            Tuple::point(-10.0, 10.0, -10.0),
            Color::WHITE,
        )),
        shapes: vec![floor, csg],
    };

    let camera = Camera::new(2000, 1000, FRAC_PI_3).with_transform(Matrix::view_transform(
        Tuple::point(0.0, 1.5, -5.0),
        Tuple::point(0.0, 1.0, 0.0),
        Tuple::vector(0.0, 1.0, 0.0),
    ));

    camera.render(&world).write_ppm(path);
}
