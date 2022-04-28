use std::f64::consts::FRAC_PI_3;

use raytracer::{
    core::{matrix::Matrix, point::Point, vector::Vector},
    graphics::{color::Color, pattern::Pattern},
    io::obj::{self, ParseError},
    raytracer::{
        camera::Camera,
        material::Material,
        point_light::PointLight,
        shapes::{Primitive, SetProperties},
        world::World,
    },
};

fn main() -> Result<(), ParseError> {
    render_teapot("images/teapot.ppm")
}

fn render_teapot(path: &str) -> Result<(), ParseError> {
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
    let teapot = obj::parse_file("exercises/resources/teapot.obj")?;

    let world = World::new(
        vec![floor, teapot],
        vec![PointLight::new(
            Point::new(-10.0, 10.0, -10.0),
            Color::WHITE,
        )],
    );

    let camera = Camera::new(500, 500, FRAC_PI_3).with_transform(Matrix::view_transform(
        Point::new(0.0, 1.5, -5.0),
        Point::new(0.0, 1.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
    ));

    camera.render(&world).write_ppm(path);

    Ok(())
}
