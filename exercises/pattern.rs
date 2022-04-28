use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_4};

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
    render_pattern("images/pattern.ppm");
}

fn render_pattern(path: &str) {
    let sub_pattern = Pattern::new_checker(
        Pattern::new_solid(Color::new(0.20, 0.20, 0.20)),
        Pattern::new_solid(Color::new(0.55, 0.55, 0.55)),
    )
    .with_transform(Matrix::scaling(0.25, 0.25, 0.25));
    let floor_pattern = Pattern::new_stripe(
        sub_pattern,
        Pattern::new_solid(Color::new(0.20, 0.20, 0.20)),
    )
    .with_transform(Matrix::rotation_y(FRAC_PI_3) * Matrix::scaling(0.5, 0.5, 0.5));
    let floor_material = Material {
        specular: 0.0,
        pattern: Some(floor_pattern),
        ..Material::default()
    };
    let floor = Primitive::new_plane()
        .with_material(floor_material.clone())
        .as_shape();

    let left_wall = Primitive::new_plane()
        .with_transform(
            Matrix::translation(0.0, 0.0, 5.0)
                * Matrix::rotation_y(-FRAC_PI_4)
                * Matrix::rotation_x(FRAC_PI_2),
        )
        .with_material(floor_material.clone())
        .as_shape();

    let right_wall = Primitive::new_plane()
        .with_transform(
            Matrix::translation(0.0, 0.0, 5.0)
                * Matrix::rotation_y(FRAC_PI_4)
                * Matrix::rotation_x(FRAC_PI_2),
        )
        .with_material(floor_material)
        .as_shape();

    let middle_pattern = Pattern::new_perturb(Pattern::new_ring(
        Pattern::new_solid(Color::new(0.0, 0.3, 0.6)),
        Pattern::new_solid(Color::new(0.1, 1.0, 0.8)),
    ))
    .with_transform(Matrix::rotation_x(-FRAC_PI_3) * Matrix::scaling(0.2, 0.2, 0.2));
    let middle = Primitive::new_sphere()
        .with_transform(Matrix::translation(-0.5, 1.0, 0.5))
        .with_material(Material {
            diffuse: 0.7,
            specular: 0.3,
            pattern: Some(middle_pattern),
            ..Material::default()
        })
        .as_shape();

    let right_pattern = Pattern::new_stripe(
        Pattern::new_solid(Color::BLACK),
        Pattern::new_solid(Color::WHITE),
    )
    .with_transform(Matrix::rotation_z(-FRAC_PI_4) * Matrix::scaling(0.2, 0.2, 0.2));
    let right = Primitive::new_sphere()
        .with_transform(Matrix::translation(1.5, 0.5, -0.5) * Matrix::scaling(0.5, 0.5, 0.5))
        .with_material(Material {
            diffuse: 0.7,
            specular: 0.3,
            pattern: Some(right_pattern),
            ..Material::default()
        })
        .as_shape();

    let left_pattern = Pattern::new_gradient(
        Pattern::new_solid(Color::new(1.0, 0.0, 0.0)),
        Pattern::new_solid(Color::new(0.0, 1.0, 1.0)),
    )
    .with_transform(Matrix::translation(-1.0, 0.0, 0.0) * Matrix::scaling(2.0, 2.0, 2.0));
    let left = Primitive::new_sphere()
        .with_transform(Matrix::translation(-1.5, 0.33, -0.75) * Matrix::scaling(0.33, 0.33, 0.33))
        .with_material(Material {
            diffuse: 0.7,
            specular: 0.3,
            pattern: Some(left_pattern),
            ..Material::default()
        })
        .as_shape();

    let world = World::new(
        vec![floor, left_wall, right_wall, middle, left, right],
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
