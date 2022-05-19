use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

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
    render_cover("images/ppm/cover.ppm");
}

fn render_cover(path: &str) {
    let white_material = Material {
        color: Color::WHITE,
        diffuse: 0.7,
        ambient: 0.1,
        specular: 0.0,
        reflective: 0.1,
        ..Material::default()
    };

    let blue_material = Material {
        color: Color::new(0.537, 0.831, 0.914),
        ..white_material.clone()
    };

    let red_material = Material {
        color: Color::new(0.941, 0.322, 0.388),
        ..white_material.clone()
    };

    let purple_material = Material {
        color: Color::new(0.373, 0.404, 0.550),
        ..white_material.clone()
    };

    let standard_transform = Matrix::scaling(0.5, 0.5, 0.5) * Matrix::translation(1.0, -1.0, 1.0);

    let large_object = Matrix::scaling(3.5, 3.5, 3.5) * standard_transform;
    let medium_object = Matrix::scaling(3.0, 3.0, 3.0) * standard_transform;
    let small_object = Matrix::scaling(2.0, 2.0, 2.0) * standard_transform;

    let backdrop = Primitive::new_plane()
        .with_material(Material {
            color: Color::WHITE,
            ambient: 1.0,
            diffuse: 0.0,
            specular: 0.0,
            ..Material::default()
        })
        .with_transform(Matrix::translation(0.0, 0.0, 500.0) * Matrix::rotation_x(FRAC_PI_2))
        .to_shape();

    let main_sphere = Primitive::new_sphere()
        .with_material(Material {
            color: Color::new(0.373, 0.404, 0.550),
            diffuse: 0.2,
            ambient: 0.0,
            specular: 1.0,
            shininess: 200.0,
            reflective: 0.7,
            transparency: 0.7,
            refractive_index: 1.5,
            ..Material::default()
        })
        .with_transform(large_object)
        .to_shape();

    let cube_specs = [
        (
            white_material.clone(),
            Matrix::translation(4.0, 0.0, 0.0),
            medium_object,
        ),
        (
            blue_material.clone(),
            Matrix::translation(8.5, 1.5, -0.5),
            large_object,
        ),
        (
            red_material.clone(),
            Matrix::translation(0.0, 0.0, 4.0),
            large_object,
        ),
        (
            white_material.clone(),
            Matrix::translation(4.0, 0.0, 4.0),
            small_object,
        ),
        (
            purple_material.clone(),
            Matrix::translation(7.5, 0.5, 4.0),
            medium_object,
        ),
        (
            white_material.clone(),
            Matrix::translation(-0.25, 0.25, 8.0),
            medium_object,
        ),
        (
            blue_material.clone(),
            Matrix::translation(4.0, 1.0, 7.5),
            large_object,
        ),
        (
            red_material.clone(),
            Matrix::translation(10.0, 2.0, 7.5),
            medium_object,
        ),
        (
            white_material.clone(),
            Matrix::translation(8.0, 2.0, 12.0),
            small_object,
        ),
        (
            white_material.clone(),
            Matrix::translation(20.0, 1.0, 9.0),
            small_object,
        ),
        (
            blue_material,
            Matrix::translation(-0.5, -5.0, 0.25),
            large_object,
        ),
        (
            red_material,
            Matrix::translation(4.0, -4.0, 0.0),
            large_object,
        ),
        (
            white_material.clone(),
            Matrix::translation(8.5, -4.0, 0.0),
            large_object,
        ),
        (
            white_material.clone(),
            Matrix::translation(0.0, -4.0, 4.0),
            large_object,
        ),
        (
            purple_material,
            Matrix::translation(-0.5, -4.5, 8.0),
            large_object,
        ),
        (
            white_material.clone(),
            Matrix::translation(0.0, -8.0, 4.0),
            large_object,
        ),
        (
            white_material,
            Matrix::translation(-0.5, -8.5, 8.0),
            large_object,
        ),
    ];

    let mut shapes = vec![backdrop, main_sphere];
    for (material, transform, base_transform) in cube_specs {
        shapes.push(
            Primitive::new_cube()
                .with_material(material)
                .with_transform(transform * base_transform)
                .to_shape(),
        );
    }

    let world = World::new(
        shapes,
        vec![
            PointLight::new(Point::new(50.0, 100.0, -50.0), Color::WHITE),
            PointLight::new(Point::new(-400.0, 50.0, -10.0), Color::new(0.7, 0.7, 0.7)),
        ],
    );

    let camera = Camera::new(2048, 1080, FRAC_PI_4).with_transform(Matrix::view_transform(
        Point::new(-6.0, 6.0, -10.0),
        Point::new(6.0, -2.0, 6.0),
        Vector::new(-0.45, 1.0, 0.0),
    ));

    camera.render(&world).write_ppm(path);
}
