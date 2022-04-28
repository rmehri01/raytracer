use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_6};

use raytracer::{
    core::{
        matrix::{Matrix, Transformation},
        point::Point,
        vector::Vector,
    },
    graphics::{color::Color, pattern::Pattern},
    raytracer::{
        camera::Camera,
        material::Material,
        point_light::PointLight,
        shapes::{Compound, Conic, Operation, Primitive, SetProperties, Shape},
        world::World,
    },
};

fn main() {
    render_csg("images/ppm/csg.ppm");
}

fn render_csg(path: &str) {
    let red_material = Material {
        color: Color::new(1.0, 0.0, 0.0),
        ambient: 0.2,
        ..Material::default()
    };
    let green_material = Material {
        color: Color::new(0.0, 1.0, 0.0),
        ambient: 0.2,
        ..Material::default()
    };
    let blue_material = Material {
        color: Color::new(0.0, 0.0, 1.0),
        ambient: 0.2,
        ..Material::default()
    };

    let dark_mirror = Material {
        color: Color::BLACK,
        ambient: 0.0,
        diffuse: 0.4,
        reflective: 0.5,
        ..Material::default()
    };
    let transparent = Material {
        color: Color::BLACK,
        ambient: 0.0,
        diffuse: 0.0,
        reflective: 0.0,
        transparency: 1.0,
        refractive_index: 1.0,
        ..Material::default()
    };

    let room_pattern = Pattern::new_checker(
        Pattern::new_solid(Color::WHITE),
        Pattern::new_solid(Color::new(0.9, 0.9, 0.9)),
    )
    .with_transform(Matrix::scaling(0.05, 0.05, 0.05));
    let room = Primitive::new_cube()
        .with_material(Material {
            ambient: 0.1,
            diffuse: 0.7,
            reflective: 0.05,
            pattern: Some(room_pattern),
            ..Material::default()
        })
        .with_transform(Matrix::scaling(5.0, 5.0, 5.0) * Matrix::translation(0.0, 1.0, 0.0))
        .as_shape();

    let left = Primitive::new_cylinder(Conic::new(-1.0, 1.0, true))
        .with_material(Material {
            ambient: 0.1,
            diffuse: 0.5,
            reflective: 0.3,
            ..red_material.clone()
        })
        .with_transform(Matrix::scaling(0.5, 1.1, 0.5))
        .as_shape();
    let right = Compound::new_csg(
        Operation::Intersection,
        Primitive::new_cylinder(Conic::new(-1.0, 1.0, true))
            .with_material(Material {
                ambient: 0.1,
                diffuse: 0.5,
                reflective: 0.3,
                ..green_material.clone()
            })
            .with_transform(Matrix::rotation_x(FRAC_PI_2) * Matrix::scaling(0.5, 1.1, 0.5))
            .as_shape(),
        Primitive::new_cylinder(Conic::new(-1.0, 1.0, true))
            .with_material(Material {
                ambient: 0.1,
                diffuse: 0.5,
                reflective: 0.3,
                ..blue_material.clone()
            })
            .with_transform(Matrix::rotation_z(FRAC_PI_2) * Matrix::scaling(0.5, 1.1, 0.5))
            .as_shape(),
    )
    .as_shape();
    let tricylinder = Compound::new_csg(Operation::Intersection, left, right)
        .with_transform(
            Matrix::translation(-1.5, 0.7, 0.0)
                * Matrix::rotation_z(-0.2)
                * Matrix::rotation_x(-0.1)
                * Matrix::rotation_y(0.4),
        )
        .as_shape();

    let sphere = Primitive::new_sphere()
        .with_material(Material {
            color: Color::new(0.1, 0.1, 0.1),
            ambient: 0.2,
            diffuse: 0.9,
            specular: 1.0,
            shininess: 50.0,
            ..Material::default()
        })
        .with_transform(Matrix::scaling(1.4, 1.4, 1.4))
        .as_shape();
    let cylinders = Compound::new_group(vec![
        Primitive::new_cylinder(Conic::new(-1.0, 1.0, true))
            .with_material(red_material.clone())
            .with_transform(Matrix::scaling(0.5, 1.1, 0.5))
            .as_shape(),
        Primitive::new_cylinder(Conic::new(-1.0, 1.0, true))
            .with_material(green_material)
            .with_transform(Matrix::rotation_x(FRAC_PI_2) * Matrix::scaling(0.5, 1.1, 0.5))
            .as_shape(),
        Primitive::new_cylinder(Conic::new(-1.0, 1.0, true))
            .with_material(blue_material)
            .with_transform(Matrix::rotation_z(FRAC_PI_2) * Matrix::scaling(0.5, 1.1, 0.5))
            .as_shape(),
    ])
    .as_shape();
    let cube_minus_cylinders = Compound::new_csg(
        Operation::Difference,
        Primitive::new_cube().with_material(dark_mirror).as_shape(),
        cylinders,
    )
    .as_shape();
    let hollowed_box = Compound::new_csg(Operation::Intersection, sphere, cube_minus_cylinders)
        .with_transform(
            Matrix::rotation_y(1.3)
                * Matrix::scaling(0.5, 0.5, 0.5)
                * Matrix::translation(0.0, 1.0, 0.0),
        )
        .as_shape();

    let inside = Primitive::new_sphere()
        .with_material(red_material)
        .as_shape();
    let children = (0..12)
        .map(|i| {
            let angle = (i as f64) * FRAC_PI_6;
            wedge(Matrix::rotation_y(angle))
        })
        .collect::<Vec<_>>();
    let outside = Compound::new_group(children)
        .with_material(transparent)
        .as_shape();
    let ball = Compound::new_csg(Operation::Intersection, inside, outside)
        .with_transform(
            Matrix::translation(1.5, 0.25, 0.0)
                * Matrix::rotation_z(0.1)
                * Matrix::rotation_x(-0.1)
                * Matrix::rotation_y(-0.5)
                * Matrix::scaling(0.5, 0.5, 0.5)
                * Matrix::translation(0.0, 1.0, 0.0),
        )
        .as_shape();

    let world = World::new(
        vec![room, tricylinder, hollowed_box, ball],
        vec![PointLight::new(Point::new(-2.0, 5.0, -2.0), Color::WHITE)],
    );

    let camera = Camera::new(2048, 1080, 0.9).with_transform(Matrix::view_transform(
        Point::new(0.0, 2.0, -4.9),
        Point::new(0.0, 0.5, 0.0),
        Vector::new(0.0, 1.0, 0.0),
    ));

    camera.render(&world).write_ppm(path);
}

fn wedge(transform: Transformation) -> Shape {
    Primitive::new_cube()
        .with_transform(
            transform
                * Matrix::scaling(1.0, 1.0, 0.15)
                * Matrix::translation(2.0_f64.sqrt(), 0.0, 0.0)
                * Matrix::rotation_y(FRAC_PI_4),
        )
        .with_shadow(false)
        .as_shape()
}
