use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_6};

use raytracer::{
    core::{matrix::Matrix, point::Point, vector::Vector},
    graphics::{color::Color, pattern::Pattern},
    raytracer::{
        camera::Camera,
        material::Material,
        point_light::PointLight,
        shapes::{Compound, Conic, Primitive, SetProperties},
        world::World,
    },
};

fn main() {
    render_hexagon("images/ppm/hexagon.ppm");
}

fn render_hexagon(path: &str) {
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
        .to_shape();
    let hexagon = hexagon()
        .with_transform(
            Matrix::translation(0.0, 1.0, 0.0)
                * Matrix::rotation_x(-FRAC_PI_2)
                * Matrix::scaling(0.75, 0.75, 0.75),
        )
        .to_shape();

    let world = World::new(
        vec![floor, hexagon],
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

fn hexagon() -> Compound {
    let mut hex = Vec::new();

    for n in 0..6 {
        let side = hexagon_side()
            .with_transform(Matrix::rotation_y(n as f64 * FRAC_PI_3))
            .to_shape();
        hex.push(side);
    }

    Compound::new_group(hex)
}

fn hexagon_side() -> Compound {
    Compound::new_group(vec![hexagon_corner().to_shape(), hexagon_edge().to_shape()])
}

fn hexagon_corner() -> Primitive {
    Primitive::new_sphere()
        .with_transform(Matrix::translation(0.0, 0.0, -1.0) * Matrix::scaling(0.25, 0.25, 0.25))
}

fn hexagon_edge() -> Primitive {
    Primitive::new_cylinder(Conic::new(0.0, 1.0, false)).with_transform(
        Matrix::translation(0.0, 0.0, -1.0)
            * Matrix::rotation_y(-FRAC_PI_6)
            * Matrix::rotation_z(-FRAC_PI_2)
            * Matrix::scaling(0.25, 1.0, 0.25),
    )
}
