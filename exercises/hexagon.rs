use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_6};

use raytracer::{
    core::{matrix::Matrix, tuple::Tuple},
    graphics::{color::Color, pattern::Pattern},
    raytracer::{
        camera::Camera,
        material::Material,
        point_light::PointLight,
        shape::{Conic, Shape},
        world::World,
    },
};

fn main() {
    render_hexagon("images/hexagon.ppm");
}

fn render_hexagon(path: &str) {
    let floor_material = Material {
        color: Color::new(1.0, 0.9, 0.9),
        specular: 0.0,
        pattern: Some(Pattern::new_checker(Color::WHITE, Color::BLACK)),
        ..Material::default()
    };

    let floor = Shape::new_plane().with_material(floor_material);
    let hexagon = hexagon().with_transform(
        Matrix::translation(0.0, 1.0, 0.0)
            * Matrix::rotation_x(-FRAC_PI_2)
            * Matrix::scaling(0.75, 0.75, 0.75),
    );

    let world = World {
        light: Some(PointLight::new(
            Tuple::point(-10.0, 10.0, -10.0),
            Color::WHITE,
        )),
        shapes: vec![floor, hexagon],
    };

    let camera = Camera::new(2000, 1000, FRAC_PI_3).with_transform(Matrix::view_transform(
        Tuple::point(0.0, 1.5, -5.0),
        Tuple::point(0.0, 1.0, 0.0),
        Tuple::vector(0.0, 1.0, 0.0),
    ));

    camera.render(world).write_ppm(path);
}

fn hexagon() -> Shape {
    let mut hex = Vec::new();

    for n in 0..6 {
        let side = hexagon_side().with_transform(Matrix::rotation_y(n as f64 * FRAC_PI_3));
        hex.push(side);
    }

    Shape::new_group(hex)
}

fn hexagon_side() -> Shape {
    Shape::new_group(vec![hexagon_corner(), hexagon_edge()])
}

fn hexagon_corner() -> Shape {
    Shape::new_sphere()
        .with_transform(Matrix::translation(0.0, 0.0, -1.0) * Matrix::scaling(0.25, 0.25, 0.25))
}

fn hexagon_edge() -> Shape {
    Shape::new_cylinder(Conic::new(0.0, 1.0, false)).with_transform(
        Matrix::translation(0.0, 0.0, -1.0)
            * Matrix::rotation_y(-FRAC_PI_6)
            * Matrix::rotation_z(-FRAC_PI_2)
            * Matrix::scaling(0.25, 1.0, 0.25),
    )
}
