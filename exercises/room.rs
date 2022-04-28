use std::f64::consts::FRAC_PI_4;

use raytracer::{
    core::{matrix::Matrix, point::Point, vector::Vector},
    graphics::{color::Color, pattern::Pattern},
    raytracer::{
        camera::Camera,
        material::Material,
        point_light::PointLight,
        shapes::{Primitive, SetProperties, Shape},
        world::World,
    },
};

fn main() {
    render_room("images/room.ppm");
}

fn render_room(path: &str) {
    let floor_ceiling_pattern = Pattern::new_checker(
        Pattern::new_solid(Color::BLACK),
        Pattern::new_solid(Color::new(0.25, 0.25, 0.25)),
    )
    .with_transform(Matrix::scaling(0.07, 0.07, 0.07));
    let floor_ceiling = Primitive::new_cube()
        .with_material(Material {
            ambient: 0.25,
            diffuse: 0.75,
            specular: 0.9,
            shininess: 300.0,
            reflective: 0.1,
            pattern: Some(floor_ceiling_pattern),
            ..Material::default()
        })
        .with_transform(Matrix::scaling(20.0, 7.0, 20.0) * Matrix::translation(0.0, 1.0, 0.0))
        .as_shape();

    let walls_pattern = Pattern::new_checker(
        Pattern::new_solid(Color::new(0.4863, 0.3765, 0.2941)),
        Pattern::new_solid(Color::new(0.3725, 0.2902, 0.2275)),
    )
    .with_transform(Matrix::scaling(0.05, 20.0, 0.05));
    let walls = Primitive::new_cube()
        .with_material(Material {
            ambient: 0.1,
            diffuse: 0.7,
            specular: 0.9,
            shininess: 300.0,
            reflective: 0.1,
            pattern: Some(walls_pattern),
            ..Material::default()
        })
        .with_transform(Matrix::scaling(10.0, 10.0, 10.0))
        .as_shape();

    let table_top_pattern = Pattern::new_stripe(
        Pattern::new_solid(Color::new(0.5529, 0.4235, 0.3255)),
        Pattern::new_solid(Color::new(0.6588, 0.5098, 0.4000)),
    )
    .with_transform(Matrix::scaling(0.05, 0.05, 0.05) * Matrix::rotation_y(0.1));
    let table_top = Primitive::new_cube()
        .with_material(Material {
            ambient: 0.1,
            diffuse: 0.7,
            specular: 0.9,
            shininess: 300.0,
            reflective: 0.2,
            pattern: Some(table_top_pattern),
            ..Material::default()
        })
        .with_transform(Matrix::translation(0.0, 3.1, 0.0) * Matrix::scaling(3.0, 0.1, 2.0))
        .as_shape();

    let leg1 = make_table_leg(2.7, -1.7);
    let leg2 = make_table_leg(2.7, 1.7);
    let leg3 = make_table_leg(-2.7, -1.7);
    let leg4 = make_table_leg(-2.7, 1.7);

    let glass_cube = Primitive::new_cube()
        .with_material(Material {
            color: Color::new(1.0, 1.0, 0.8),
            ambient: 0.0,
            diffuse: 0.3,
            specular: 0.9,
            shininess: 300.0,
            reflective: 0.7,
            transparency: 0.7,
            refractive_index: 1.5,
            ..Material::default()
        })
        .with_transform(
            Matrix::translation(0.0, 3.450001, 0.0)
                * Matrix::rotation_y(0.2)
                * Matrix::scaling(0.25, 0.25, 0.25),
        )
        .as_shape();

    let little_cube1 = Primitive::new_cube()
        .with_material(Material {
            color: Color::new(1.0, 0.5, 0.5),
            diffuse: 0.4,
            reflective: 0.6,
            ..Material::default()
        })
        .with_transform(
            Matrix::translation(1.0, 3.35, -0.9)
                * Matrix::rotation_y(-0.4)
                * Matrix::scaling(0.15, 0.15, 0.15),
        )
        .as_shape();

    let little_cube2 = Primitive::new_cube()
        .with_material(Material {
            color: Color::new(1.0, 0.5, 0.5),
            ..Material::default()
        })
        .with_transform(
            Matrix::translation(-1.5, 3.27, 0.3)
                * Matrix::rotation_y(0.4)
                * Matrix::scaling(0.15, 0.17, 0.15),
        )
        .as_shape();

    let little_cube3 = Primitive::new_cube()
        .with_material(Material {
            color: Color::new(0.5, 1.0, 0.5),
            ..Material::default()
        })
        .with_transform(
            Matrix::translation(0.0, 3.25, 1.0)
                * Matrix::rotation_y(0.4)
                * Matrix::scaling(0.2, 0.05, 0.05),
        )
        .as_shape();

    let little_cube4 = Primitive::new_cube()
        .with_material(Material {
            color: Color::new(0.5, 0.5, 1.0),
            ..Material::default()
        })
        .with_transform(
            Matrix::translation(-0.6, 3.4, -1.0)
                * Matrix::rotation_y(0.8)
                * Matrix::scaling(0.05, 0.2, 0.05),
        )
        .as_shape();

    let little_cube5 = Primitive::new_cube()
        .with_material(Material {
            color: Color::new(0.5, 1.0, 1.0),
            ..Material::default()
        })
        .with_transform(
            Matrix::translation(2.0, 3.4, 1.0)
                * Matrix::rotation_y(0.8)
                * Matrix::scaling(0.05, 0.2, 0.05),
        )
        .as_shape();

    let frame1 = Primitive::new_cube()
        .with_material(Material {
            color: Color::new(0.7098, 0.2471, 0.2196),
            diffuse: 0.6,
            ..Material::default()
        })
        .with_transform(Matrix::translation(-10.0, 4.0, 1.0) * Matrix::scaling(0.05, 1.0, 1.0))
        .as_shape();

    let frame2 = Primitive::new_cube()
        .with_material(Material {
            color: Color::new(0.2667, 0.2706, 0.6902),
            diffuse: 0.6,
            ..Material::default()
        })
        .with_transform(Matrix::translation(-10.0, 3.4, 2.7) * Matrix::scaling(0.05, 0.4, 0.4))
        .as_shape();

    let frame3 = Primitive::new_cube()
        .with_material(Material {
            color: Color::new(0.3098, 0.5961, 0.3098),
            diffuse: 0.6,
            ..Material::default()
        })
        .with_transform(Matrix::translation(-10.0, 4.6, 2.7) * Matrix::scaling(0.05, 0.4, 0.4))
        .as_shape();

    let mirror_frame = Primitive::new_cube()
        .with_material(Material {
            color: Color::new(0.3882, 0.2627, 0.1882),
            diffuse: 0.7,
            ..Material::default()
        })
        .with_transform(Matrix::translation(-2.0, 3.5, 9.95) * Matrix::scaling(5.0, 1.5, 0.05))
        .as_shape();

    let mirror = Primitive::new_cube()
        .with_material(Material {
            color: Color::BLACK,
            ambient: 0.0,
            diffuse: 0.0,
            specular: 1.0,
            shininess: 300.0,
            reflective: 1.0,
            ..Material::default()
        })
        .with_transform(Matrix::translation(-2.0, 3.5, 9.95) * Matrix::scaling(4.8, 1.4, 0.06))
        .as_shape();

    let world = World::new(
        vec![
            floor_ceiling,
            walls,
            table_top,
            leg1,
            leg2,
            leg3,
            leg4,
            glass_cube,
            little_cube1,
            little_cube2,
            little_cube3,
            little_cube4,
            little_cube5,
            frame1,
            frame2,
            frame3,
            mirror_frame,
            mirror,
        ],
        vec![PointLight::new(Point::new(0.0, 6.9, -5.0), Color::WHITE)],
    );

    let camera = Camera::new(2048, 1080, FRAC_PI_4).with_transform(Matrix::view_transform(
        Point::new(8.0, 6.0, -8.0),
        Point::new(0.0, 3.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
    ));

    camera.render(&world).write_ppm(path);
}

fn make_table_leg(x: f64, z: f64) -> Shape {
    Primitive::new_cube()
        .with_material(Material {
            color: Color::new(0.5529, 0.4235, 0.3255),
            ambient: 0.2,
            diffuse: 0.7,
            ..Material::default()
        })
        .with_transform(Matrix::translation(x, 1.5, z) * Matrix::scaling(0.1, 1.5, 0.1))
        .as_shape()
}
