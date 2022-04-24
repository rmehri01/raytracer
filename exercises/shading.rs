use raytracer::{
    core::point::Point,
    graphics::{canvas::Canvas, color::Color},
    raytracer::{
        intersection::Intersection,
        material::Material,
        point_light::PointLight,
        ray::Ray,
        shapes::{Intersect, SetProperties, Single},
    },
};

fn main() {
    render_shaded_sphere("images/shading.ppm");
}

fn render_shaded_sphere(path: &str) {
    let side_len = 500;
    let mut canvas = Canvas::empty(side_len, side_len);

    let ray_origin = Point::new(0.0, 0.0, -5.0);
    let wall_z = 10.0;
    let wall_size = 7.0;

    let pixel_size = wall_size / side_len as f64;
    let half = wall_size / 2.0;

    let sphere = Single::new_sphere().with_material(Material {
        color: Color::new(1.0, 0.2, 1.0),
        ..Material::default()
    });

    let light_position = Point::new(-10.0, 10.0, -10.0);
    let light_color = Color::WHITE;
    let light = PointLight::new(light_position, light_color);

    for y in 0..side_len {
        let world_y = half - pixel_size * y as f64;

        for x in 0..side_len {
            let world_x = -half + pixel_size * x as f64;
            let position = Point::new(world_x, world_y, wall_z);

            let r = Ray::new(ray_origin, (position - ray_origin).normalize());
            let xs = sphere.intersect(&r, &im_rc::Vector::new());

            if let Some(hit) = xs.hit() {
                let point = r.position(hit.t);
                let normal = sphere.normal_at(
                    &point,
                    &Intersection::new(0.0, &sphere, im_rc::Vector::new()),
                );
                let eye = -r.direction;
                let object_point = sphere.world_to_object(&point, &im_rc::Vector::new());
                let color = hit.shape.properties.material.lighting(
                    &object_point,
                    &point,
                    &light,
                    &eye,
                    &normal,
                    false,
                );
                canvas.write_pixel(x, y, color);
            }
        }
    }

    canvas.write_ppm(path);
}
