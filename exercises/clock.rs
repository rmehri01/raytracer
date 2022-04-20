use std::f64::consts::FRAC_PI_6;

use raytracer::{
    core::{matrix::Matrix, point::Point},
    graphics::{canvas::Canvas, color::Color},
};

fn main() {
    plot_clock(400, "images/clock.ppm");
}

fn plot_clock(side_len: usize, path: &str) {
    let mut canvas = Canvas::empty(side_len, side_len);

    let twelve = Point::new(0.0, 0.0, 1.0);
    let angle = FRAC_PI_6;
    let pixel_color = Color::WHITE;
    let radius = 3.0 / 8.0 * canvas.width as f64;

    for i in 0..12 {
        let rotation = Matrix::rotation_y(angle * i as f64);
        let point = rotation * twelve;
        let x = point.x * radius + canvas.width as f64 / 2.0;
        let z = point.z * radius + canvas.height as f64 / 2.0;
        canvas.write_pixel(x as usize, z as usize, pixel_color);
    }

    canvas.write_ppm(path);
}
