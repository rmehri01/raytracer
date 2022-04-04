#![allow(dead_code)]

use std::f64::consts::PI;

use crate::core::{canvas::Canvas, color::Color, matrix::Matrix, tuple::Tuple};

pub fn plot_clock(side_len: usize, path: &str) {
    let mut canvas = Canvas::new(side_len, side_len);

    let twelve = Tuple::point(0.0, 0.0, 1.0);
    let angle = PI / 6.0;
    let pixel_color = Color::new(1.0, 1.0, 1.0);
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
