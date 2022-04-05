#![warn(clippy::use_self, clippy::semicolon_if_nothing_returned)]

use exercises::sphere::plot_sphere;

mod core {
    pub mod canvas;
    pub mod color;
    pub mod intersection;
    pub mod matrix;
    pub mod object;
    pub mod ray;
    pub mod tuple;
}

mod exercises {
    pub mod clock;
    pub mod projectile;
    pub mod sphere;
}

fn main() {
    plot_sphere("images/sphere.ppm");
}
