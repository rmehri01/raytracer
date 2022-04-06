#![warn(clippy::use_self, clippy::semicolon_if_nothing_returned)]

use exercises::sphere::plot_sphere;

mod core {
    pub mod matrix;
    pub mod tuple;
}

mod exercises {
    pub mod clock;
    pub mod projectile;
    pub mod sphere;
}

mod graphics {
    pub mod canvas;
    pub mod color;
}

mod raytracer {
    pub mod intersection;
    pub mod material;
    pub mod object;
    pub mod point_light;
    pub mod ray;
}

fn main() {
    plot_sphere("images/sphere.ppm");
}
