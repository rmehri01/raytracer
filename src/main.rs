#![warn(clippy::use_self, clippy::semicolon_if_nothing_returned)]

use exercises::projectile::{plot_trajectory, Environment, Projectile};

use crate::core::tuple::Tuple;

mod core {
    pub mod canvas;
    pub mod color;
    pub mod matrix;
    pub mod tuple;
}

mod exercises {
    pub mod projectile;
}

fn main() {
    let start = Tuple::point(0.0, 1.0, 0.0);
    let velocity = Tuple::vector(1.0, 1.8, 0.0).normalize() * 11.25;
    let p = Projectile::new(start, velocity);

    let gravity = Tuple::vector(0.0, -0.1, 0.0);
    let wind = Tuple::vector(-0.01, 0.0, 0.0);
    let e = Environment::new(gravity, wind);

    plot_trajectory(p, &e, "images/trajectory.ppm");
}
