#![warn(clippy::use_self, clippy::semicolon_if_nothing_returned)]

use exercises::clock::plot_clock;

mod core {
    pub mod canvas;
    pub mod color;
    pub mod matrix;
    pub mod tuple;
}

mod exercises {
    pub mod clock;
    pub mod projectile;
}

fn main() {
    plot_clock(400, "images/clock.ppm");
}
