use exercises::projectile::{compute_trajectory, Environment, Projectile};

use crate::core::tuple::Tuple;

mod core {
    pub mod canvas;
    pub mod color;
    pub mod tuple;
}

mod exercises {
    pub mod projectile;
}

fn main() {
    let p = Projectile::new(
        Tuple::point(0.0, 1.0, 0.0),
        Tuple::vector(1.0, 1.0, 0.0).normalize(),
    );

    let e = Environment::new(
        Tuple::vector(0.0, -0.1, 0.0),
        Tuple::vector(-0.01, 0.0, 0.0),
    );

    compute_trajectory(p, &e);
}
