use raytracer::{
    core::tuple::Tuple,
    graphics::{canvas::Canvas, color::Color},
};

struct Projectile {
    position: Tuple,
    velocity: Tuple,
}

impl Projectile {
    fn new(position: Tuple, velocity: Tuple) -> Self {
        Self { position, velocity }
    }
}

struct Environment {
    gravity: Tuple,
    wind: Tuple,
}

impl Environment {
    fn new(gravity: Tuple, wind: Tuple) -> Self {
        Self { gravity, wind }
    }
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

fn plot_trajectory(mut p: Projectile, e: &Environment, path: &str) {
    let mut canvas = Canvas::new(900, 550);
    let pixel_color = Color::new(0.85, 0.35, 0.40);

    while p.position.y > 0.0 {
        canvas.write_pixel(
            p.position.x.round() as usize,
            (550_f64 - p.position.y).round() as usize,
            pixel_color,
        );
        p = tick(p, e);
    }

    canvas.write_ppm(path);
}

fn tick(proj: Projectile, env: &Environment) -> Projectile {
    Projectile {
        position: proj.position + proj.velocity,
        velocity: proj.velocity + env.gravity + env.wind,
    }
}
