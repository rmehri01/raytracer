#![allow(dead_code)]

use crate::core::{canvas::Canvas, color::Color, tuple::Tuple};

pub struct Projectile {
    position: Tuple,
    velocity: Tuple,
}

impl Projectile {
    pub fn new(position: Tuple, velocity: Tuple) -> Self {
        Self { position, velocity }
    }
}

pub struct Environment {
    gravity: Tuple,
    wind: Tuple,
}

impl Environment {
    pub fn new(gravity: Tuple, wind: Tuple) -> Self {
        Self { gravity, wind }
    }
}

pub fn plot_trajectory(mut p: Projectile, e: &Environment, path: &str) {
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

pub fn print_trajectory(mut p: Projectile, e: &Environment) {
    let mut ticks = 0;
    while p.position.y > 0.0 {
        p = tick(p, e);
        ticks += 1;
        println!("{:#?}", p.position);
        println!("{}", ticks);
    }

    println!("final position {:#?}", p.position);
    println!("total {} ticks", ticks);
}

fn tick(proj: Projectile, env: &Environment) -> Projectile {
    Projectile {
        position: proj.position + proj.velocity,
        velocity: proj.velocity + env.gravity + env.wind,
    }
}
