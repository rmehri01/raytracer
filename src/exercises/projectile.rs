use crate::core::tuple::Tuple;

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

pub fn compute_trajectory(mut p: Projectile, e: &Environment) {
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
