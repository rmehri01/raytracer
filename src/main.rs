use tuple::Tuple;

mod canvas;
mod color;
mod tuple;

fn main() {
    let mut p = Projectile {
        position: Tuple::point(0.0, 1.0, 0.0),
        velocity: Tuple::vector(1.0, 1.0, 0.0).normalize(),
    };

    let e = Environment {
        gravity: Tuple::vector(0.0, -0.1, 0.0),
        wind: Tuple::vector(-0.01, 0.0, 0.0),
    };

    let mut ticks = 0;
    while p.position.y > 0.0 {
        p = tick(p, &e);
        ticks += 1;
        println!("{:#?}", p.position);
        println!("{}", ticks);
    }

    println!("final position {:#?}", p.position);
    println!("total {} ticks", ticks);
}

struct Projectile {
    position: Tuple,
    velocity: Tuple,
}

struct Environment {
    gravity: Tuple,
    wind: Tuple,
}

fn tick(proj: Projectile, env: &Environment) -> Projectile {
    Projectile {
        position: proj.position + proj.velocity,
        velocity: proj.velocity + env.gravity + env.wind,
    }
}
