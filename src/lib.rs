#![warn(clippy::use_self, clippy::semicolon_if_nothing_returned)]

pub mod core {
    pub mod matrix;
    pub mod tuple;
}

pub mod graphics {
    pub mod canvas;
    pub mod color;
}

pub mod raytracer {
    pub mod camera;
    pub mod intersection;
    pub mod material;
    pub mod object;
    pub mod point_light;
    pub mod ray;
    pub mod world;

    pub mod shapes {
        pub mod plane;
        pub mod shape;
        pub mod sphere;
    }
}
