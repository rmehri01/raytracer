#![warn(
    clippy::use_self,
    clippy::semicolon_if_nothing_returned,
    clippy::needless_pass_by_value,
    clippy::inconsistent_struct_constructor,
    clippy::trivially_copy_pass_by_ref,
    clippy::float_cmp,
    clippy::match_same_arms
)]

pub mod core {
    pub mod matrix;
    pub mod point;
    pub mod vector;
}

pub mod graphics {
    pub mod canvas;
    pub mod color;
    pub mod pattern;
}

pub mod io {
    pub mod obj;
}

pub mod raytracer {
    mod bounds;
    pub mod camera;
    pub mod intersection;
    pub mod material;
    pub mod point_light;
    pub mod ray;
    pub mod world;

    pub mod shapes {
        pub use compound::*;
        pub use shape::*;
        pub use single::*;

        mod compound;
        mod shape;
        mod single;
    }
}
