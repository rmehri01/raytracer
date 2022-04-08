use crate::{
    core::{matrix::Matrix, tuple::Tuple},
    graphics::canvas::Canvas,
};

use super::{ray::Ray, world::World};

pub struct Camera {
    hsize: usize,
    vsize: usize,
    field_of_view: f64,
    transform: Matrix<4>,
    half_width: f64,
    half_height: f64,
    pixel_size: f64,
}

impl Camera {
    pub fn new(hsize: usize, vsize: usize, field_of_view: f64) -> Self {
        let half_view = (field_of_view / 2.0).tan();
        let aspect = hsize as f64 / vsize as f64;

        let half_width;
        let half_height;

        if aspect >= 1.0 {
            half_width = half_view;
            half_height = half_view / aspect;
        } else {
            half_width = half_view * aspect;
            half_height = half_view;
        }

        Self {
            hsize,
            vsize,
            field_of_view,
            transform: Matrix::identity(),
            half_width,
            half_height,
            pixel_size: (half_width * 2.0) / hsize as f64,
        }
    }

    pub fn render(&self, world: World) -> Canvas {
        let mut canvas = Canvas::new(self.hsize, self.vsize);

        for y in 0..self.vsize {
            for x in 0..self.hsize {
                let ray = self.ray_for_pixel(x, y);
                let color = world.color_at(&ray);
                canvas.write_pixel(x, y, color);
            }
        }

        canvas
    }

    fn ray_for_pixel(&self, px: usize, py: usize) -> Ray {
        let x_offset = (px as f64 + 0.5) * self.pixel_size;
        let y_offset = (py as f64 + 0.5) * self.pixel_size;

        let world_x = self.half_width - x_offset;
        let world_y = self.half_height - y_offset;

        let pixel = self.transform.inverse() * Tuple::point(world_x, world_y, -1.0);
        let origin = self.transform.inverse() * Tuple::point(0.0, 0.0, 0.0);
        let direction = (pixel - origin).normalize();

        Ray::new(origin, direction)
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

    use approx::assert_abs_diff_eq;

    use crate::graphics::color::Color;

    use super::*;

    #[test]
    fn test_camera() {
        let c = Camera::new(160, 120, FRAC_PI_2);

        assert_eq!(c.hsize, 160);
        assert_eq!(c.vsize, 120);
        assert_abs_diff_eq!(c.field_of_view, FRAC_PI_2);
        assert_eq!(c.transform, Matrix::identity());
    }

    #[test]
    fn pixel_size_horizontal_canvas() {
        let c = Camera::new(200, 125, FRAC_PI_2);

        assert_abs_diff_eq!(c.pixel_size, 0.01);
    }

    #[test]
    fn pixel_size_vertical_canvas() {
        let c = Camera::new(125, 200, FRAC_PI_2);

        assert_abs_diff_eq!(c.pixel_size, 0.01);
    }

    #[test]
    fn ray_through_center_of_canvas() {
        let c = Camera::new(201, 101, FRAC_PI_2);

        let r = c.ray_for_pixel(100, 50);

        assert_abs_diff_eq!(r.origin, Tuple::point(0.0, 0.0, 0.0));
        assert_abs_diff_eq!(r.direction, Tuple::vector(0.0, 0.0, -1.0));
    }

    #[test]
    fn ray_through_corner_of_canvas() {
        let c = Camera::new(201, 101, FRAC_PI_2);

        let r = c.ray_for_pixel(0, 0);

        assert_abs_diff_eq!(r.origin, Tuple::point(0.0, 0.0, 0.0));
        assert_abs_diff_eq!(r.direction, Tuple::vector(0.66519, 0.33259, -0.66851));
    }

    #[test]
    fn ray_when_camera_is_transformed() {
        let mut c = Camera::new(201, 101, FRAC_PI_2);
        c.transform = Matrix::rotation_y(FRAC_PI_4) * Matrix::translation(0.0, -2.0, 5.0);

        let r = c.ray_for_pixel(100, 50);

        assert_abs_diff_eq!(r.origin, Tuple::point(0.0, 2.0, -5.0));
        assert_abs_diff_eq!(
            r.direction,
            Tuple::vector(2.0_f64.sqrt() / 2.0, 0.0, -(2.0_f64.sqrt()) / 2.0)
        );
    }

    #[test]
    fn render_world_with_camera() {
        let w = World::default();
        let mut c = Camera::new(11, 11, FRAC_PI_2);

        let from = Tuple::point(0.0, 0.0, -5.0);
        let to = Tuple::point(0.0, 0.0, 0.0);
        let up = Tuple::vector(0.0, 1.0, 0.0);

        c.transform = Matrix::view_transform(from, to, up);

        let image = c.render(w);

        assert_abs_diff_eq!(image.pixel_at(5, 5), Color::new(0.38066, 0.47583, 0.2855));
    }
}