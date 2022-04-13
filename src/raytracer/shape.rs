use std::collections::BTreeSet;

use approx::AbsDiffEq;

use crate::core::{matrix::Matrix, tuple::Tuple};

use super::{
    intersection::{Intersection, Intersections},
    material::Material,
    ray::Ray,
};

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Shape {
    // TODO: better way of setting/using this
    pub transform: Matrix<4>,
    pub material: Material,
    kind: ShapeKind,
}

impl Shape {
    pub fn new(kind: ShapeKind) -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
            kind,
        }
    }

    pub fn new_sphere() -> Self {
        Self::new(ShapeKind::Sphere)
    }

    pub fn new_glass_sphere() -> Self {
        let mut sphere = Self::new_sphere();
        sphere.material.transparency = 1.0;
        sphere.material.refractive_index = 1.5;

        sphere
    }

    pub fn new_plane() -> Self {
        Self::new(ShapeKind::Plane)
    }

    pub fn intersect(&self, ray: &Ray) -> Intersections {
        let local_ray = ray.transform(&self.transform.inverse());
        let ts = self.kind.intersect(&local_ray);
        let intersections = ts.iter().map(|t| Intersection::new(*t, *self));

        Intersections(BTreeSet::from_iter(intersections))
    }

    pub fn normal_at(&self, point: &Tuple) -> Tuple {
        let local_point = self.transform.inverse() * *point;
        let local_normal = self.kind.normal_at(&local_point);

        let mut world_normal = self.transform.inverse().transpose() * local_normal;
        world_normal.w = 0.0;

        world_normal.normalize()
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum ShapeKind {
    Sphere,
    Plane,
}

impl ShapeKind {
    pub fn intersect(&self, ray: &Ray) -> Vec<f64> {
        match self {
            Self::Sphere => Self::sphere_intersect(ray),
            Self::Plane => Self::plane_intersect(ray),
        }
    }

    fn sphere_intersect(ray: &Ray) -> Vec<f64> {
        let sphere_to_ray = ray.origin - Tuple::point(0.0, 0.0, 0.0);

        let a = ray.direction.dot(&ray.direction);
        let b = 2.0 * sphere_to_ray.dot(&ray.direction);
        let c = sphere_to_ray.dot(&sphere_to_ray) - 1.0;

        let discriminant = b * b - 4.0 * a * c;

        if discriminant < 0.0 {
            Vec::new()
        } else {
            let t1 = (-b - discriminant.sqrt()) / (2.0 * a);
            let t2 = (-b + discriminant.sqrt()) / (2.0 * a);

            vec![t1, t2]
        }
    }

    fn plane_intersect(ray: &Ray) -> Vec<f64> {
        if ray.direction.y.abs() < Tuple::default_epsilon() {
            Vec::new()
        } else {
            vec![-ray.origin.y / ray.direction.y]
        }
    }

    pub fn normal_at(&self, point: &Tuple) -> Tuple {
        match self {
            Self::Sphere => Self::sphere_normal_at(point),
            Self::Plane => Self::plane_normal_at(),
        }
    }

    fn sphere_normal_at(object_point: &Tuple) -> Tuple {
        *object_point - Tuple::point(0.0, 0.0, 0.0)
    }

    fn plane_normal_at() -> Tuple {
        Tuple::vector(0.0, 1.0, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_4};

    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn default_transformation() {
        let shape = Shape::new_sphere();

        assert_eq!(shape.transform, Matrix::identity());
    }

    #[test]
    fn change_transformation() {
        let mut shape = Shape::new_sphere();
        shape.transform = Matrix::translation(2.0, 3.0, 4.0);

        assert_eq!(shape.transform, Matrix::translation(2.0, 3.0, 4.0));
    }

    #[test]
    fn default_material() {
        let shape = Shape::new_sphere();

        assert_eq!(shape.material, Material::default());
    }

    #[test]
    fn assign_material() {
        let material = Material {
            ambient: 1.0,
            ..Material::default()
        };

        let mut shape = Shape::new_sphere();
        shape.material = material;

        assert_eq!(shape.material, material);
    }

    #[test]
    fn glass_sphere() {
        let shape = Shape::new_glass_sphere();

        assert_eq!(shape.transform, Matrix::identity());
        assert_eq!(shape.material.transparency, 1.0);
        assert_eq!(shape.material.refractive_index, 1.5);
    }

    #[test]
    fn intersect_scaled_sphere_with_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut shape = Shape::new_sphere();
        shape.transform = Matrix::scaling(2.0, 2.0, 2.0);

        let xs = shape
            .intersect(&r)
            .0
            .iter()
            .map(|i| i.t)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], 3.0);
        assert_abs_diff_eq!(xs[1], 7.0);
    }

    #[test]
    fn intersect_translated_sphere_with_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut shape = Shape::new_sphere();
        shape.transform = Matrix::translation(5.0, 0.0, 0.0);

        assert!(shape.intersect(&r).0.is_empty());
    }

    #[test]
    fn ray_intersects_sphere() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = ShapeKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], 4.0);
        assert_abs_diff_eq!(xs[1], 6.0);
    }

    #[test]
    fn ray_tangent_to_sphere() {
        let r = Ray::new(Tuple::point(0.0, 1.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = ShapeKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], 5.0);
        assert_abs_diff_eq!(xs[1], 5.0);
    }

    #[test]
    fn ray_misses_sphere() {
        let r = Ray::new(Tuple::point(0.0, 2.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = ShapeKind::Sphere;

        let intersects = s.intersect(&r);

        assert!(intersects.is_empty());
    }

    #[test]
    fn ray_originates_inside_sphere() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = ShapeKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], -1.0);
        assert_abs_diff_eq!(xs[1], 1.0);
    }

    #[test]
    fn sphere_behind_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, 5.0), Tuple::vector(0.0, 0.0, 1.0));
        let s = ShapeKind::Sphere;

        let xs = s.intersect(&r);

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], -6.0);
        assert_abs_diff_eq!(xs[1], -4.0);
    }

    #[test]
    fn intersecting_a_ray_parallel_to_the_plane() {
        let plane = ShapeKind::Plane;
        let r = Ray::new(Tuple::point(0.0, 10.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = plane.intersect(&r);

        assert!(xs.is_empty());
    }

    #[test]
    fn intersecting_a_coplanar_ray() {
        let plane = ShapeKind::Plane;
        let r = Ray::new(Tuple::point(0.0, 0.0, 0.0), Tuple::vector(0.0, 0.0, 1.0));

        let xs = plane.intersect(&r);

        assert!(xs.is_empty());
    }

    #[test]
    fn ray_intersecting_from_above() {
        let plane = ShapeKind::Plane;
        let r = Ray::new(Tuple::point(0.0, 1.0, 0.0), Tuple::vector(0.0, -1.0, 0.0));

        let xs = plane.intersect(&r);

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 1.0);
    }

    #[test]
    fn ray_intersecting_from_below() {
        let plane = ShapeKind::Plane;
        let r = Ray::new(Tuple::point(0.0, -1.0, 0.0), Tuple::vector(0.0, 1.0, 0.0));

        let xs = plane.intersect(&r);

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 1.0);
    }

    #[test]
    fn intersect_sets_shape() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let shape = Shape::new_sphere();

        let shapes = shape
            .intersect(&r)
            .0
            .iter()
            .map(|i| i.shape)
            .collect::<Vec<_>>();

        assert_eq!(shapes[0], shape);
        assert_eq!(shapes[1], shape);
    }

    #[test]
    fn normal_on_translated_sphere() {
        let mut shape = Shape::new_sphere();
        shape.transform = Matrix::translation(0.0, 1.0, 0.0);

        let n = shape.normal_at(&Tuple::point(0.0, 1.70711, -FRAC_1_SQRT_2));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, FRAC_1_SQRT_2, -FRAC_1_SQRT_2));
    }

    #[test]
    fn normal_on_transformed_sphere() {
        let mut shape = Shape::new_sphere();
        shape.transform = Matrix::scaling(1.0, 0.5, 1.0) * Matrix::rotation_z(FRAC_PI_4);

        let n = shape.normal_at(&Tuple::point(
            0.0,
            2.0_f64.sqrt() / 2.0,
            -(2.0_f64.sqrt()) / 2.0,
        ));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 0.97014, -0.24254));
    }

    #[test]
    fn sphere_normal_x_axis() {
        let s = ShapeKind::Sphere;
        let n = s.normal_at(&Tuple::point(1.0, 0.0, 0.0));

        assert_abs_diff_eq!(n, Tuple::vector(1.0, 0.0, 0.0));
    }

    #[test]
    fn sphere_normal_y_axis() {
        let s = ShapeKind::Sphere;
        let n = s.normal_at(&Tuple::point(0.0, 1.0, 0.0));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 1.0, 0.0));
    }

    #[test]
    fn sphere_normal_z_axis() {
        let s = ShapeKind::Sphere;
        let n = s.normal_at(&Tuple::point(0.0, 0.0, 1.0));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 0.0, 1.0));
    }

    #[test]
    fn sphere_normal_non_axial_point() {
        let s = ShapeKind::Sphere;
        let n = s.normal_at(&Tuple::point(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
        ));

        assert_abs_diff_eq!(
            n,
            Tuple::vector(
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0
            )
        );
    }

    #[test]
    fn normal_is_normalized() {
        let s = ShapeKind::Sphere;
        let n = s.normal_at(&Tuple::point(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
        ));

        assert_abs_diff_eq!(n, n.normalize());
    }

    #[test]
    fn normal_of_plane_is_constant() {
        let plane = ShapeKind::Plane;
        let n1 = plane.normal_at(&Tuple::point(0.0, 0.0, 0.0));
        let n2 = plane.normal_at(&Tuple::point(10.0, 0.0, -10.0));
        let n3 = plane.normal_at(&Tuple::point(-5.0, 0.0, 150.0));

        assert_abs_diff_eq!(n1, Tuple::vector(0.0, 1.0, 0.0));
        assert_abs_diff_eq!(n2, Tuple::vector(0.0, 1.0, 0.0));
        assert_abs_diff_eq!(n3, Tuple::vector(0.0, 1.0, 0.0));
    }
}