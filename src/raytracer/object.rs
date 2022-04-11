use std::collections::BTreeSet;

use crate::core::{matrix::Matrix, tuple::Tuple};

use super::{
    intersection::{Intersection, Intersections},
    material::Material,
    ray::Ray,
    shapes::{plane::Plane, shape::Shape, sphere::Sphere},
};

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Object {
    // TODO: better way of setting/using this
    pub transform: Matrix<4>,
    pub material: Material,
    shape: Shape,
}

impl Object {
    pub fn new(shape: Shape) -> Self {
        Self {
            transform: Matrix::identity(),
            material: Material::default(),
            shape,
        }
    }

    pub fn new_sphere() -> Self {
        Self::new(Shape::Sphere(Sphere::default()))
    }

    pub fn new_plane() -> Self {
        Self::new(Shape::Plane(Plane::default()))
    }

    pub fn intersect(&self, ray: &Ray) -> Intersections {
        let local_ray = ray.transform(&self.transform.inverse());
        let ts = self.shape.intersect(&local_ray);
        let intersections = ts.iter().map(|t| Intersection::new(*t, *self));

        Intersections(BTreeSet::from_iter(intersections))
    }

    pub fn normal_at(&self, point: &Tuple) -> Tuple {
        let local_point = self.transform.inverse() * *point;
        let local_normal = self.shape.normal_at(&local_point);

        let mut world_normal = self.transform.inverse().transpose() * local_normal;
        world_normal.w = 0.0;

        world_normal.normalize()
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_4};

    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn default_transformation() {
        let obj = Object::new_sphere();

        assert_eq!(obj.transform, Matrix::identity());
    }

    #[test]
    fn change_transformation() {
        let mut obj = Object::new_sphere();
        obj.transform = Matrix::translation(2.0, 3.0, 4.0);

        assert_eq!(obj.transform, Matrix::translation(2.0, 3.0, 4.0));
    }

    #[test]
    fn default_material() {
        let obj = Object::new_sphere();

        assert_eq!(obj.material, Material::default());
    }

    #[test]
    fn assign_material() {
        let material = Material {
            ambient: 1.0,
            ..Material::default()
        };

        let mut obj = Object::new_sphere();
        obj.material = material;

        assert_eq!(obj.material, material);
    }

    #[test]
    fn intersect_scaled_sphere_with_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut obj = Object::new_sphere();
        obj.transform = Matrix::scaling(2.0, 2.0, 2.0);

        let xs = obj.intersect(&r).0.iter().map(|i| i.t).collect::<Vec<_>>();

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], 3.0);
        assert_abs_diff_eq!(xs[1], 7.0);
    }

    #[test]
    fn intersect_translated_sphere_with_ray() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let mut obj = Object::new_sphere();
        obj.transform = Matrix::translation(5.0, 0.0, 0.0);

        assert!(obj.intersect(&r).0.is_empty());
    }

    #[test]
    fn normal_on_translated_sphere() {
        let mut obj = Object::new_sphere();
        obj.transform = Matrix::translation(0.0, 1.0, 0.0);

        let n = obj.normal_at(&Tuple::point(0.0, 1.70711, -FRAC_1_SQRT_2));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, FRAC_1_SQRT_2, -FRAC_1_SQRT_2));
    }

    #[test]
    fn normal_on_transformed_sphere() {
        let mut obj = Object::new_sphere();
        obj.transform = Matrix::scaling(1.0, 0.5, 1.0) * Matrix::rotation_z(FRAC_PI_4);

        let n = obj.normal_at(&Tuple::point(
            0.0,
            2.0_f64.sqrt() / 2.0,
            -(2.0_f64.sqrt()) / 2.0,
        ));

        assert_abs_diff_eq!(n, Tuple::vector(0.0, 0.97014, -0.24254));
    }

    #[test]
    fn intersect_sets_object() {
        let r = Ray::new(Tuple::point(0.0, 0.0, -5.0), Tuple::vector(0.0, 0.0, 1.0));
        let obj = Object::new_sphere();

        let objects = obj
            .intersect(&r)
            .0
            .iter()
            .map(|i| i.object)
            .collect::<Vec<_>>();

        assert_eq!(objects[0], obj);
        assert_eq!(objects[1], obj);
    }
}
