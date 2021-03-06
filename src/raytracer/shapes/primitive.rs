use approx::AbsDiffEq;

use crate::{
    core::{matrix::Transformation, point::Point, vector::Vector},
    graphics::color::Color,
    raytracer::{
        bounds::Bounds,
        intersection::{Intersection, Intersections},
        material::Material,
        point_light::PointLight,
        ray::Ray,
    },
};

use super::{HasProperties, Intersect, Properties, SetProperties, Shape};

#[derive(Debug, PartialEq, Clone)]
pub struct Primitive {
    properties: Properties,
    pub(crate) has_shadow: bool,
    kind: Kind,
}

impl Primitive {
    fn new(kind: Kind) -> Self {
        Self {
            properties: Properties {
                bounds: Self::compute_bounds(&kind),
                ..Properties::default()
            },
            has_shadow: true,
            kind,
        }
    }

    fn compute_bounds(kind: &Kind) -> Bounds {
        match kind {
            Kind::Cube | Kind::Sphere => {
                Bounds::new(Point::new(-1.0, -1.0, -1.0), Point::new(1.0, 1.0, 1.0))
            }
            Kind::Plane => Bounds::new(
                Point::new(f64::NEG_INFINITY, 0.0, f64::NEG_INFINITY),
                Point::new(f64::INFINITY, 0.0, f64::INFINITY),
            ),
            Kind::Cylinder(conic) | Kind::Cone(conic) => conic.bounds(),
            Kind::Triangle { triangular, .. } | Kind::SmoothTriangle { triangular, .. } => {
                triangular.bounds()
            }
        }
    }

    pub fn new_sphere() -> Self {
        Self::new(Kind::Sphere)
    }

    pub fn new_glass_sphere() -> Self {
        Self::new_sphere().with_material(Material {
            transparency: 1.0,
            refractive_index: 1.5,
            ..Material::default()
        })
    }

    pub fn new_plane() -> Self {
        Self::new(Kind::Plane)
    }

    pub fn new_cube() -> Self {
        Self::new(Kind::Cube)
    }

    pub fn new_cylinder(conic: Conic) -> Self {
        Self::new(Kind::Cylinder(conic))
    }

    pub fn new_cone(conic: Conic) -> Self {
        Self::new(Kind::Cone(conic))
    }

    pub fn new_triangle(p1: Point, p2: Point, p3: Point) -> Self {
        let triangular = Triangular::new(p1, p2, p3);
        let normal = triangular.e1.cross(&triangular.e2).normalize();

        Self::new(Kind::Triangle { triangular, normal })
    }

    pub fn new_smooth_triangle(
        p1: Point,
        p2: Point,
        p3: Point,
        n1: Vector,
        n2: Vector,
        n3: Vector,
    ) -> Self {
        let triangular = Triangular::new(p1, p2, p3);

        Self::new(Kind::SmoothTriangle {
            triangular,
            n1,
            n2,
            n3,
        })
    }

    pub fn with_shadow(mut self, has_shadow: bool) -> Self {
        self.has_shadow = has_shadow;
        self
    }

    /// Returns the color of the material using the Phuong Reflection Model.
    pub fn lighting(
        &self,
        position: &Point,
        light: &PointLight,
        eye_v: &Vector,
        normal_v: &Vector,
        in_shadow: bool,
        trail: &im_rc::Vector<Transformation>,
    ) -> Color {
        let material = &self.properties.material;
        let color = material.pattern.as_ref().map_or(material.color, |pattern| {
            pattern.pattern_at_shape(self, position, trail)
        });

        let effective_color = color * light.intensity;
        let light_v = (light.position - *position).normalize();
        let ambient = effective_color * material.ambient;

        let diffuse;
        let specular;

        let light_dot_normal = light_v.dot(normal_v);
        if light_dot_normal < 0.0 || in_shadow {
            diffuse = Color::BLACK;
            specular = Color::BLACK;
        } else {
            diffuse = effective_color * material.diffuse * light_dot_normal;

            let reflect_v = (-light_v).reflect(normal_v);
            let reflect_dot_eye = reflect_v.dot(eye_v);

            specular = if reflect_dot_eye <= 0.0 {
                Color::BLACK
            } else {
                let factor = reflect_dot_eye.powf(material.shininess);
                light.intensity * material.specular * factor
            };
        }

        ambient + diffuse + specular
    }

    pub fn normal_at(&self, point: &Point, hit: &Intersection) -> Vector {
        let local_point = self.world_to_object(point, &hit.trail);
        let local_normal = match &self.kind {
            Kind::Sphere => Kind::sphere_normal_at(&local_point),
            Kind::Plane => Kind::plane_normal_at(),
            Kind::Cube => Kind::cube_normal_at(&local_point),
            Kind::Cylinder(conic) => Kind::cylinder_normal_at(&local_point, conic),
            Kind::Cone(conic) => Kind::cone_normal_at(&local_point, conic),
            Kind::Triangle { normal, .. } => *normal,
            Kind::SmoothTriangle { n1, n2, n3, .. } => {
                Kind::smooth_triangle_normal_at(n1, n2, n3, hit)
            }
        };

        self.normal_to_world(&local_normal, &hit.trail)
    }

    pub(crate) fn world_to_object(
        &self,
        world_point: &Point,
        trail: &im_rc::Vector<Transformation>,
    ) -> Point {
        let trail_point = trail.iter().rev().fold(*world_point, |acc, &mat| mat * acc);

        self.properties.inverse_transform * trail_point
    }

    pub(crate) fn normal_to_world(
        &self,
        normal: &Vector,
        trail: &im_rc::Vector<Transformation>,
    ) -> Vector {
        let normal = self.properties.inverse_transform.transpose() * *normal;
        let normal = normal.normalize();

        trail.iter().fold(normal, |acc, mat| {
            let normal = mat.transpose() * acc;
            normal.normalize()
        })
    }

    pub fn to_shape(self) -> Shape {
        Shape::Primitive(self)
    }
}

impl HasProperties for Primitive {
    fn properties(&self) -> &Properties {
        &self.properties
    }

    fn properties_mut(&mut self) -> &mut Properties {
        &mut self.properties
    }
}

impl Intersect for Primitive {
    fn local_intersect(&self, ray: &Ray, trail: &im_rc::Vector<Transformation>) -> Intersections {
        let mut u_v = None;
        let ts = match &self.kind {
            Kind::Sphere => Kind::sphere_intersect(ray),
            Kind::Plane => Kind::plane_intersect(ray),
            Kind::Cube => Kind::cube_intersect(ray),
            Kind::Cylinder(conic) => Kind::cylinder_intersect(ray, conic),
            Kind::Cone(conic) => Kind::cone_intersect(ray, conic),
            Kind::Triangle { triangular, .. } => triangular.intersect(ray),
            Kind::SmoothTriangle { triangular, .. } => {
                let (ts, u, v) = triangular.intersect_uv(ray);

                u_v = Some((u, v));
                ts
            }
        };

        let intersections = ts
            .iter()
            .map(|t| Intersection::new_with_uv(*t, self, u_v, trail.clone()))
            .collect();

        Intersections(intersections)
    }
}

#[derive(Debug, PartialEq, Clone)]
enum Kind {
    /// A unit sphere with center at the origin.
    Sphere,
    /// A perfectly flat surface that extends infinitely in x and z.
    Plane,
    /// An axis-aligned bounding box that is centered at the origin and
    /// extends from -1 to 1 along each axis.
    Cube,
    /// A cylinder that is centered at the origin with radius 1 and extends
    /// from minimum to maximum exclusive along the y axis.
    Cylinder(Conic),
    /// A double-napped cone that is centered at the origin and extends
    /// from minimum to maximum exclusive along the y axis.
    Cone(Conic),
    /// A triangle with vertices at p1, p2, and p3, along with two edge
    /// vectors and a normal vector to optimize intersection calculations.
    Triangle {
        triangular: Triangular,
        normal: Vector,
    },
    /// A triangle with vertices at p1, p2, and p3 and a normal at each of
    /// the vertices. Uses normal interpolation to calculate the normal at
    /// any point on the triangle.
    SmoothTriangle {
        triangular: Triangular,
        n1: Vector,
        n2: Vector,
        n3: Vector,
    },
}

impl Kind {
    fn sphere_intersect(ray: &Ray) -> Vec<f64> {
        let sphere_to_ray = ray.origin - Point::new(0.0, 0.0, 0.0);

        let a = ray.direction.dot(&ray.direction);
        let b = 2.0 * sphere_to_ray.dot(&ray.direction);
        let c = sphere_to_ray.dot(&sphere_to_ray) - 1.0;

        let discriminant = b.powi(2) - 4.0 * a * c;

        if discriminant < 0.0 {
            Vec::new()
        } else {
            let t1 = (-b - discriminant.sqrt()) / (2.0 * a);
            let t2 = (-b + discriminant.sqrt()) / (2.0 * a);

            vec![t1, t2]
        }
    }

    fn plane_intersect(ray: &Ray) -> Vec<f64> {
        if ray.direction.y.abs() < Vector::default_epsilon() {
            Vec::new()
        } else {
            vec![-ray.origin.y / ray.direction.y]
        }
    }

    fn cube_intersect(ray: &Ray) -> Vec<f64> {
        let min = -1.0;
        let max = 1.0;

        let (x_t_min, x_t_max) = check_axis(ray.origin.x, ray.direction.x, min, max);
        let (y_t_min, y_t_max) = check_axis(ray.origin.y, ray.direction.y, min, max);
        let (z_t_min, z_t_max) = check_axis(ray.origin.z, ray.direction.z, min, max);

        let t_min = x_t_min.max(y_t_min).max(z_t_min);
        let t_max = x_t_max.min(y_t_max).min(z_t_max);

        if t_min > t_max {
            Vec::new()
        } else {
            vec![t_min, t_max]
        }
    }

    fn cylinder_intersect(ray: &Ray, conic: &Conic) -> Vec<f64> {
        let a = ray.direction.x.powi(2) + ray.direction.z.powi(2);
        let b = 2.0 * ray.origin.x * ray.direction.x + 2.0 * ray.origin.z * ray.direction.z;
        let c = ray.origin.x.powi(2) + ray.origin.z.powi(2) - 1.0;

        conic.intersect(ray, a, b, c, |_| 1.0)
    }

    fn cone_intersect(ray: &Ray, conic: &Conic) -> Vec<f64> {
        let a = ray.direction.x.powi(2) - ray.direction.y.powi(2) + ray.direction.z.powi(2);
        let b = 2.0 * ray.origin.x * ray.direction.x - 2.0 * ray.origin.y * ray.direction.y
            + 2.0 * ray.origin.z * ray.direction.z;
        let c = ray.origin.x.powi(2) - ray.origin.y.powi(2) + ray.origin.z.powi(2);

        conic.intersect(ray, a, b, c, f64::abs)
    }

    fn sphere_normal_at(object_point: &Point) -> Vector {
        *object_point - Point::new(0.0, 0.0, 0.0)
    }

    fn plane_normal_at() -> Vector {
        Vector::new(0.0, 1.0, 0.0)
    }

    fn cube_normal_at(object_point: &Point) -> Vector {
        let max_c = object_point
            .x
            .abs()
            .max(object_point.y.abs())
            .max(object_point.z.abs());

        if max_c.abs_diff_eq(&object_point.x.abs(), Point::default_epsilon()) {
            Vector::new(object_point.x, 0.0, 0.0)
        } else if max_c.abs_diff_eq(&object_point.y.abs(), Point::default_epsilon()) {
            Vector::new(0.0, object_point.y, 0.0)
        } else {
            Vector::new(0.0, 0.0, object_point.z)
        }
    }

    fn cylinder_normal_at(object_point: &Point, conic: &Conic) -> Vector {
        conic.normal_at(object_point, 0.0)
    }

    fn cone_normal_at(object_point: &Point, conic: &Conic) -> Vector {
        let y = object_point.x.hypot(object_point.z);
        let y = if object_point.y > 0.0 { -y } else { y };

        conic.normal_at(object_point, y)
    }

    fn smooth_triangle_normal_at(
        n1: &Vector,
        n2: &Vector,
        n3: &Vector,
        hit: &Intersection,
    ) -> Vector {
        let (u, v) = hit.u_v.expect("u and v must be set for smooth triangle");

        *n2 * u + *n3 * v + *n1 * (1.0 - u - v)
    }
}

/// Checks if a ray intersects an axis plane and returns the minimum and
/// maximum t values.
pub(crate) fn check_axis(origin: f64, direction: f64, min: f64, max: f64) -> (f64, f64) {
    let t_min_numerator = min - origin;
    let t_max_numerator = max - origin;

    let t_min = t_min_numerator / direction;
    let t_max = t_max_numerator / direction;

    if t_min > t_max {
        (t_max, t_min)
    } else {
        (t_min, t_max)
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
struct Triangular {
    p1: Point,
    p2: Point,
    p3: Point,
    e1: Vector,
    e2: Vector,
}

impl Triangular {
    fn new(p1: Point, p2: Point, p3: Point) -> Self {
        let e1 = p2 - p1;
        let e2 = p3 - p1;

        Self { p1, p2, p3, e1, e2 }
    }

    fn intersect(&self, ray: &Ray) -> Vec<f64> {
        let (ts, _, _) = self.intersect_uv(ray);
        ts
    }

    fn intersect_uv(&self, ray: &Ray) -> (Vec<f64>, f64, f64) {
        let dir_cross_e2 = ray.direction.cross(&self.e2);
        let det = self.e1.dot(&dir_cross_e2);

        let f = 1.0 / det;
        let p1_to_origin = ray.origin - self.p1;
        let u = f * p1_to_origin.dot(&dir_cross_e2);

        let origin_cross_e1 = p1_to_origin.cross(&self.e1);
        let v = f * ray.direction.dot(&origin_cross_e1);

        let ts = if det.abs().abs_diff_eq(&0.0, Point::default_epsilon())
            || !(0.0..=1.0).contains(&u)
            || v < 0.0
            || (u + v) > 1.0
        {
            Vec::new()
        } else {
            let t = f * self.e2.dot(&origin_cross_e1);

            vec![t]
        };

        (ts, u, v)
    }

    fn bounds(&self) -> Bounds {
        let mut bounds = Bounds::default();

        bounds.add_point(&self.p1);
        bounds.add_point(&self.p2);
        bounds.add_point(&self.p3);

        bounds
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Conic {
    minimum: f64,
    maximum: f64,
    closed: bool,
}

impl Conic {
    pub fn new(minimum: f64, maximum: f64, closed: bool) -> Self {
        Self {
            minimum,
            maximum,
            closed,
        }
    }

    fn normal_at(&self, object_point: &Point, y: f64) -> Vector {
        let dist = object_point.x.powi(2) + object_point.z.powi(2);

        if dist < 1.0 && object_point.y >= self.maximum - Point::default_epsilon() {
            Vector::new(0.0, 1.0, 0.0)
        } else if dist < 1.0 && object_point.y <= self.minimum + Point::default_epsilon() {
            Vector::new(0.0, -1.0, 0.0)
        } else {
            Vector::new(object_point.x, y, object_point.z)
        }
    }

    fn intersect(&self, ray: &Ray, a: f64, b: f64, c: f64, radius_at: fn(f64) -> f64) -> Vec<f64> {
        let discriminant = b.powi(2) - 4.0 * a * c;
        let mut xs = self.intersect_caps(ray, radius_at);

        if discriminant >= 0.0 {
            if a.abs_diff_eq(&0.0, Point::default_epsilon()) {
                if !b.abs_diff_eq(&0.0, Point::default_epsilon()) {
                    let t = -c / (2.0 * b);
                    xs.push(t);
                }
            } else {
                let t0 = (-b - discriminant.sqrt()) / (2.0 * a);
                let t1 = (-b + discriminant.sqrt()) / (2.0 * a);

                let y0 = ray.origin.y + t0 * ray.direction.y;
                if self.minimum < y0 && y0 < self.maximum {
                    xs.push(t0);
                }

                let y1 = ray.origin.y + t1 * ray.direction.y;
                if self.minimum < y1 && y1 < self.maximum {
                    xs.push(t1);
                }
            }
        }

        xs
    }

    fn intersect_caps(&self, ray: &Ray, radius_at: fn(f64) -> f64) -> Vec<f64> {
        let mut xs = Vec::new();

        if self.closed && !ray.direction.y.abs_diff_eq(&0.0, Vector::default_epsilon()) {
            let t = (self.minimum - ray.origin.y) / ray.direction.y;
            if Self::check_cap(ray, t, radius_at(self.minimum)) {
                xs.push(t);
            }

            let t = (self.maximum - ray.origin.y) / ray.direction.y;
            if Self::check_cap(ray, t, radius_at(self.maximum)) {
                xs.push(t);
            }
        }

        xs
    }

    fn check_cap(ray: &Ray, t: f64, radius: f64) -> bool {
        let x = ray.origin.x + t * ray.direction.x;
        let z = ray.origin.z + t * ray.direction.z;

        (x.powi(2) + z.powi(2)) <= radius.powi(2)
    }

    fn bounds(&self) -> Bounds {
        Bounds::new(
            Point::new(-1.0, self.minimum, -1.0),
            Point::new(1.0, self.maximum, 1.0),
        )
    }
}

impl Default for Conic {
    fn default() -> Self {
        Self {
            minimum: f64::NEG_INFINITY,
            maximum: f64::INFINITY,
            closed: false,
        }
    }
}

impl PartialEq for Conic {
    fn eq(&self, other: &Self) -> bool {
        self.abs_diff_eq(other, Self::default_epsilon())
    }
}

impl AbsDiffEq for Conic {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-5
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.minimum.abs_diff_eq(&other.minimum, epsilon)
            && self.maximum.abs_diff_eq(&other.maximum, epsilon)
            && self.closed == other.closed
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, FRAC_PI_4};

    use approx::assert_abs_diff_eq;

    use crate::{core::matrix::Matrix, graphics::pattern::Pattern, raytracer::shapes::Compound};

    use super::*;

    #[test]
    fn glass_sphere() {
        let shape = Primitive::new_glass_sphere();

        assert_eq!(shape.properties.transform, Matrix::identity());
        assert_abs_diff_eq!(shape.properties.material.transparency, 1.0);
        assert_abs_diff_eq!(shape.properties.material.refractive_index, 1.5);
    }

    #[test]
    fn constructing_a_triangle() {
        let p1 = Point::new(0.0, 1.0, 0.0);
        let p2 = Point::new(-1.0, 0.0, 0.0);
        let p3 = Point::new(1.0, 0.0, 0.0);
        let triangle = Primitive::new_triangle(p1, p2, p3);

        if let Kind::Triangle { triangular, normal } = triangle.kind {
            assert_eq!(p1, triangular.p1);
            assert_eq!(p2, triangular.p2);
            assert_eq!(p3, triangular.p3);
            assert_eq!(triangular.e1, Vector::new(-1.0, -1.0, 0.0));
            assert_eq!(triangular.e2, Vector::new(1.0, -1.0, 0.0));
            assert_eq!(normal, Vector::new(0.0, 0.0, 1.0));
        } else {
            panic!("expected a triangle");
        }
    }

    #[test]
    fn intersect_scaled_sphere_with_ray() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let shape = Primitive::new_sphere().with_transform(Matrix::scaling(2.0, 2.0, 2.0));

        let xs = shape
            .intersect(&r, &im_rc::Vector::new())
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
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let shape = Primitive::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));

        assert!(shape.intersect(&r, &im_rc::Vector::new()).0.is_empty());
    }

    #[test]
    fn ray_intersects_sphere() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let s = Primitive::new_sphere();

        let xs = s
            .intersect(&r, &im_rc::Vector::new())
            .0
            .iter()
            .map(|i| i.t)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], 4.0);
        assert_abs_diff_eq!(xs[1], 6.0);
    }

    #[test]
    fn ray_tangent_to_sphere() {
        let r = Ray::new(Point::new(0.0, 1.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let s = Primitive::new_sphere();

        let xs = s
            .intersect(&r, &im_rc::Vector::new())
            .0
            .iter()
            .map(|i| i.t)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 5.0);
    }

    #[test]
    fn ray_misses_sphere() {
        let r = Ray::new(Point::new(0.0, 2.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let s = Primitive::new_sphere();

        let intersects = s.intersect(&r, &im_rc::Vector::new());

        assert!(intersects.0.is_empty());
    }

    #[test]
    fn ray_originates_inside_sphere() {
        let r = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));
        let s = Primitive::new_sphere();

        let xs = s
            .intersect(&r, &im_rc::Vector::new())
            .0
            .iter()
            .map(|i| i.t)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], -1.0);
        assert_abs_diff_eq!(xs[1], 1.0);
    }

    #[test]
    fn sphere_behind_ray() {
        let r = Ray::new(Point::new(0.0, 0.0, 5.0), Vector::new(0.0, 0.0, 1.0));
        let s = Primitive::new_sphere();

        let xs = s
            .intersect(&r, &im_rc::Vector::new())
            .0
            .iter()
            .map(|i| i.t)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 2);
        assert_abs_diff_eq!(xs[0], -6.0);
        assert_abs_diff_eq!(xs[1], -4.0);
    }

    #[test]
    fn intersecting_a_ray_parallel_to_the_plane() {
        let plane = Primitive::new_plane();
        let r = Ray::new(Point::new(0.0, 10.0, 0.0), Vector::new(0.0, 0.0, 1.0));

        let xs = plane.intersect(&r, &im_rc::Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn intersecting_a_coplanar_ray() {
        let plane = Primitive::new_plane();
        let r = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));

        let xs = plane.intersect(&r, &im_rc::Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_intersecting_from_above() {
        let plane = Primitive::new_plane();
        let r = Ray::new(Point::new(0.0, 1.0, 0.0), Vector::new(0.0, -1.0, 0.0));

        let xs = plane
            .intersect(&r, &im_rc::Vector::new())
            .0
            .iter()
            .map(|i| i.t)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 1.0);
    }

    #[test]
    fn ray_intersecting_from_below() {
        let plane = Primitive::new_plane();
        let r = Ray::new(Point::new(0.0, -1.0, 0.0), Vector::new(0.0, 1.0, 0.0));

        let xs = plane
            .intersect(&r, &im_rc::Vector::new())
            .0
            .iter()
            .map(|i| i.t)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 1.0);
    }

    #[test]
    fn intersect_sets_shape() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let shape = Primitive::new_sphere();

        let shape_clone = shape.clone();
        let shapes = shape_clone
            .intersect(&r, &im_rc::Vector::new())
            .0
            .iter()
            .map(|i| i.shape)
            .collect::<Vec<_>>();

        assert_eq!(shapes[0], &shape);
        assert_eq!(shapes[1], &shape);
    }

    #[test]
    fn ray_intersects_cube() {
        let scenarios = [
            (
                Point::new(5.0, 0.5, 0.0),
                Vector::new(-1.0, 0.0, 0.0),
                4.0,
                6.0,
            ),
            (
                Point::new(-5.0, 0.5, 0.0),
                Vector::new(1.0, 0.0, 0.0),
                4.0,
                6.0,
            ),
            (
                Point::new(0.5, 5.0, 0.0),
                Vector::new(0.0, -1.0, 0.0),
                4.0,
                6.0,
            ),
            (
                Point::new(0.5, -5.0, 0.0),
                Vector::new(0.0, 1.0, 0.0),
                4.0,
                6.0,
            ),
            (
                Point::new(0.5, 0.0, 5.0),
                Vector::new(0.0, 0.0, -1.0),
                4.0,
                6.0,
            ),
            (
                Point::new(0.5, 0.0, -5.0),
                Vector::new(0.0, 0.0, 1.0),
                4.0,
                6.0,
            ),
            (
                Point::new(0.0, 0.5, 0.0),
                Vector::new(0.0, 0.0, 1.0),
                -1.0,
                1.0,
            ),
        ];

        for (origin, direction, t0, t1) in scenarios {
            let c = Primitive::new_cube();
            let r = Ray::new(origin, direction);

            let xs = c
                .intersect(&r, &im_rc::Vector::new())
                .0
                .iter()
                .map(|i| i.t)
                .collect::<Vec<_>>();

            assert_eq!(xs.len(), 2);
            assert_abs_diff_eq!(xs[0], t0);
            assert_abs_diff_eq!(xs[1], t1);
        }
    }

    #[test]
    fn ray_misses_cube() {
        let scenarios = [
            (
                Point::new(-2.0, 0.0, 0.0),
                Vector::new(0.2673, 0.5345, 0.8018),
            ),
            (
                Point::new(0.0, -2.0, 0.0),
                Vector::new(0.8018, 0.2673, 0.5345),
            ),
            (
                Point::new(0.0, 0.0, -2.0),
                Vector::new(0.5345, 0.8018, 0.2673),
            ),
            (Point::new(2.0, 0.0, 2.0), Vector::new(0.0, 0.0, -1.0)),
            (Point::new(0.0, 2.0, 2.0), Vector::new(0.0, -1.0, 0.0)),
            (Point::new(2.0, 2.0, 0.0), Vector::new(-1.0, 0.0, 0.0)),
        ];

        for (origin, direction) in scenarios {
            let c = Primitive::new_cube();
            let r = Ray::new(origin, direction);

            let xs = c.intersect(&r, &im_rc::Vector::new());

            assert!(xs.0.is_empty());
        }
    }

    #[test]
    fn ray_misses_cylinder() {
        let scenarios = [
            (Point::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
            (Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
            (Point::new(0.0, 0.0, -5.0), Vector::new(1.0, 1.0, 1.0)),
        ];

        for (origin, direction) in scenarios {
            let c = Primitive::new_cylinder(Conic::default());
            let r = Ray::new(origin, direction);

            let xs = c.intersect(&r, &im_rc::Vector::new());

            assert!(xs.0.is_empty());
        }
    }

    #[test]
    fn ray_intersects_cylinder() {
        let scenarios = [
            (
                Point::new(1.0, 0.0, -5.0),
                Vector::new(0.0, 0.0, 1.0),
                vec![5.0],
            ),
            (
                Point::new(0.0, 0.0, -5.0),
                Vector::new(0.0, 0.0, 1.0),
                vec![4.0, 6.0],
            ),
            (
                Point::new(0.5, 0.0, -5.0),
                Vector::new(0.1, 1.0, 1.0),
                vec![6.80798, 7.08872],
            ),
        ];

        for (origin, direction, ts) in scenarios {
            let c = Primitive::new_cylinder(Conic::default());
            let r = Ray::new(origin, direction.normalize());

            let xs = c
                .intersect(&r, &im_rc::Vector::new())
                .0
                .iter()
                .map(|i| i.t)
                .collect::<Vec<_>>();

            assert_eq!(xs.len(), ts.len());
            for (x_t, t) in xs.iter().zip(ts.iter()) {
                assert_abs_diff_eq!(x_t, t, epsilon = Point::default_epsilon());
            }
        }
    }

    #[test]
    fn intersecting_constrained_cylinder() {
        let scenarios = [
            (Point::new(0.0, 1.5, 0.0), Vector::new(0.1, 1.0, 0.0), 0),
            (Point::new(0.0, 3.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
            (Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
            (Point::new(0.0, 2.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
            (Point::new(0.0, 1.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
            (Point::new(0.0, 1.5, -2.0), Vector::new(0.0, 0.0, 1.0), 2),
        ];

        for (point, direction, count) in scenarios {
            let cyl = Primitive::new_cylinder(Conic::new(1.0, 2.0, false));

            let r = Ray::new(point, direction.normalize());
            let xs = cyl.intersect(&r, &im_rc::Vector::new());

            assert_eq!(xs.0.len(), count);
        }
    }

    #[test]
    fn intersecting_caps_of_closed_cylinder() {
        let scenarios = [
            (Point::new(0.0, 3.0, 0.0), Vector::new(0.0, -1.0, 0.0), 2),
            (Point::new(0.0, 3.0, -2.0), Vector::new(0.0, -1.0, 2.0), 2),
            (Point::new(0.0, 4.0, -2.0), Vector::new(0.0, -1.0, 1.0), 2),
            (Point::new(0.0, 0.0, -2.0), Vector::new(0.0, 1.0, 2.0), 2),
            (Point::new(0.0, -1.0, -2.0), Vector::new(0.0, 1.0, 1.0), 2),
        ];

        for (point, direction, count) in scenarios {
            let cyl = Primitive::new_cylinder(Conic::new(1.0, 2.0, true));

            let r = Ray::new(point, direction.normalize());
            let xs = cyl.intersect(&r, &im_rc::Vector::new());

            assert_eq!(xs.0.len(), count);
        }
    }

    #[test]
    fn intersect_cone_with_ray() {
        let scenarios = [
            (
                Point::new(0.0, 0.0, -5.0),
                Vector::new(0.0, 0.0, 1.0),
                vec![5.0],
            ),
            (
                Point::new(0.0, 0.0, -5.0),
                Vector::new(1.0, 1.0, 1.0),
                vec![8.66025],
            ),
            (
                Point::new(1.0, 1.0, -5.0),
                Vector::new(-0.5, -1.0, 1.0),
                vec![4.55006, 49.44994],
            ),
        ];

        for (origin, direction, ts) in scenarios {
            let shape = Primitive::new_cone(Conic::default());

            let r = Ray::new(origin, direction.normalize());
            let xs = shape
                .intersect(&r, &im_rc::Vector::new())
                .0
                .iter()
                .map(|i| i.t)
                .collect::<Vec<_>>();

            assert_eq!(xs.len(), ts.len());
            for (x_t, t) in xs.iter().zip(ts.iter()) {
                assert_abs_diff_eq!(x_t, t, epsilon = Point::default_epsilon());
            }
        }
    }

    #[test]
    fn intersect_cone_with_ray_parallel_to_one_half() {
        let shape = Primitive::new_cone(Conic::default());

        let direction = Vector::new(0.0, 1.0, 1.0).normalize();
        let r = Ray::new(Point::new(0.0, 0.0, -1.0), direction);

        let xs = shape
            .intersect(&r, &im_rc::Vector::new())
            .0
            .iter()
            .map(|i| i.t)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 0.35355, epsilon = Point::default_epsilon());
    }

    #[test]
    fn intersect_cone_end_caps() {
        let scenarios = [
            (Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 1.0, 0.0), 0),
            (Point::new(0.0, 0.0, -0.25), Vector::new(0.0, 1.0, 1.0), 2),
            (Point::new(0.0, 0.0, -0.25), Vector::new(0.0, 1.0, 0.0), 4),
        ];

        for (origin, direction, count) in scenarios {
            let shape = Primitive::new_cone(Conic::new(-0.5, 0.5, true));

            let r = Ray::new(origin, direction.normalize());
            let xs = shape.intersect(&r, &im_rc::Vector::new());

            assert_eq!(xs.0.len(), count);
        }
    }

    #[test]
    fn intersect_ray_parallel_to_triangle() {
        let t = Primitive::new_triangle(
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Point::new(0.0, -1.0, -2.0), Vector::new(0.0, 1.0, 0.0));

        let xs = t.intersect(&r, &im_rc::Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_misses_p1_p3_edge() {
        let t = Primitive::new_triangle(
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Point::new(1.0, 1.0, -2.0), Vector::new(0.0, 0.0, 1.0));

        let xs = t.intersect(&r, &im_rc::Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_misses_p1_p2_edge() {
        let t = Primitive::new_triangle(
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Point::new(-1.0, 1.0, -2.0), Vector::new(0.0, 0.0, 1.0));

        let xs = t.intersect(&r, &im_rc::Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_misses_p2_p3_edge() {
        let t = Primitive::new_triangle(
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Point::new(0.0, -1.0, -2.0), Vector::new(0.0, 0.0, 1.0));

        let xs = t.intersect(&r, &im_rc::Vector::new());

        assert!(xs.0.is_empty());
    }

    #[test]
    fn ray_strikes_triangle() {
        let t = Primitive::new_triangle(
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        );
        let r = Ray::new(Point::new(0.0, 0.5, -2.0), Vector::new(0.0, 0.0, 1.0));

        let xs = t
            .intersect(&r, &im_rc::Vector::new())
            .0
            .iter()
            .map(|i| i.t)
            .collect::<Vec<_>>();

        assert_eq!(xs.len(), 1);
        assert_abs_diff_eq!(xs[0], 2.0);
    }

    #[test]
    fn intersection_with_smooth_triangle_stores_uv() {
        let p1 = Point::new(0.0, 1.0, 0.0);
        let p2 = Point::new(-1.0, 0.0, 0.0);
        let p3 = Point::new(1.0, 0.0, 0.0);
        let n1 = Vector::new(0.0, 1.0, 0.0);
        let n2 = Vector::new(-1.0, 0.0, 0.0);
        let n3 = Vector::new(1.0, 0.0, 0.0);

        let tri = Primitive::new_smooth_triangle(p1, p2, p3, n1, n2, n3);
        let r = Ray::new(Point::new(-0.2, 0.3, -2.0), Vector::new(0.0, 0.0, 1.0));
        let intersections = tri.intersect(&r, &im_rc::Vector::new());
        let xs = intersections.0.iter().collect::<Vec<_>>();

        assert_eq!(xs.len(), 1);
        let (u, v) = xs[0].u_v.expect("u and v to be set");
        assert_abs_diff_eq!(u, 0.45);
        assert_abs_diff_eq!(v, 0.25);
    }

    #[test]
    fn normal_on_translated_sphere() {
        let shape = Primitive::new_sphere().with_transform(Matrix::translation(0.0, 1.0, 0.0));
        let n = shape.normal_at(
            &Point::new(0.0, 1.70711, -FRAC_1_SQRT_2),
            &Intersection::new(0.0, &shape, im_rc::Vector::new()),
        );

        assert_abs_diff_eq!(n, Vector::new(0.0, FRAC_1_SQRT_2, -FRAC_1_SQRT_2));
    }

    #[test]
    fn normal_on_transformed_sphere() {
        let shape = Primitive::new_sphere()
            .with_transform(Matrix::scaling(1.0, 0.5, 1.0) * Matrix::rotation_z(FRAC_PI_4));
        let n = shape.normal_at(
            &Point::new(0.0, 2.0_f64.sqrt() / 2.0, -(2.0_f64.sqrt()) / 2.0),
            &Intersection::new(0.0, &shape, im_rc::Vector::new()),
        );

        assert_abs_diff_eq!(n, Vector::new(0.0, 0.97014, -0.24254));
    }

    #[test]
    fn sphere_normal_x_axis() {
        let n = Kind::sphere_normal_at(&Point::new(1.0, 0.0, 0.0));

        assert_abs_diff_eq!(n, Vector::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn sphere_normal_y_axis() {
        let n = Kind::sphere_normal_at(&Point::new(0.0, 1.0, 0.0));

        assert_abs_diff_eq!(n, Vector::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn sphere_normal_z_axis() {
        let n = Kind::sphere_normal_at(&Point::new(0.0, 0.0, 1.0));

        assert_abs_diff_eq!(n, Vector::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn sphere_normal_non_axial_point() {
        let n = Kind::sphere_normal_at(&Point::new(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
        ));

        assert_abs_diff_eq!(
            n,
            Vector::new(
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0
            )
        );
    }

    #[test]
    fn normal_is_normalized() {
        let n = Kind::sphere_normal_at(&Point::new(
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
            3.0_f64.sqrt() / 3.0,
        ));

        assert_abs_diff_eq!(n, n.normalize());
    }

    #[test]
    fn normal_of_plane_is_constant() {
        let plane = Primitive::new_plane();
        let hit = Intersection::new(0.0, &plane, im_rc::Vector::new());
        let n1 = plane.normal_at(&Point::new(0.0, 0.0, 0.0), &hit);
        let n2 = plane.normal_at(&Point::new(10.0, 0.0, -10.0), &hit);
        let n3 = plane.normal_at(&Point::new(-5.0, 0.0, 150.0), &hit);

        assert_abs_diff_eq!(n1, Vector::new(0.0, 1.0, 0.0));
        assert_abs_diff_eq!(n2, Vector::new(0.0, 1.0, 0.0));
        assert_abs_diff_eq!(n3, Vector::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn normal_on_surface_of_cube() {
        let scenarios = [
            (Point::new(1.0, 0.5, -0.8), Vector::new(1.0, 0.0, 0.0)),
            (Point::new(-1.0, -0.2, 0.9), Vector::new(-1.0, 0.0, 0.0)),
            (Point::new(-0.4, 1.0, -0.1), Vector::new(0.0, 1.0, 0.0)),
            (Point::new(0.3, -1.0, -0.7), Vector::new(0.0, -1.0, 0.0)),
            (Point::new(-0.6, 0.3, 1.0), Vector::new(0.0, 0.0, 1.0)),
            (Point::new(0.4, 0.4, -1.0), Vector::new(0.0, 0.0, -1.0)),
            (Point::new(1.0, 1.0, 1.0), Vector::new(1.0, 0.0, 0.0)),
            (Point::new(-1.0, -1.0, -1.0), Vector::new(-1.0, 0.0, 0.0)),
        ];

        for (point, normal) in scenarios {
            let n = Kind::cube_normal_at(&point);

            assert_abs_diff_eq!(n, normal);
        }
    }

    #[test]
    fn normal_on_cylinder() {
        let scenarios = [
            (Point::new(1.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
            (Point::new(0.0, 5.0, -1.0), Vector::new(0.0, 0.0, -1.0)),
            (Point::new(0.0, -2.0, 1.0), Vector::new(0.0, 0.0, 1.0)),
            (Point::new(-1.0, 1.0, 0.0), Vector::new(-1.0, 0.0, 0.0)),
        ];

        for (point, normal) in scenarios {
            let n = Kind::cylinder_normal_at(&point, &Conic::default());

            assert_abs_diff_eq!(n, normal);
        }
    }

    #[test]
    fn normal_on_cylinder_end_caps() {
        let scenarios = [
            (Point::new(0.0, 1.0, 0.03), Vector::new(0.0, -1.0, 0.0)),
            (Point::new(0.5, 1.0, 0.0), Vector::new(0.0, -1.0, 0.0)),
            (Point::new(0.0, 1.0, 0.5), Vector::new(0.0, -1.0, 0.0)),
            (Point::new(0.0, 2.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
            (Point::new(0.5, 2.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
            (Point::new(0.0, 2.0, 0.5), Vector::new(0.0, 1.0, 0.0)),
        ];

        for (point, normal) in scenarios {
            let n = Kind::cylinder_normal_at(&point, &Conic::new(1.0, 2.0, true));

            assert_abs_diff_eq!(n, normal);
        }
    }

    #[test]
    fn normal_on_cone() {
        let scenarios = [
            (Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 0.0)),
            (
                Point::new(1.0, 1.0, 1.0),
                Vector::new(1.0, -(2.0_f64.sqrt()), 1.0),
            ),
            (Point::new(-1.0, -1.0, 0.0), Vector::new(-1.0, 1.0, 0.0)),
        ];

        for (point, normal) in scenarios {
            let n = Kind::cone_normal_at(&point, &Conic::default());

            assert_abs_diff_eq!(n, normal);
        }
    }

    #[test]
    fn normal_on_child_object() {
        let s = Primitive::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));
        let g2 = Compound::new_group(vec![s.clone().to_shape()])
            .with_transform(Matrix::scaling(1.0, 2.0, 3.0));
        let g1 = Compound::new_group(vec![g2.clone().to_shape()])
            .with_transform(Matrix::rotation_y(FRAC_PI_2));

        let n = s.normal_at(
            &Point::new(1.7321, 1.1547, -5.5774),
            &Intersection::new(
                0.0,
                &s,
                im_rc::vector![
                    g2.properties().inverse_transform,
                    g1.properties().inverse_transform
                ],
            ),
        );

        assert_abs_diff_eq!(n, Vector::new(0.2857, 0.4286, -0.8571));
    }

    #[test]
    fn normal_on_triangle() {
        let t = Primitive::new_triangle(
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        );
        let hit = Intersection::new(0.0, &t, im_rc::Vector::new());

        assert_eq!(
            t.normal_at(&Point::new(0.0, 0.5, 0.0), &hit),
            Vector::new(0.0, 0.0, 1.0)
        );
        assert_eq!(
            t.normal_at(&Point::new(-0.5, 0.75, 0.0), &hit),
            Vector::new(0.0, 0.0, 1.0)
        );
        assert_eq!(
            t.normal_at(&Point::new(0.5, 0.25, 0.0), &hit),
            Vector::new(0.0, 0.0, 1.0)
        );
    }

    #[test]
    fn smooth_triangle_uses_uv_to_interpolate_normal() {
        let p1 = Point::new(0.0, 1.0, 0.0);
        let p2 = Point::new(-1.0, 0.0, 0.0);
        let p3 = Point::new(1.0, 0.0, 0.0);
        let n1 = Vector::new(0.0, 1.0, 0.0);
        let n2 = Vector::new(-1.0, 0.0, 0.0);
        let n3 = Vector::new(1.0, 0.0, 0.0);

        let tri = Primitive::new_smooth_triangle(p1, p2, p3, n1, n2, n3);
        let i = Intersection::new_with_uv(1.0, &tri, Some((0.45, 0.25)), im_rc::Vector::new());
        let n = tri.normal_at(&Point::new(0.0, 0.0, 0.0), &i);

        assert_eq!(n, Vector::new(-0.5547, 0.83205, 0.0));
    }

    #[test]
    fn convert_point_world_to_object() {
        let s = Primitive::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));
        let g2 = Compound::new_group(vec![s.clone().to_shape()])
            .with_transform(Matrix::scaling(2.0, 2.0, 2.0));
        let g1 = Compound::new_group(vec![g2.clone().to_shape()])
            .with_transform(Matrix::rotation_y(FRAC_PI_2));

        assert_eq!(
            s.world_to_object(
                &Point::new(-2.0, 0.0, -10.0),
                &im_rc::vector![
                    g2.properties().inverse_transform,
                    g1.properties().inverse_transform
                ]
            ),
            Point::new(0.0, 0.0, -1.0)
        );
    }

    #[test]
    fn convert_normal_object_to_world() {
        let s = Primitive::new_sphere().with_transform(Matrix::translation(5.0, 0.0, 0.0));
        let g2 = Compound::new_group(vec![s.clone().to_shape()])
            .with_transform(Matrix::scaling(1.0, 2.0, 3.0));
        let g1 = Compound::new_group(vec![g2.clone().to_shape()])
            .with_transform(Matrix::rotation_y(FRAC_PI_2));

        let n = s.normal_to_world(
            &Vector::new(
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0,
                3.0_f64.sqrt() / 3.0,
            ),
            &im_rc::vector![
                g2.properties().inverse_transform,
                g1.properties().inverse_transform
            ],
        );

        assert_eq!(n, Vector::new(0.2857, 0.4286, -0.8571));
    }

    #[test]
    fn lighting_with_eye_between_light_and_surface() {
        let position = Point::new(0.0, 0.0, 0.0);

        let eye_v = Vector::new(0.0, 0.0, -1.0);
        let normal_v = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 0.0, -10.0), Color::WHITE);
        let shape = Primitive::new_sphere();

        let result = shape.lighting(
            &position,
            &light,
            &eye_v,
            &normal_v,
            false,
            &im_rc::Vector::new(),
        );

        assert_abs_diff_eq!(result, Color::new(1.9, 1.9, 1.9));
    }

    #[test]
    fn lighting_with_eye_between_light_and_surface_eye_offset_45_degrees() {
        let position = Point::new(0.0, 0.0, 0.0);

        let eye_v = Vector::new(0.0, 2_f64.sqrt() / 2.0, -(2_f64.sqrt()) / 2.0);
        let normal_v = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 0.0, -10.0), Color::WHITE);
        let shape = Primitive::new_sphere();

        let result = shape.lighting(
            &position,
            &light,
            &eye_v,
            &normal_v,
            false,
            &im_rc::Vector::new(),
        );

        assert_abs_diff_eq!(result, Color::WHITE);
    }

    #[test]
    fn lighting_with_eye_opposite_surface_light_offset_45_degrees() {
        let position = Point::new(0.0, 0.0, 0.0);

        let eye_v = Vector::new(0.0, 0.0, -1.0);
        let normal_v = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 10.0, -10.0), Color::WHITE);
        let shape = Primitive::new_sphere();

        let result = shape.lighting(
            &position,
            &light,
            &eye_v,
            &normal_v,
            false,
            &im_rc::Vector::new(),
        );

        assert_abs_diff_eq!(result, Color::new(0.7364, 0.7364, 0.7364));
    }

    #[test]
    fn lighting_with_eye_in_path_of_reflection_vector() {
        let position = Point::new(0.0, 0.0, 0.0);

        let eye_v = Vector::new(0.0, -(2_f64.sqrt()) / 2.0, -(2_f64.sqrt()) / 2.0);
        let normal_v = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 10.0, -10.0), Color::WHITE);
        let shape = Primitive::new_sphere();

        let result = shape.lighting(
            &position,
            &light,
            &eye_v,
            &normal_v,
            false,
            &im_rc::Vector::new(),
        );

        assert_abs_diff_eq!(result, Color::new(1.6364, 1.6364, 1.6364));
    }

    #[test]
    fn lighting_with_light_behind_surface() {
        let position = Point::new(0.0, 0.0, 0.0);

        let eye_v = Vector::new(0.0, 0.0, -1.0);
        let normal_v = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 0.0, 10.0), Color::WHITE);
        let shape = Primitive::new_sphere();

        let result = shape.lighting(
            &position,
            &light,
            &eye_v,
            &normal_v,
            false,
            &im_rc::Vector::new(),
        );

        assert_abs_diff_eq!(result, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn lighting_with_surface_in_shadow() {
        let position = Point::new(0.0, 0.0, 0.0);

        let eye_v = Vector::new(0.0, 0.0, -1.0);
        let normal_v = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 0.0, -10.0), Color::WHITE);
        let shape = Primitive::new_sphere();

        assert_abs_diff_eq!(
            shape.lighting(
                &position,
                &light,
                &eye_v,
                &normal_v,
                true,
                &im_rc::Vector::new()
            ),
            Color::new(0.1, 0.1, 0.1)
        );
    }

    #[test]
    fn lighting_with_pattern_applied() {
        let stripe = Pattern::new_stripe(
            Pattern::new_solid(Color::WHITE),
            Pattern::new_solid(Color::BLACK),
        );
        let material = Material {
            ambient: 1.0,
            diffuse: 0.0,
            specular: 0.0,
            pattern: Some(stripe),
            ..Material::default()
        };

        let eye_v = Vector::new(0.0, 0.0, -1.0);
        let normal_v = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 0.0, -10.0), Color::WHITE);
        let shape = Primitive::new_sphere().with_material(material);

        let position = Point::new(0.9, 0.0, 0.0);
        let c1 = shape.lighting(
            &position,
            &light,
            &eye_v,
            &normal_v,
            false,
            &im_rc::Vector::new(),
        );

        let position = Point::new(1.1, 0.0, 0.0);
        let c2 = shape.lighting(
            &position,
            &light,
            &eye_v,
            &normal_v,
            false,
            &im_rc::Vector::new(),
        );

        assert_abs_diff_eq!(c1, Color::WHITE);
        assert_abs_diff_eq!(c2, Color::BLACK);
    }
}
