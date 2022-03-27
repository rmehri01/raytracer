/// A tuple using left-hand coordinates where:
/// * `w = 0.0` is a vector
/// * `w = 1.0` is a point
#[derive(Debug, PartialEq)]
pub struct Tuple {
    x: f64,
    y: f64,
    z: f64,
    w: f64,
}

impl Tuple {
    pub fn new_point(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z, w: 1.0 }
    }

    pub fn new_vector(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z, w: 0.0 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tuple_point() {
        let p = Tuple {
            x: 4.3,
            y: -4.2,
            z: 3.1,
            w: 1.0,
        };

        assert_eq!(p.x, 4.3);
        assert_eq!(p.y, -4.2);
        assert_eq!(p.z, 3.1);
        assert_eq!(p.w, 1.0);
    }

    #[test]
    fn test_tuple_vector() {
        let v = Tuple {
            x: 4.3,
            y: -4.2,
            z: 3.1,
            w: 0.0,
        };

        assert_eq!(v.x, 4.3);
        assert_eq!(v.y, -4.2);
        assert_eq!(v.z, 3.1);
        assert_eq!(v.w, 0.0);
    }

    #[test]
    fn test_tuple_point_equals_point() {
        let p1 = Tuple::new_point(4.3, -4.2, 3.1);
        let p2 = Tuple {
            x: 4.3,
            y: -4.2,
            z: 3.1,
            w: 1.0,
        };

        assert_eq!(p1, p2);
    }

    #[test]
    fn test_tuple_vector_equals_vector() {
        let v1 = Tuple::new_vector(4.3, -4.2, 3.1);
        let v2 = Tuple {
            x: 4.3,
            y: -4.2,
            z: 3.1,
            w: 0.0,
        };

        assert_eq!(v1, v2);
    }
}
