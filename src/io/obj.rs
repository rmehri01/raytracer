use std::{
    collections::HashMap,
    fs, io,
    num::{ParseFloatError, ParseIntError},
};

use crate::{core::tuple::Tuple, raytracer::shape::Shape};

pub fn parse_obj_file(path: &str) -> io::Result<Shape> {
    let obj_string = fs::read_to_string(path)?;

    Ok(parse_obj_string(obj_string))
}

/// Parses a string representation of an OBJ file, deliberately ignoring
/// any lines that aren't recognized.
fn parse_obj_string(obj_string: String) -> Shape {
    let mut triangles = Vec::new();
    let mut vertices = Vec::new();
    let mut normals = Vec::new();
    let mut groups = HashMap::new();
    let mut current_group = None;

    for line in obj_string.lines() {
        match &line.trim().split(' ').collect::<Vec<&str>>()[..] {
            ["v", v1, v2, v3] => {
                if let Ok(p) = parse_vertex(v1, v2, v3) {
                    vertices.push(p);
                }
            }
            ["vn", v1, v2, v3] => {
                if let Ok(p) = parse_normal(v1, v2, v3) {
                    normals.push(p);
                }
            }
            ["f", vs @ ..] if vs.len() >= 3 => {
                if let Ok(mut t) = fan_triangulation(vs, &vertices, &normals) {
                    match current_group {
                        Some(group) => {
                            groups.entry(group).or_insert_with(Vec::new).append(&mut t);
                        }
                        None => triangles.append(&mut t),
                    }
                }
            }
            ["g", group_name] => {
                current_group = Some(*group_name);
            }
            _ => {}
        }
    }

    for (_, group) in groups {
        triangles.push(Shape::new_group(group));
    }

    Shape::new_group(triangles)
}

fn parse_vertex(v1: &str, v2: &str, v3: &str) -> Result<Tuple, ParseFloatError> {
    let v1 = v1.parse::<f64>()?;
    let v2 = v2.parse::<f64>()?;
    let v3 = v3.parse::<f64>()?;

    Ok(Tuple::point(v1, v2, v3))
}

fn parse_normal(v1: &str, v2: &str, v3: &str) -> Result<Tuple, ParseFloatError> {
    let v1 = v1.parse::<f64>()?;
    let v2 = v2.parse::<f64>()?;
    let v3 = v3.parse::<f64>()?;

    Ok(Tuple::vector(v1, v2, v3))
}

fn fan_triangulation(
    vs: &[&str],
    vertices: &[Tuple],
    normals: &[Tuple],
) -> Result<Vec<Shape>, ParseIntError> {
    (2..vs.len())
        .map(|i| parse_triangle(vs[0], vs[i - 1], vs[i], vertices, normals))
        .collect()
}

fn parse_triangle(
    v1: &str,
    v2: &str,
    v3: &str,
    vertices: &[Tuple],
    normals: &[Tuple],
) -> Result<Shape, ParseIntError> {
    let (v1, n1) = parse_vertex_ref(v1, vertices, normals)?;
    let (v2, n2) = parse_vertex_ref(v2, vertices, normals)?;
    let (v3, n3) = parse_vertex_ref(v3, vertices, normals)?;

    match (n1, n2, n3) {
        (Some(n1), Some(n2), Some(n3)) => Ok(Shape::new_smooth_triangle(v1, v2, v3, n1, n2, n3)),
        (None, None, None) => Ok(Shape::new_triangle(v1, v2, v3)),
        // TODO: real error type
        _ => todo!("other error"),
    }
}

fn parse_vertex_ref(
    v: &str,
    vertices: &[Tuple],
    normals: &[Tuple],
) -> Result<(Tuple, Option<Tuple>), ParseIntError> {
    match v.split('/').collect::<Vec<&str>>()[..] {
        [v, _, n] => {
            let v = v.parse::<usize>()? - 1;
            let n = n.parse::<usize>()? - 1;

            Ok((vertices[v], Some(normals[n])))
        }
        [v] | [v, _] => {
            let v = v.parse::<usize>()? - 1;

            Ok((vertices[v], None))
        }
        _ => todo!("other error"),
    }
}

// TODO: helper or newtype for 1-based indexing

#[cfg(test)]
mod tests {
    use crate::raytracer::shape::ShapeKind;

    use super::*;

    #[test]
    fn ignore_unrecognized_lines() {
        let obj_string = [
            "There was a young lady named Bright",
            "who traveled much faster than light.",
            "She set out one day",
            "in a relative way,",
            "and came back the previous night.",
        ]
        .join("\n");

        parse_obj_string(obj_string);
    }

    #[test]
    fn vertex_records() {
        let v1 = "1.0";
        let v2 = "2.0";
        let v3 = "3.0";

        let result = parse_vertex(v1, v2, v3);

        assert_eq!(result.unwrap(), Tuple::point(1.0, 2.0, 3.0));
    }

    #[test]
    fn parsing_triangle_faces() {
        let obj_string = [
            "v -1 1 0", "v -1 0 0", "v 1 0 0", "v 1 1 0", "f 1 2 3", "f 1 3 4",
        ]
        .join("\n");

        let shape = parse_obj_string(obj_string);

        if let ShapeKind::Group(group) = shape.kind {
            assert_eq!(group.len(), 2);
            assert_eq!(
                group[0],
                Shape::new_triangle(
                    Tuple::point(-1.0, 1.0, 0.0),
                    Tuple::point(-1.0, 0.0, 0.0),
                    Tuple::point(1.0, 0.0, 0.0),
                )
            );
            assert_eq!(
                group[1],
                Shape::new_triangle(
                    Tuple::point(-1.0, 1.0, 0.0),
                    Tuple::point(1.0, 0.0, 0.0),
                    Tuple::point(1.0, 1.0, 0.0),
                )
            );
        } else {
            panic!("expected a group of two triangles");
        }
    }

    #[test]
    fn triangulating_polygons() {
        let obj_string = [
            "v -1 1 0",
            "v -1 0 0",
            "v 1 0 0",
            "v 1 1 0",
            "v 0 2 0",
            "f 1 2 3 4 5",
        ]
        .join("\n");

        let shape = parse_obj_string(obj_string);

        if let ShapeKind::Group(group) = shape.kind {
            assert_eq!(group.len(), 3);
            assert_eq!(
                group[0],
                Shape::new_triangle(
                    Tuple::point(-1.0, 1.0, 0.0),
                    Tuple::point(-1.0, 0.0, 0.0),
                    Tuple::point(1.0, 0.0, 0.0),
                )
            );
            assert_eq!(
                group[1],
                Shape::new_triangle(
                    Tuple::point(-1.0, 1.0, 0.0),
                    Tuple::point(1.0, 0.0, 0.0),
                    Tuple::point(1.0, 1.0, 0.0),
                )
            );
            assert_eq!(
                group[2],
                Shape::new_triangle(
                    Tuple::point(-1.0, 1.0, 0.0),
                    Tuple::point(1.0, 1.0, 0.0),
                    Tuple::point(0.0, 2.0, 0.0),
                )
            );
        } else {
            panic!("expected a group of three triangles");
        }
    }

    #[test]
    fn triangles_in_named_groups() {
        let obj_string = [
            "v -1 1 0",
            "v -1 0 0",
            "v 1 0 0",
            "v 1 1 0",
            "g FirstGroup",
            "f 1 2 3",
            "g SecondGroup",
            "f 1 3 4",
        ]
        .join("\n");

        let shape = parse_obj_string(obj_string);

        if let ShapeKind::Group(group) = shape.kind {
            assert_eq!(group.len(), 2);
            let t1 = Shape::new_triangle(
                Tuple::point(-1.0, 1.0, 0.0),
                Tuple::point(-1.0, 0.0, 0.0),
                Tuple::point(1.0, 0.0, 0.0),
            );
            let g1 = Shape::new_group(vec![t1]);
            assert!(group.contains(&g1));

            let t2 = Shape::new_triangle(
                Tuple::point(-1.0, 1.0, 0.0),
                Tuple::point(1.0, 0.0, 0.0),
                Tuple::point(1.0, 1.0, 0.0),
            );
            let g2 = Shape::new_group(vec![t2]);
            assert!(group.contains(&g2));
        } else {
            panic!("expected a group of two triangles");
        }
    }

    #[test]
    fn faces_with_normals() {
        let obj_string = [
            "v 0 1 0",
            "v -1 0 0",
            "v 1 0 0",
            "vn -1 0 0",
            "vn 1 0 0",
            "vn 0 1 0",
            "f 1//3 2//1 3//2",
            "f 1/0/3 2/102/1 3/14/2",
        ]
        .join("\n");

        let shape = parse_obj_string(obj_string);
        // ShapeKind::SmoothTriangle {triangular, n1, n2, n3}
        if let ShapeKind::Group(group) = shape.kind {
            let t1 = &group[0];
            if let ShapeKind::SmoothTriangle {
                triangular,
                n1,
                n2,
                n3,
            } = t1.kind
            {
                assert_eq!(triangular.p1, Tuple::point(0.0, 1.0, 0.0));
                assert_eq!(triangular.p2, Tuple::point(-1.0, 0.0, 0.0));
                assert_eq!(triangular.p3, Tuple::point(1.0, 0.0, 0.0));
                assert_eq!(n1, Tuple::vector(0.0, 1.0, 0.0));
                assert_eq!(n2, Tuple::vector(-1.0, 0.0, 0.0));
                assert_eq!(n3, Tuple::vector(1.0, 0.0, 0.0));
            } else {
                panic!("expected a smooth triangle");
            }

            let t2 = &group[1];
            assert_eq!(t1, t2);
        } else {
            panic!("expected a group of two triangles");
        }
    }
}
