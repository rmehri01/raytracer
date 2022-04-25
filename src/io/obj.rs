use std::{
    collections::HashMap,
    fs, io,
    num::{ParseFloatError, ParseIntError},
};

use crate::{
    core::{point::Point, vector::Vector},
    raytracer::shapes::{Compound, Shape, Single},
};

pub fn parse_file(path: &str) -> io::Result<Shape> {
    let obj_string = fs::read_to_string(path)?;

    Ok(parse_string(&obj_string))
}

/// Parses a string representation of an OBJ file, deliberately ignoring
/// any lines that aren't recognized.
fn parse_string(obj_string: &str) -> Shape {
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
                        None => {
                            triangles.append(&mut t.into_iter().map(|t| t.as_shape()).collect());
                        }
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
        triangles.push(
            Compound::new_group(group.into_iter().map(|t| t.as_shape()).collect()).as_shape(),
        );
    }

    Compound::new_group(triangles).as_shape()
}

fn parse_vertex(v1: &str, v2: &str, v3: &str) -> Result<Point, ParseFloatError> {
    let v1 = v1.parse::<f64>()?;
    let v2 = v2.parse::<f64>()?;
    let v3 = v3.parse::<f64>()?;

    Ok(Point::new(v1, v2, v3))
}

fn parse_normal(v1: &str, v2: &str, v3: &str) -> Result<Vector, ParseFloatError> {
    let v1 = v1.parse::<f64>()?;
    let v2 = v2.parse::<f64>()?;
    let v3 = v3.parse::<f64>()?;

    Ok(Vector::new(v1, v2, v3))
}

fn fan_triangulation(
    vs: &[&str],
    vertices: &[Point],
    normals: &[Vector],
) -> Result<Vec<Single>, ParseIntError> {
    (2..vs.len())
        .map(|i| parse_triangle(vs[0], vs[i - 1], vs[i], vertices, normals))
        .collect()
}

fn parse_triangle(
    v1: &str,
    v2: &str,
    v3: &str,
    vertices: &[Point],
    normals: &[Vector],
) -> Result<Single, ParseIntError> {
    let (v1, n1) = parse_vertex_ref(v1, vertices, normals)?;
    let (v2, n2) = parse_vertex_ref(v2, vertices, normals)?;
    let (v3, n3) = parse_vertex_ref(v3, vertices, normals)?;

    match (n1, n2, n3) {
        (Some(n1), Some(n2), Some(n3)) => Ok(Single::new_smooth_triangle(v1, v2, v3, n1, n2, n3)),
        (None, None, None) => Ok(Single::new_triangle(v1, v2, v3)),
        // TODO: real error type
        _ => todo!("other error"),
    }
}

fn parse_vertex_ref(
    v: &str,
    vertices: &[Point],
    normals: &[Vector],
) -> Result<(Point, Option<Vector>), ParseIntError> {
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

        parse_string(&obj_string);
    }

    #[test]
    fn vertex_records() {
        let v1 = "1.0";
        let v2 = "2.0";
        let v3 = "3.0";

        let result = parse_vertex(v1, v2, v3);

        assert_eq!(result.unwrap(), Point::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn parsing_triangle_faces() {
        let obj_string = [
            "v -1 1 0", "v -1 0 0", "v 1 0 0", "v 1 1 0", "f 1 2 3", "f 1 3 4",
        ]
        .join("\n");

        let shape = parse_string(&obj_string);

        assert_eq!(
            shape,
            Compound::new_group(vec![
                Single::new_triangle(
                    Point::new(-1.0, 1.0, 0.0),
                    Point::new(-1.0, 0.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                )
                .as_shape(),
                Single::new_triangle(
                    Point::new(-1.0, 1.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                    Point::new(1.0, 1.0, 0.0),
                )
                .as_shape()
            ])
            .as_shape()
        );
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

        let shape = parse_string(&obj_string);

        assert_eq!(
            shape,
            Compound::new_group(vec![
                Single::new_triangle(
                    Point::new(-1.0, 1.0, 0.0),
                    Point::new(-1.0, 0.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                )
                .as_shape(),
                Single::new_triangle(
                    Point::new(-1.0, 1.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                    Point::new(1.0, 1.0, 0.0),
                )
                .as_shape(),
                Single::new_triangle(
                    Point::new(-1.0, 1.0, 0.0),
                    Point::new(1.0, 1.0, 0.0),
                    Point::new(0.0, 2.0, 0.0),
                )
                .as_shape()
            ])
            .as_shape()
        );
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

        let shape = parse_string(&obj_string);

        let t1 = Single::new_triangle(
            Point::new(-1.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        )
        .as_shape();
        let g1 = Compound::new_group(vec![t1]).as_shape();

        let t2 = Single::new_triangle(
            Point::new(-1.0, 1.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
        )
        .as_shape();
        let g2 = Compound::new_group(vec![t2]).as_shape();

        assert_eq!(shape, Compound::new_group(vec![g1, g2]).as_shape());
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

        let shape = parse_string(&obj_string);

        let t1 = Single::new_smooth_triangle(
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Vector::new(0.0, 1.0, 0.0),
            Vector::new(-1.0, 0.0, 0.0),
            Vector::new(1.0, 0.0, 0.0),
        )
        .as_shape();
        assert_eq!(shape, Compound::new_group(vec![t1.clone(), t1]).as_shape());
    }
}
