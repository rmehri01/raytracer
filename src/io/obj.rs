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
    let mut groups = HashMap::new();
    let mut current_group = None;

    for line in obj_string.lines() {
        match &line.split(' ').collect::<Vec<&str>>()[..] {
            ["v", v1, v2, v3] => {
                if let Ok(p) = parse_vertex(v1, v2, v3) {
                    vertices.push(p);
                }
            }
            ["f", vs @ ..] if vs.len() >= 3 => {
                if let Ok(mut t) = fan_triangulation(vs, &vertices) {
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

fn fan_triangulation(vs: &[&str], vertices: &[Tuple]) -> Result<Vec<Shape>, ParseIntError> {
    (2..vs.len())
        .map(|i| parse_triangle(vs[0], vs[i - 1], vs[i], vertices))
        .collect()
}

fn parse_triangle(
    v1: &str,
    v2: &str,
    v3: &str,
    vertices: &[Tuple],
) -> Result<Shape, ParseIntError> {
    let v1 = v1.parse::<usize>()? - 1;
    let v2 = v2.parse::<usize>()? - 1;
    let v3 = v3.parse::<usize>()? - 1;

    Ok(Shape::new_triangle(
        vertices[v1],
        vertices[v2],
        vertices[v3],
    ))
}

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
}
