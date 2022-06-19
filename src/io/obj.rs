use core::fmt;
use std::{collections::HashMap, error::Error, fs, io};

use crate::{
    core::{point::Point, vector::Vector},
    raytracer::shapes::{Compound, Primitive, Shape},
};

#[derive(Debug)]
pub enum ParseError {
    Logic(String),
    Syntax(String),
    Io(io::Error),
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::Logic(s) => write!(f, "Logic error, {}", s),
            Self::Syntax(s) => write!(f, "Syntax error, {}", s),
            Self::Io(e) => write!(f, "IO error, {}", e),
        }
    }
}

type Result<T> = std::result::Result<T, ParseError>;

impl Error for ParseError {}

#[derive(Debug, Default)]
struct Parser<'a> {
    ignored_lines: usize,
    shapes: Vec<Shape>,
    vertices: Vec<Point>,
    normals: Vec<Vector>,
    groups: HashMap<&'a str, Vec<Primitive>>,
}

impl Parser<'_> {
    fn make_shape(mut self) -> Result<Shape> {
        for (_, group) in self.groups {
            self.shapes.push(
                Compound::new_group(group.into_iter().map(|t| t.to_shape()).collect()).to_shape(),
            );
        }

        match self.shapes.len() {
            0 => Err(ParseError::Logic("no shapes found".to_string())),
            1 => Ok(self.shapes.pop().unwrap()),
            _ => Ok(Compound::new_group(self.shapes).to_shape()),
        }
    }
}

pub fn parse_file(path: &str) -> Result<Shape> {
    let obj_string = fs::read_to_string(path).map_err(ParseError::Io)?;
    let shape = parse_shape(&obj_string)?;

    Ok(shape)
}

fn parse_shape(obj_string: &str) -> Result<Shape> {
    let parser = parse_string(obj_string)?;
    let shape = parser.make_shape()?;

    Ok(shape)
}

/// Parses a string representation of an OBJ file, deliberately ignoring
/// any lines that aren't recognized.
fn parse_string(obj_string: &str) -> Result<Parser> {
    let mut parser = Parser::default();
    let mut current_group = None;

    for line in obj_string.lines() {
        match &line.trim().split(' ').collect::<Vec<_>>()[..] {
            ["v", x, y, z] => {
                let p = parse_vertex(x, y, z)?;
                parser.vertices.push(p);
            }
            ["vn", x, y, z] => {
                let v = parse_normal(x, y, z)?;
                parser.normals.push(v);
            }
            ["f", vs @ ..] if vs.len() >= 3 => {
                let mut ts = fan_triangulation(vs, &parser.vertices, &parser.normals)?;

                match current_group {
                    Some(group) => {
                        parser
                            .groups
                            .entry(group)
                            .or_insert_with(Vec::new)
                            .append(&mut ts);
                    }
                    None => {
                        let mut triangle_shapes = ts.into_iter().map(|t| t.to_shape()).collect();
                        parser.shapes.append(&mut triangle_shapes);
                    }
                };
            }
            ["g", group_name] => {
                current_group = Some(*group_name);
            }
            _ => {
                parser.ignored_lines += 1;
            }
        }
    }

    Ok(parser)
}

fn parse_vertex(x: &str, y: &str, z: &str) -> Result<Point> {
    let err_fn = |_| ParseError::Syntax(format!("invalid vertex: {} {} {}", x, y, z));
    let x = x.parse().map_err(err_fn)?;
    let y = y.parse().map_err(err_fn)?;
    let z = z.parse().map_err(err_fn)?;

    Ok(Point::new(x, y, z))
}

fn parse_normal(x: &str, y: &str, z: &str) -> Result<Vector> {
    let err_fn = |_| ParseError::Syntax(format!("invalid normal: {} {} {}", x, y, z));
    let x = x.parse().map_err(err_fn)?;
    let y = y.parse().map_err(err_fn)?;
    let z = z.parse().map_err(err_fn)?;

    Ok(Vector::new(x, y, z))
}

fn fan_triangulation(
    vs: &[&str],
    vertices: &[Point],
    normals: &[Vector],
) -> Result<Vec<Primitive>> {
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
) -> Result<Primitive> {
    let (p1, n1) = parse_vertex_ref(v1, vertices, normals)?;
    let (p2, n2) = parse_vertex_ref(v2, vertices, normals)?;
    let (p3, n3) = parse_vertex_ref(v3, vertices, normals)?;

    match (n1, n2, n3) {
        (Some(n1), Some(n2), Some(n3)) => {
            Ok(Primitive::new_smooth_triangle(p1, p2, p3, n1, n2, n3))
        }
        (None, None, None) => Ok(Primitive::new_triangle(p1, p2, p3)),
        _ => Err(ParseError::Syntax(format!(
            "invalid triangle: {} {} {}",
            v1, v2, v3
        ))),
    }
}

fn parse_vertex_ref(
    v: &str,
    vertices: &[Point],
    normals: &[Vector],
) -> Result<(Point, Option<Vector>)> {
    let err_fn = |_| ParseError::Syntax(format!("indices must be natural numbers, given: {}", v));

    match v.split('/').collect::<Vec<_>>()[..] {
        [v, _, n] => {
            let v = v.parse::<usize>().map_err(err_fn)?;
            let n = n.parse::<usize>().map_err(err_fn)?;

            Ok((vertices[v - 1], Some(normals[n - 1])))
        }
        [v] | [v, _] => {
            let v = v.parse::<usize>().map_err(err_fn)?;

            Ok((vertices[v - 1], None))
        }
        _ => Err(ParseError::Syntax(format!("invalid vertex ref: {}", v))),
    }
}

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

        assert_eq!(parse_string(&obj_string).unwrap().ignored_lines, 5);
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

        let shape = parse_shape(&obj_string).unwrap();

        assert_eq!(
            shape,
            Compound::new_group(vec![
                Primitive::new_triangle(
                    Point::new(-1.0, 1.0, 0.0),
                    Point::new(-1.0, 0.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                )
                .to_shape(),
                Primitive::new_triangle(
                    Point::new(-1.0, 1.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                    Point::new(1.0, 1.0, 0.0),
                )
                .to_shape()
            ])
            .to_shape()
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

        let shape = parse_shape(&obj_string).unwrap();

        assert_eq!(
            shape,
            Compound::new_group(vec![
                Primitive::new_triangle(
                    Point::new(-1.0, 1.0, 0.0),
                    Point::new(-1.0, 0.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                )
                .to_shape(),
                Primitive::new_triangle(
                    Point::new(-1.0, 1.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                    Point::new(1.0, 1.0, 0.0),
                )
                .to_shape(),
                Primitive::new_triangle(
                    Point::new(-1.0, 1.0, 0.0),
                    Point::new(1.0, 1.0, 0.0),
                    Point::new(0.0, 2.0, 0.0),
                )
                .to_shape()
            ])
            .to_shape()
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

        let shape = parse_shape(&obj_string).unwrap();

        let t1 = Primitive::new_triangle(
            Point::new(-1.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        )
        .to_shape();
        let g1 = Compound::new_group(vec![t1]).to_shape();

        let t2 = Primitive::new_triangle(
            Point::new(-1.0, 1.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
        )
        .to_shape();
        let g2 = Compound::new_group(vec![t2]).to_shape();

        assert!(
            shape == Compound::new_group(vec![g1.clone(), g2.clone()]).to_shape()
                || shape == Compound::new_group(vec![g2, g1]).to_shape()
        );
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

        let shape = parse_shape(&obj_string).unwrap();

        let t1 = Primitive::new_smooth_triangle(
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Vector::new(0.0, 1.0, 0.0),
            Vector::new(-1.0, 0.0, 0.0),
            Vector::new(1.0, 0.0, 0.0),
        )
        .to_shape();
        assert_eq!(shape, Compound::new_group(vec![t1.clone(), t1]).to_shape());
    }
}
