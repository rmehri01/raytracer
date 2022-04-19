use std::fs;

use itertools::Itertools;

use super::color::Color;

pub struct Canvas {
    pub width: usize,
    pub height: usize,
    pixels: Vec<Color>,
}

impl Canvas {
    const PPM_IDENTIFIER: &'static str = "P3";
    const PPM_MAX_COLOR_VALUE: u8 = 255;
    const PPM_MAX_LINE_LEN: u8 = 70;

    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            pixels: vec![Color::BLACK; width * height],
        }
    }

    pub fn write_pixel(&mut self, x: usize, y: usize, color: Color) {
        let index = y * self.width + x;
        self.pixels[index] = color;
    }

    pub fn pixel_at(&self, x: usize, y: usize) -> Color {
        let index = y * self.width + x;
        self.pixels[index]
    }

    pub fn write_ppm(&self, filename: &str) {
        let ppm = self.to_ppm();

        fs::write(filename, ppm).expect("file should be written successfully");
    }

    fn to_ppm(&self) -> String {
        let mut ppm = String::new();

        ppm.push_str(&format!(
            "{}\n{} {}\n{}\n",
            Self::PPM_IDENTIFIER,
            self.width,
            self.height,
            Self::PPM_MAX_COLOR_VALUE
        ));

        for row in self.pixels.chunks(self.width) {
            let mut colors = row
                .iter()
                .flat_map(|color| {
                    vec![
                        Self::scale_to_ppm_data(color.red),
                        Self::scale_to_ppm_data(color.green),
                        Self::scale_to_ppm_data(color.blue),
                    ]
                })
                .map(|color| color.to_string())
                .peekable();

            while colors.peek().is_some() {
                let mut line_len = 0;
                let line = colors
                    .take_while_ref(|color| {
                        line_len += color.len() + 1;
                        line_len <= Self::PPM_MAX_LINE_LEN as usize
                    })
                    .join(" ");
                ppm.push_str(&line);
                ppm.push('\n');
            }
        }

        ppm
    }

    fn scale_to_ppm_data(color_scale: f64) -> u8 {
        let scaled_data = color_scale * f64::from(Self::PPM_MAX_COLOR_VALUE);

        scaled_data
            .clamp(0.0, f64::from(Self::PPM_MAX_COLOR_VALUE))
            .round() as u8
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn new() {
        let canvas = Canvas::new(10, 20);

        assert_eq!(canvas.width, 10);
        assert_eq!(canvas.height, 20);
        for pixel in canvas.pixels {
            assert_relative_eq!(pixel.red, 0.0);
            assert_relative_eq!(pixel.blue, 0.0);
            assert_relative_eq!(pixel.green, 0.0);
        }
    }

    #[test]
    fn write_pixel() {
        let mut canvas = Canvas::new(10, 20);
        let red = Color::new(1.0, 0.0, 0.0);

        canvas.write_pixel(2, 3, red);

        let pixel = canvas.pixel_at(2, 3);
        assert_relative_eq!(pixel.red, 1.0);
        assert_relative_eq!(pixel.blue, 0.0);
        assert_relative_eq!(pixel.green, 0.0);
    }

    #[test]
    fn ppm_header() {
        let canvas = Canvas::new(5, 3);

        let ppm = canvas.to_ppm();

        let mut lines = ppm.lines();
        assert_eq!(lines.next().unwrap(), "P3");
        assert_eq!(lines.next().unwrap(), "5 3");
        assert_eq!(lines.next().unwrap(), "255");
    }

    #[test]
    fn ppm_pixel_data() {
        let mut canvas = Canvas::new(5, 3);
        let c1 = Color::new(1.5, 0.0, 0.0);
        let c2 = Color::new(0.0, 0.5, 0.0);
        let c3 = Color::new(-0.5, 0.0, 1.0);

        canvas.write_pixel(0, 0, c1);
        canvas.write_pixel(2, 1, c2);
        canvas.write_pixel(4, 2, c3);

        let ppm = canvas.to_ppm();

        let mut lines = ppm.lines();
        lines.next();
        lines.next();
        lines.next();

        assert_eq!(lines.next().unwrap(), "255 0 0 0 0 0 0 0 0 0 0 0 0 0 0");
        assert_eq!(lines.next().unwrap(), "0 0 0 0 0 0 0 128 0 0 0 0 0 0 0");
        assert_eq!(lines.next().unwrap(), "0 0 0 0 0 0 0 0 0 0 0 0 0 0 255");
        assert!(lines.next().is_none());
    }

    #[test]
    fn ppm_split_long_lines() {
        let mut canvas = Canvas::new(10, 2);

        for row in 0..canvas.height {
            for col in 0..canvas.width {
                canvas.write_pixel(col, row, Color::new(1.0, 0.8, 0.6));
            }
        }

        let ppm = canvas.to_ppm();

        let mut lines = ppm.lines();
        lines.next();
        lines.next();
        lines.next();

        assert_eq!(
            lines.next().unwrap(),
            "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204"
        );
        assert_eq!(
            lines.next().unwrap(),
            "153 255 204 153 255 204 153 255 204 153 255 204 153"
        );
        assert_eq!(
            lines.next().unwrap(),
            "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204"
        );
        assert_eq!(
            lines.next().unwrap(),
            "153 255 204 153 255 204 153 255 204 153 255 204 153"
        );
    }

    #[test]
    fn ppm_terminate_with_newline() {
        let canvas = Canvas::new(5, 3);

        let ppm = canvas.to_ppm();

        assert!(ppm.ends_with('\n'));
    }
}
