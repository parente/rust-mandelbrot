use crossbeam;
use image::png::PNGEncoder;
use image::ColorType;
use log::{debug, info};
use num::Complex;
use std::fs::File;
use std::str::FromStr;
use structopt::StructOpt;

/// Compute whether the given complex number is in the mandelbrot set after at most limit
/// iterations.
fn escape_time(c: Complex<f64>, limit: usize) -> Option<usize> {
    let mut z = Complex { re: 0.0, im: 0.0 };
    for i in 0..limit {
        if z.norm_sqr() > 4.0 {
            return Some(i);
        }
        z = z * z + c
    }
    None
}

/// Parse the string `s` as a coordinate pair, like `"400x600"` or `"1.0,0.5"`.
///
/// Specifically, `s` should have the form <left><sep><right>, where <sep> is
/// the character given by the `separator` argument, and <left> and <right> are
/// both strings that can be parsed by `T::from_str`. `separator` must be an
/// ASCII character.
///
/// If `s` has the proper form, return `Some<(x, y)>`. If it doesn't parse
/// correctly, return `None`.
fn parse_pair<T: FromStr>(s: &str, separator: char) -> Option<(T, T)> {
    match s.find(separator) {
        None => None,
        Some(index) => match (T::from_str(&s[..index]), T::from_str(&s[index + 1..])) {
            (Ok(l), Ok(r)) => Some((l, r)),
            _ => None,
        },
    }
}

#[test]
fn test_parse_pair() {
    assert_eq!(parse_pair::<i32>("", ','), None);
    assert_eq!(parse_pair::<i32>("10,", ','), None);
    assert_eq!(parse_pair::<i32>(",10", ','), None);
    assert_eq!(parse_pair::<i32>("10,20", ','), Some((10, 20)));
    assert_eq!(parse_pair::<i32>("10,20xy", ','), None);
    assert_eq!(parse_pair::<f64>("0.5", 'x'), None);
    assert_eq!(parse_pair::<f64>("0.5x1.5", 'x'), Some((0.5, 1.5)));
}

/// Parse a pair of integers separated by an x as pixel dimensions.
fn parse_dimensions(s: &str) -> Result<(usize, usize), &'static str> {
    match parse_pair(s, 'x') {
        Some((w, h)) => Ok((w, h)),
        None => Err("Could not parse dimensions"),
    }
}

#[test]
fn test_parse_dimensions() {
    assert_eq!(parse_dimensions("100x200"), Ok((100, 200)));
    assert_eq!(
        parse_dimensions(",-0.0625"),
        Err("Could not parse dimensions")
    );
}

/// Parse a pair of floating-point numbers separated by a comma as a complex number.
fn parse_complex(s: &str) -> Result<Complex<f64>, &'static str> {
    match parse_pair(s, ',') {
        Some((re, im)) => Ok(Complex { re, im }),
        None => Err("Could not parse complex value"),
    }
}

#[test]
fn test_parse_complex() {
    assert_eq!(
        parse_complex("1.25,-0.0625"),
        Ok(Complex {
            re: 1.25,
            im: -0.0625
        })
    );
    assert_eq!(
        parse_complex(",-0.0625"),
        Err("Could not parse complex value")
    );
}

/// Given the row and column of a pixel in the output image, return the corresponding point on the
/// complex plane.assert_eq!
///
/// `bounds` is a pair givin the width and height of the image in pixels.
/// `pixel` is a (column, row) pair indicating a particular pixel.
/// The `upper_left` and `lower_right` parameters are points on the complex plane designating the
/// area our image covers.
fn pixel_to_point(
    bounds: (usize, usize),
    pixel: (usize, usize),
    upper_left: Complex<f64>,
    lower_right: Complex<f64>,
) -> Complex<f64> {
    let (width, height) = (
        lower_right.re - upper_left.re,
        upper_left.im - lower_right.im,
    );
    Complex {
        re: upper_left.re + pixel.0 as f64 * width / bounds.0 as f64,
        im: upper_left.im - pixel.1 as f64 * height / bounds.1 as f64,
    }
}

#[test]
fn test_pixel_to_point() {
    assert_eq!(
        pixel_to_point(
            (100, 200),
            (25, 175),
            Complex { re: -1.0, im: 1.0 },
            Complex { re: 1.0, im: -1.0 }
        ),
        Complex {
            re: -0.5,
            im: -0.75
        }
    );
}

/// Render a rectangle of the Mandelbrot set into a buffer of pixels.
///
/// The `bounds` argument gives the width and height of the buffer `pixels`,
/// which holds one grayscale pixel per byte. The `upper_left` and `lower_right`
/// arguments specify points on the complex plane corresponding to the upper-
/// left and lower-right corners of the pixel buffer.
fn render(
    pixels: &mut [u8],
    bounds: (usize, usize),
    upper_left: Complex<f64>,
    lower_right: Complex<f64>,
) {
    assert!(pixels.len() == bounds.0 * bounds.1);
    for row in 0..bounds.1 {
        for column in 0..bounds.0 {
            let point = pixel_to_point(bounds, (column, row), upper_left, lower_right);
            pixels[row * bounds.0 + column] = match escape_time(point, 255) {
                None => 0,
                Some(count) => 255 - count as u8,
            }
        }
    }
}

/// Write the buffer `pixels`, whose dimensions are given by `bounds`, to the
/// file named `filename`.
fn write_image(
    filename: &str,
    pixels: &[u8],
    bounds: (usize, usize),
) -> Result<(), std::io::Error> {
    let output = File::create(filename)?;
    let encoder = PNGEncoder::new(output);
    encoder.encode(pixels, bounds.0 as u32, bounds.1 as u32, ColorType::Gray(8))?;
    Ok(())
}

#[derive(StructOpt)]
struct Opt {
    /// Output PNG filename
    filename: String,
    /// Pixel dimensions of the image, <width>x<height>
    #[structopt(long, parse(try_from_str = parse_dimensions))]
    bounds: (usize, usize),
    /// Upper-left bound of the complex space, <real>,<imaginary>
    #[structopt(long, parse(try_from_str = parse_complex))]
    upper_left: Complex<f64>,
    /// Lower-right bound of the complex space, <real>,<imaginary>
    #[structopt(long, parse(try_from_str = parse_complex))]
    lower_right: Complex<f64>,
    /// Log in structured JSON format
    #[structopt(long)]
    json: bool,
}

fn main() {
    let opt = Opt::from_args();

    if opt.json {
        json_env_logger::init();
    } else {
        env_logger::init();
    }

    let mut pixels = vec![0; opt.bounds.0 * opt.bounds.1];

    let threads = 16;
    let rows_per_band = opt.bounds.1 / threads + 1;
    debug!("Assigning {} rows per thread", rows_per_band);

    info!("Starting render");
    {
        let bands: Vec<&mut [u8]> = pixels.chunks_mut(rows_per_band * opt.bounds.0).collect();
        crossbeam::scope(|spawner| {
            for (i, band) in bands.into_iter().enumerate() {
                let top = rows_per_band * i;
                let height = band.len() / opt.bounds.0;
                let band_bounds = (opt.bounds.0, height);
                let band_upper_left =
                    pixel_to_point(opt.bounds, (0, top), opt.upper_left, opt.lower_right);
                let band_lower_right = pixel_to_point(
                    opt.bounds,
                    (opt.bounds.0, top + height),
                    opt.upper_left,
                    opt.lower_right,
                );
                spawner.spawn(move |_| {
                    render(band, band_bounds, band_upper_left, band_lower_right);
                    debug!("Finished band {} render", i);
                });
            }
        })
        .unwrap();
    }
    write_image(&opt.filename[..], &pixels, opt.bounds).expect("error writing PNG file");
    info!("Wrote file: {}", opt.filename);
}
