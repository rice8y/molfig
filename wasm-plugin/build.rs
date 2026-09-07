use std::env;
use std::fs;
use std::io::Cursor;
use std::path::{Path, PathBuf};

fn main() {
    println!("cargo:rerun-if-changed=src/render/postprocessing/smaa-data/area.png.b64");
    println!("cargo:rerun-if-changed=src/render/postprocessing/smaa-data/search.png.b64");

    let output = PathBuf::from(env::var_os("OUT_DIR").expect("OUT_DIR is set by Cargo"));
    decode_lookup(
        "src/render/postprocessing/smaa-data/area.png.b64",
        &output.join("smaa-area.rgb"),
        160,
        560,
        png::ColorType::Rgb,
    );
    decode_lookup(
        "src/render/postprocessing/smaa-data/search.png.b64",
        &output.join("smaa-search.r"),
        66,
        33,
        png::ColorType::Grayscale,
    );
}

fn decode_lookup(
    source: &str,
    destination: &Path,
    expected_width: u32,
    expected_height: u32,
    expected_color: png::ColorType,
) {
    let encoded = fs::read_to_string(source)
        .unwrap_or_else(|error| panic!("failed to read {source}: {error}"));
    let png_bytes = decode_base64(encoded.trim())
        .unwrap_or_else(|error| panic!("failed to decode {source}: {error}"));
    let decoder = png::Decoder::new(Cursor::new(png_bytes));
    let mut reader = decoder
        .read_info()
        .unwrap_or_else(|error| panic!("failed to parse {source}: {error}"));
    let mut pixels = vec![
        0;
        reader
            .output_buffer_size()
            .expect("bounded SMAA lookup texture")
    ];
    let info = reader
        .next_frame(&mut pixels)
        .unwrap_or_else(|error| panic!("failed to decode {source}: {error}"));
    assert_eq!((info.width, info.height), (expected_width, expected_height));
    assert_eq!(info.bit_depth, png::BitDepth::Eight);
    assert_eq!(info.color_type, expected_color);
    pixels.truncate(info.buffer_size());
    fs::write(destination, pixels)
        .unwrap_or_else(|error| panic!("failed to write {}: {error}", destination.display()));
}

fn decode_base64(input: &str) -> Result<Vec<u8>, &'static str> {
    if !input.len().is_multiple_of(4) {
        return Err("length is not divisible by four");
    }
    let mut output = Vec::with_capacity(input.len() / 4 * 3);
    for chunk in input.as_bytes().chunks_exact(4) {
        let a = base64_value(chunk[0])?;
        let b = base64_value(chunk[1])?;
        let c = if chunk[2] == b'=' {
            0
        } else {
            base64_value(chunk[2])?
        };
        let d = if chunk[3] == b'=' {
            0
        } else {
            base64_value(chunk[3])?
        };
        output.push((a << 2) | (b >> 4));
        if chunk[2] != b'=' {
            output.push((b << 4) | (c >> 2));
        }
        if chunk[3] != b'=' {
            output.push((c << 6) | d);
        }
    }
    Ok(output)
}

fn base64_value(byte: u8) -> Result<u8, &'static str> {
    match byte {
        b'A'..=b'Z' => Ok(byte - b'A'),
        b'a'..=b'z' => Ok(byte - b'a' + 26),
        b'0'..=b'9' => Ok(byte - b'0' + 52),
        b'+' => Ok(62),
        b'/' => Ok(63),
        _ => Err("invalid base64 character"),
    }
}
