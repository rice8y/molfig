use std::io::Cursor;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RenderOutputFormat {
    Svg,
    Png,
}

impl RenderOutputFormat {
    pub(crate) fn parse(bytes: &[u8]) -> Result<Self, String> {
        let value =
            std::str::from_utf8(bytes).map_err(|_| "output-format must be UTF-8".to_string())?;
        match value {
            "svg" => Ok(Self::Svg),
            "png" => Ok(Self::Png),
            other => Err(format!(
                "output-format must be one of \"svg\" or \"png\"; got {other}"
            )),
        }
    }
}

pub(crate) fn encode_render_output(
    format: RenderOutputFormat,
    width: u32,
    height: u32,
    rgba: &[u8],
) -> Result<Vec<u8>, String> {
    let png = encode_png(width, height, rgba)?;
    match format {
        RenderOutputFormat::Png => Ok(png),
        RenderOutputFormat::Svg => Ok(encode_svg(width, height, &png)),
    }
}

fn encode_png(width: u32, height: u32, rgba: &[u8]) -> Result<Vec<u8>, String> {
    let expected = (width as usize)
        .checked_mul(height as usize)
        .and_then(|pixels| pixels.checked_mul(4))
        .ok_or_else(|| "render target dimensions overflow".to_string())?;
    if rgba.len() != expected {
        return Err(format!(
            "render target byte length mismatch: expected {expected}, got {}",
            rgba.len()
        ));
    }

    let mut encoded = Vec::new();
    {
        let mut encoder = png::Encoder::new(Cursor::new(&mut encoded), width, height);
        encoder.set_color(png::ColorType::Rgba);
        encoder.set_depth(png::BitDepth::Eight);
        let mut writer = encoder
            .write_header()
            .map_err(|error| format!("failed to encode PNG header: {error}"))?;
        writer
            .write_image_data(rgba)
            .map_err(|error| format!("failed to encode PNG pixels: {error}"))?;
    }
    Ok(encoded)
}

fn encode_svg(width: u32, height: u32, png: &[u8]) -> Vec<u8> {
    let encoded = base64(png);
    format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" viewBox=\"0 0 {width} {height}\"><image width=\"{width}\" height=\"{height}\" preserveAspectRatio=\"none\" href=\"data:image/png;base64,{encoded}\"/></svg>"
    )
    .into_bytes()
}

fn base64(bytes: &[u8]) -> String {
    const TABLE: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut encoded = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        let a = u32::from(chunk[0]);
        let b = u32::from(*chunk.get(1).unwrap_or(&0));
        let c = u32::from(*chunk.get(2).unwrap_or(&0));
        let value = (a << 16) | (b << 8) | c;
        encoded.push(TABLE[((value >> 18) & 0x3f) as usize] as char);
        encoded.push(TABLE[((value >> 12) & 0x3f) as usize] as char);
        encoded.push(if chunk.len() > 1 {
            TABLE[((value >> 6) & 0x3f) as usize] as char
        } else {
            '='
        });
        encoded.push(if chunk.len() > 2 {
            TABLE[(value & 0x3f) as usize] as char
        } else {
            '='
        });
    }
    encoded
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_format_is_strictly_svg_or_png() {
        assert_eq!(
            RenderOutputFormat::parse(b"svg"),
            Ok(RenderOutputFormat::Svg)
        );
        assert_eq!(
            RenderOutputFormat::parse(b"png"),
            Ok(RenderOutputFormat::Png)
        );
        assert!(RenderOutputFormat::parse(b"obj").is_err());
        assert!(RenderOutputFormat::parse(b"SVG").is_err());
    }

    #[test]
    fn png_encodes_rgba8_without_changing_pixels() {
        let rgba = [255, 0, 0, 255, 0, 128, 255, 64];
        let png = encode_render_output(RenderOutputFormat::Png, 2, 1, &rgba).unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");

        let decoder = png::Decoder::new(Cursor::new(png));
        let mut reader = decoder.read_info().unwrap();
        let mut decoded = vec![0; reader.output_buffer_size().unwrap()];
        let info = reader.next_frame(&mut decoded).unwrap();
        assert_eq!(info.width, 2);
        assert_eq!(info.height, 1);
        assert_eq!(info.color_type, png::ColorType::Rgba);
        assert_eq!(&decoded[..info.buffer_size()], &rgba);
    }

    #[test]
    fn svg_is_a_dimensioned_lossless_png_container() {
        let rgba = [1, 2, 3, 4];
        let svg = encode_render_output(RenderOutputFormat::Svg, 1, 1, &rgba).unwrap();
        let text = std::str::from_utf8(&svg).unwrap();
        assert!(
            text.starts_with("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1\" height=\"1\"")
        );
        assert!(text.contains("href=\"data:image/png;base64,iVBORw0KGgo"));
        assert!(text.ends_with("/></svg>"));
    }

    #[test]
    fn base64_matches_rfc_4648_vectors() {
        assert_eq!(base64(b""), "");
        assert_eq!(base64(b"f"), "Zg==");
        assert_eq!(base64(b"fo"), "Zm8=");
        assert_eq!(base64(b"foo"), "Zm9v");
        assert_eq!(base64(b"foobar"), "Zm9vYmFy");
    }
}
