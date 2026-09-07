pub(super) fn parse_color(value: &str) -> Result<u32, String> {
    let hex = value.strip_prefix('#').unwrap_or(value);
    if hex.len() != 6 || !hex.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        return Err(format!(
            "color must be a six-digit hexadecimal value; got {value}"
        ));
    }
    u32::from_str_radix(hex, 16).map_err(|_| format!("invalid color: {value}"))
}

pub(super) fn color_f32(color: u32) -> [f32; 3] {
    [
        ((color >> 16) & 0xff) as f32 / 255.0,
        ((color >> 8) & 0xff) as f32 / 255.0,
        (color & 0xff) as f32 / 255.0,
    ]
}

pub(super) fn quantize(value: f32) -> u8 {
    // RGBA8 conversion rounds the normalized fragment value, without an
    // intermediate float32 multiplication. For example, f32(0.7) is below
    // 178.5 / 255 even though multiplying it by 255 in f32 produces 178.5.
    (f64::from(value.clamp(0.0, 1.0)) * 255.0).round() as u8
}

#[cfg(test)]
mod tests {
    use super::quantize;

    #[test]
    fn unorm8_conversion_does_not_round_the_scaled_float_twice() {
        assert_eq!(quantize(0.3), 77);
        assert_eq!(quantize(0.7), 178);
        assert_eq!(quantize(0.9), 229);
        assert_eq!(quantize(f32::from_bits(0.7f32.to_bits() + 1)), 179);
        assert_eq!(quantize(f32::from_bits(0.9f32.to_bits() + 1)), 230);
    }

    #[test]
    fn unorm8_conversion_preserves_bytes_and_clamps_endpoints() {
        for byte in 0..=255u8 {
            assert_eq!(quantize(f32::from(byte) / 255.0), byte);
        }
        assert_eq!(quantize(-1.0), 0);
        assert_eq!(quantize(0.5), 128);
        assert_eq!(quantize(2.0), 255);
    }
}
