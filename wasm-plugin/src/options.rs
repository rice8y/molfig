#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InputFormat {
    Auto,
    Pdb,
    Cif,
    BinaryCif,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Representation {
    Molstar,
    Spacefill,
    BallAndStick,
    Cartoon,
    Ribbon,
    Backbone,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PolymerProfile {
    Elliptical,
    Rounded,
    Square,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum VisualQuality {
    Auto,
    Custom,
    Highest,
    Higher,
    High,
    Medium,
    Low,
    Lower,
    Lowest,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ColorTheme {
    ChainId,
}

#[derive(Clone, Debug)]
pub struct MeshOptions {
    pub format: InputFormat,
    pub(crate) representation: Representation,
    pub(crate) color_theme: ColorTheme,
    pub sphere_detail: usize,
    pub radius_scale: f32,
    pub atom_radius: f32,
    pub bond_radius: f32,
    pub infer_bonds: bool,
    pub center: bool,
    pub assembly: Option<String>,
    pub alt_loc: String,
    pub(crate) block_index: Option<usize>,
    pub(crate) block_header: Option<String>,
    pub ribbon_radius: f32,
    pub ribbon_width: f32,
    pub(crate) helix_profile: PolymerProfile,
    pub(crate) round_cap: bool,
    pub(crate) sheet_arrow_factor: f32,
    pub(crate) tubular_helices: bool,
    pub(crate) linear_segments: usize,
    pub(crate) radial_segments: usize,
    pub(crate) quality: Option<VisualQuality>,
    pub(crate) visuals: Vec<String>,
    pub(crate) obj_basename: Option<String>,
    pub(crate) include_operator_metadata: bool,
    pub(crate) obj_groups: bool,
}

impl Default for MeshOptions {
    fn default() -> Self {
        Self {
            format: InputFormat::Auto,
            representation: Representation::Molstar,
            color_theme: ColorTheme::ChainId,
            sphere_detail: 2,
            radius_scale: 1.0,
            atom_radius: 0.28,
            bond_radius: 0.12,
            infer_bonds: true,
            center: true,
            assembly: Some("1".to_string()),
            alt_loc: String::new(),
            block_index: None,
            block_header: None,
            ribbon_radius: 0.2,
            ribbon_width: 0.55,
            helix_profile: PolymerProfile::Elliptical,
            round_cap: false,
            sheet_arrow_factor: 1.5,
            tubular_helices: false,
            linear_segments: 8,
            radial_segments: 16,
            quality: None,
            visuals: Vec::new(),
            obj_basename: None,
            include_operator_metadata: true,
            obj_groups: true,
        }
    }
}

impl MeshOptions {
    pub(crate) fn from_json(bytes: &[u8]) -> Result<Self, String> {
        let text =
            std::str::from_utf8(bytes).map_err(|_| "options must be UTF-8 JSON".to_string())?;
        let mut options = MeshOptions::default();

        if let Some(value) = json_string(text, "format") {
            options.format = match value.to_ascii_lowercase().as_str() {
                "auto" => InputFormat::Auto,
                "pdb" => InputFormat::Pdb,
                "cif" | "mmcif" => InputFormat::Cif,
                "bcif" | "binarycif" | "binary-cif" => InputFormat::BinaryCif,
                other => {
                    return Err(format!(
                        "unsupported format: {other}; expected one of \"auto\", \"pdb\", \"cif\", \"mmcif\", \"bcif\", \"binarycif\", or \"binary-cif\""
                    ))
                }
            };
        }

        if let Some(value) = json_string(text, "representation") {
            options.representation = match value.to_ascii_lowercase().as_str() {
                "spacefill" | "space-fill" => Representation::Spacefill,
                "ball-and-stick" | "ball_and_stick" | "balls" => Representation::BallAndStick,
                "cartoon" => Representation::Cartoon,
                "ribbon" => Representation::Ribbon,
                "backbone" => Representation::Backbone,
                "molstar" | "default" | "auto" => Representation::Molstar,
                other => return Err(format!("unsupported representation: {other}")),
            };
        }

        if let Some(value) = json_string(text, "color-theme") {
            options.color_theme = parse_color_theme(&value)?;
        }
        if let Some(value) = json_string(text, "color_theme") {
            options.color_theme = parse_color_theme(&value)?;
        }

        if let Some(value) = json_number(text, "sphere-detail") {
            options.sphere_detail = value.clamp(1.0, 5.0) as usize;
        }
        if let Some(value) = json_number(text, "sphere_detail") {
            options.sphere_detail = value.clamp(1.0, 5.0) as usize;
        }
        if let Some(value) = json_number(text, "radius-scale") {
            options.radius_scale = value.clamp(0.05, 5.0);
        }
        if let Some(value) = json_number(text, "radius_scale") {
            options.radius_scale = value.clamp(0.05, 5.0);
        }
        if let Some(value) = json_number(text, "atom-radius") {
            options.atom_radius = value.clamp(0.02, 5.0);
        }
        if let Some(value) = json_number(text, "atom_radius") {
            options.atom_radius = value.clamp(0.02, 5.0);
        }
        if let Some(value) = json_number(text, "bond-radius") {
            options.bond_radius = value.clamp(0.01, 2.0);
        }
        if let Some(value) = json_number(text, "bond_radius") {
            options.bond_radius = value.clamp(0.01, 2.0);
        }
        if let Some(value) = json_bool(text, "infer-bonds") {
            options.infer_bonds = value;
        }
        if let Some(value) = json_bool(text, "infer_bonds") {
            options.infer_bonds = value;
        }
        if let Some(value) = json_bool(text, "center") {
            options.center = value;
        }
        if let Some(value) = json_string(text, "assembly") {
            options.assembly = (!value.is_empty() && value != "none" && value != "asymmetric-unit")
                .then_some(value);
        }
        if let Some(value) = json_string(text, "alt-loc") {
            options.alt_loc = value;
        }
        if let Some(value) = json_string(text, "alt_loc") {
            options.alt_loc = value;
        }
        if let Some(value) = json_number(text, "block-index") {
            options.block_index = Some(value.max(0.0) as usize);
        }
        if let Some(value) = json_number(text, "block_index") {
            options.block_index = Some(value.max(0.0) as usize);
        }
        if let Some(value) = json_string(text, "block-header") {
            options.block_header = (!value.is_empty()).then_some(value);
        }
        if let Some(value) = json_string(text, "block_header") {
            options.block_header = (!value.is_empty()).then_some(value);
        }
        if let Some(value) = json_number(text, "ribbon-radius") {
            options.ribbon_radius = value.clamp(0.03, 2.0);
        }
        if let Some(value) = json_number(text, "ribbon_radius") {
            options.ribbon_radius = value.clamp(0.03, 2.0);
        }
        if let Some(value) = json_number(text, "ribbon-width") {
            options.ribbon_width = value.clamp(0.05, 4.0);
        }
        if let Some(value) = json_number(text, "ribbon_width") {
            options.ribbon_width = value.clamp(0.05, 4.0);
        }
        if let Some(value) = json_string(text, "helix-profile") {
            options.helix_profile = parse_profile(&value)?;
        }
        if let Some(value) = json_string(text, "helix_profile") {
            options.helix_profile = parse_profile(&value)?;
        }
        if let Some(value) = json_bool(text, "round-cap") {
            options.round_cap = value;
        }
        if let Some(value) = json_bool(text, "round_cap") {
            options.round_cap = value;
        }
        if let Some(value) = json_number(text, "sheet-arrow-factor") {
            options.sheet_arrow_factor = value.clamp(0.0, 3.0);
        }
        if let Some(value) = json_number(text, "sheet_arrow_factor") {
            options.sheet_arrow_factor = value.clamp(0.0, 3.0);
        }
        if let Some(value) = json_bool(text, "tubular-helices") {
            options.tubular_helices = value;
        }
        if let Some(value) = json_bool(text, "tubular_helices") {
            options.tubular_helices = value;
        }
        if let Some(value) = json_number(text, "linear-segments") {
            options.linear_segments = value.clamp(1.0, 48.0) as usize;
        }
        if let Some(value) = json_number(text, "linear_segments") {
            options.linear_segments = value.clamp(1.0, 48.0) as usize;
        }
        if let Some(value) = json_number(text, "radial-segments") {
            options.radial_segments = clamp_radial_segments(value);
        }
        if let Some(value) = json_number(text, "radial_segments") {
            options.radial_segments = clamp_radial_segments(value);
        }
        if let Some(value) = json_string(text, "quality") {
            options.quality = Some(parse_quality(&value)?);
        }
        if let Some(values) = json_string_array(text, "visuals") {
            options.visuals = values
                .into_iter()
                .map(|visual| visual.trim().to_ascii_lowercase())
                .filter(|visual| !visual.is_empty())
                .collect();
        }
        if let Some(value) = json_string(text, "obj-basename") {
            options.obj_basename = normalize_obj_basename(value);
        }
        if let Some(value) = json_string(text, "obj_basename") {
            options.obj_basename = normalize_obj_basename(value);
        }
        if let Some(value) = json_bool(text, "operator-metadata") {
            options.include_operator_metadata = value;
        }
        if let Some(value) = json_bool(text, "operator_metadata") {
            options.include_operator_metadata = value;
        }
        if let Some(value) = json_bool(text, "obj-groups") {
            options.obj_groups = value;
        }
        if let Some(value) = json_bool(text, "obj_groups") {
            options.obj_groups = value;
        }
        if options.quality.is_none() {
            options.quality = Some(VisualQuality::Auto);
        }

        Ok(options)
    }

    pub(crate) fn resolved_for_quality(&self, element_count: usize, is_coarse: bool) -> Self {
        let Some(quality) = self.quality else {
            return self.clone();
        };
        let quality = match quality {
            VisualQuality::Auto => auto_quality(element_count, is_coarse),
            VisualQuality::Custom => return self.clone(),
            quality => quality,
        };
        let mut options = self.clone();
        apply_quality(&mut options, quality);
        options
    }
}

fn clamp_radial_segments(value: f32) -> usize {
    let value = value.clamp(2.0, 56.0) as usize;
    if value == 2 {
        2
    } else {
        value.max(3)
    }
}

fn parse_profile(value: &str) -> Result<PolymerProfile, String> {
    match value.to_ascii_lowercase().as_str() {
        "elliptical" | "ellipse" => Ok(PolymerProfile::Elliptical),
        "rounded" | "round" => Ok(PolymerProfile::Rounded),
        "square" | "sheet" => Ok(PolymerProfile::Square),
        other => Err(format!(
            "unsupported helix-profile: {other}; expected one of \"elliptical\", \"rounded\", or \"square\""
        )),
    }
}

fn parse_quality(value: &str) -> Result<VisualQuality, String> {
    match value.to_ascii_lowercase().as_str() {
        "auto" => Ok(VisualQuality::Auto),
        "custom" => Ok(VisualQuality::Custom),
        "highest" => Ok(VisualQuality::Highest),
        "higher" => Ok(VisualQuality::Higher),
        "high" => Ok(VisualQuality::High),
        "medium" => Ok(VisualQuality::Medium),
        "low" => Ok(VisualQuality::Low),
        "lower" => Ok(VisualQuality::Lower),
        "lowest" => Ok(VisualQuality::Lowest),
        other => Err(format!(
            "unsupported quality: {other}; expected one of \"auto\", \"custom\", \"highest\", \"higher\", \"high\", \"medium\", \"low\", \"lower\", or \"lowest\""
        )),
    }
}

fn parse_color_theme(value: &str) -> Result<ColorTheme, String> {
    match value.to_ascii_lowercase().as_str() {
        "chain-id" | "chain_id" | "chain" | "molstar" | "default" | "auto" => {
            Ok(ColorTheme::ChainId)
        }
        other => Err(format!(
            "unsupported color-theme: {other}; expected \"chain-id\""
        )),
    }
}

fn normalize_obj_basename(value: String) -> Option<String> {
    let mut value = value.trim().to_string();
    if value.to_ascii_lowercase().ends_with(".mtl") {
        value.truncate(value.len().saturating_sub(4));
    }
    (!value.is_empty()).then_some(value)
}

fn auto_quality(element_count: usize, is_coarse: bool) -> VisualQuality {
    let score = if is_coarse {
        element_count.saturating_mul(10)
    } else {
        element_count
    };
    if score > 1_000_000 {
        VisualQuality::Lowest
    } else if score > 500_000 {
        VisualQuality::Lower
    } else if score > 100_000 {
        VisualQuality::Low
    } else if score > 20_000 {
        VisualQuality::Medium
    } else if score > 2_000 {
        VisualQuality::High
    } else {
        VisualQuality::Higher
    }
}

fn apply_quality(options: &mut MeshOptions, quality: VisualQuality) {
    match quality {
        VisualQuality::Highest => {
            options.sphere_detail = 3;
            options.radial_segments = 36;
            options.linear_segments = 18;
        }
        VisualQuality::Higher => {
            options.sphere_detail = 3;
            options.radial_segments = 28;
            options.linear_segments = 14;
        }
        VisualQuality::High => {
            options.sphere_detail = 2;
            options.radial_segments = 20;
            options.linear_segments = 10;
        }
        VisualQuality::Medium => {
            options.sphere_detail = 1;
            options.radial_segments = 12;
            options.linear_segments = 8;
        }
        VisualQuality::Low => {
            options.sphere_detail = 0;
            options.radial_segments = 8;
            options.linear_segments = 3;
        }
        VisualQuality::Lower => {
            options.sphere_detail = 0;
            options.radial_segments = 4;
            options.linear_segments = 2;
        }
        VisualQuality::Lowest => {
            options.sphere_detail = 0;
            options.radial_segments = 2;
            options.linear_segments = 1;
        }
        VisualQuality::Auto | VisualQuality::Custom => {}
    }
}

fn json_string(text: &str, key: &str) -> Option<String> {
    let marker = format!("\"{key}\"");
    let mut tail = text.split_once(&marker)?.1;
    tail = tail.split_once(':')?.1.trim_start();
    if !tail.starts_with('"') {
        return None;
    }
    let mut value = String::new();
    let mut escape = false;
    for ch in tail[1..].chars() {
        if escape {
            value.push(ch);
            escape = false;
        } else if ch == '\\' {
            escape = true;
        } else if ch == '"' {
            return Some(value);
        } else {
            value.push(ch);
        }
    }
    None
}

fn json_number(text: &str, key: &str) -> Option<f32> {
    let marker = format!("\"{key}\"");
    let mut tail = text.split_once(&marker)?.1;
    tail = tail.split_once(':')?.1.trim_start();
    let end = tail
        .find(|c: char| !(c.is_ascii_digit() || matches!(c, '-' | '+' | '.' | 'e' | 'E')))
        .unwrap_or(tail.len());
    tail[..end].parse().ok()
}

fn json_bool(text: &str, key: &str) -> Option<bool> {
    let marker = format!("\"{key}\"");
    let mut tail = text.split_once(&marker)?.1;
    tail = tail.split_once(':')?.1.trim_start();
    if tail.starts_with("true") {
        Some(true)
    } else if tail.starts_with("false") {
        Some(false)
    } else {
        None
    }
}

fn json_string_array(text: &str, key: &str) -> Option<Vec<String>> {
    let marker = format!("\"{key}\"");
    let mut tail = text.split_once(&marker)?.1;
    tail = tail.split_once(':')?.1.trim_start();
    if !tail.starts_with('[') {
        return None;
    }
    let mut values = Vec::new();
    let mut in_string = false;
    let mut escape = false;
    let mut current = String::new();
    for ch in tail[1..].chars() {
        if in_string {
            if escape {
                current.push(ch);
                escape = false;
            } else if ch == '\\' {
                escape = true;
            } else if ch == '"' {
                values.push(std::mem::take(&mut current));
                in_string = false;
            } else {
                current.push(ch);
            }
            continue;
        }
        match ch {
            '"' => in_string = true,
            ']' => return Some(values),
            _ => {}
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::{MeshOptions, VisualQuality};

    #[test]
    fn named_quality_matches_molstar_quality_props() {
        let options = MeshOptions::from_json(br#"{"quality":"higher","linear-segments":1}"#)
            .unwrap()
            .resolved_for_quality(100, false);

        assert_eq!(options.quality, Some(VisualQuality::Higher));
        assert_eq!(options.sphere_detail, 3);
        assert_eq!(options.radial_segments, 28);
        assert_eq!(options.linear_segments, 14);
    }

    #[test]
    fn auto_quality_uses_molstar_structure_thresholds() {
        let options = MeshOptions::from_json(br#"{"quality":"auto"}"#)
            .unwrap()
            .resolved_for_quality(2_870, false);

        assert_eq!(options.sphere_detail, 2);
        assert_eq!(options.radial_segments, 20);
        assert_eq!(options.linear_segments, 10);
    }

    #[test]
    fn omitted_quality_defaults_to_molstar_auto_for_public_options() {
        let options = MeshOptions::from_json(br#"{"sphere-detail":1,"linear-segments":1}"#)
            .unwrap()
            .resolved_for_quality(2_870, false);

        assert_eq!(options.quality, Some(VisualQuality::Auto));
        assert_eq!(options.sphere_detail, 2);
        assert_eq!(options.radial_segments, 20);
        assert_eq!(options.linear_segments, 10);
    }

    #[test]
    fn custom_quality_keeps_explicit_custom_segments() {
        let options = MeshOptions::from_json(
            br#"{"quality":"custom","sphere-detail":1,"linear-segments":6,"radial-segments":18}"#,
        )
        .unwrap()
        .resolved_for_quality(10, false);

        assert_eq!(options.sphere_detail, 1);
        assert_eq!(options.radial_segments, 18);
        assert_eq!(options.linear_segments, 6);
    }
}
