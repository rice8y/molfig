use std::ops::Range;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DiffReport {
    pub passed: bool,
    pub message: String,
    pub details: Vec<(String, String)>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DiffBytesOptions {
    pub stl_delta_scan: bool,
    pub stl_facet_range: Option<Range<usize>>,
}

impl Default for DiffBytesOptions {
    fn default() -> Self {
        Self {
            stl_delta_scan: true,
            stl_facet_range: None,
        }
    }
}

pub fn diff_bytes(reference: &[u8], generated: &[u8], label: &str) -> DiffReport {
    diff_bytes_with_options(reference, generated, label, DiffBytesOptions::default())
}

pub fn diff_bytes_with_options(
    reference: &[u8],
    generated: &[u8],
    label: &str,
    options: DiffBytesOptions,
) -> DiffReport {
    if reference == generated {
        return DiffReport {
            passed: true,
            message: format!(
                "PASS {label}: byte-for-byte match ({} bytes)",
                reference.len()
            ),
            details: vec![
                ("kind".to_string(), "bytes".to_string()),
                ("reference_len".to_string(), reference.len().to_string()),
                ("generated_len".to_string(), generated.len().to_string()),
            ],
        };
    }

    if let Some(range) = options.stl_facet_range.as_ref() {
        if let Some(report) = diff_stl_facet_range(reference, generated, label, range, &options) {
            return report;
        }
    }

    let first_diff = first_diff_byte(
        reference,
        generated,
        0,
        reference.len().min(generated.len()),
    );
    let reference_byte = reference
        .get(first_diff)
        .map(|byte| format!("0x{byte:02x}"))
        .unwrap_or_else(|| "<eof>".to_string());
    let generated_byte = generated
        .get(first_diff)
        .map(|byte| format!("0x{byte:02x}"))
        .unwrap_or_else(|| "<eof>".to_string());

    let stl_context = stl_diff_context(reference, generated, first_diff, &options);
    let stl_context_message = stl_context
        .as_ref()
        .map(|context| format!("; {context}"))
        .unwrap_or_default();
    let mut details = vec![
        ("kind".to_string(), "bytes".to_string()),
        ("first_byte".to_string(), first_diff.to_string()),
        ("reference_byte".to_string(), reference_byte.clone()),
        ("generated_byte".to_string(), generated_byte.clone()),
        ("reference_len".to_string(), reference.len().to_string()),
        ("generated_len".to_string(), generated.len().to_string()),
    ];
    if let Some(context) = stl_context {
        if let Some((_, delta_scan)) = context.split_once("stl_delta_scan=") {
            details.push(("stl_delta_scan".to_string(), delta_scan.to_string()));
        }
        details.push(("stl_context".to_string(), context));
    }

    DiffReport {
        passed: false,
        message: format!(
            "FAIL {label}: first difference at byte {first_diff}; reference={reference_byte}, generated={generated_byte}; reference_len={}, generated_len={}{stl_context_message}",
            reference.len(),
            generated.len()
        ),
        details,
    }
}

fn first_diff_byte(reference: &[u8], generated: &[u8], start: usize, end: usize) -> usize {
    reference[start..end]
        .iter()
        .zip(generated[start..end].iter())
        .position(|(a, b)| a != b)
        .map(|offset| start + offset)
        .unwrap_or_else(|| {
            if reference.len() == generated.len() {
                end
            } else {
                reference.len().min(generated.len())
            }
        })
}

fn diff_stl_facet_range(
    reference: &[u8],
    generated: &[u8],
    label: &str,
    range: &Range<usize>,
    options: &DiffBytesOptions,
) -> Option<DiffReport> {
    if !looks_like_binary_stl(reference) || !looks_like_binary_stl(generated) {
        return None;
    }
    let reference_count = stl_facet_count(reference)?;
    let generated_count = stl_facet_count(generated)?;
    let common_count = reference_count.min(generated_count);
    if range.start >= range.end || range.end > common_count {
        return Some(DiffReport {
            passed: false,
            message: format!(
                "FAIL {label}: invalid STL facet range {}..{} for reference_facet_count={reference_count}, generated_facet_count={generated_count}",
                range.start, range.end
            ),
            details: vec![
                ("kind".to_string(), "bytes".to_string()),
                ("stl_facet_range".to_string(), format!("{}..{}", range.start, range.end)),
                ("reference_facet_count".to_string(), reference_count.to_string()),
                ("generated_facet_count".to_string(), generated_count.to_string()),
            ],
        });
    }

    let start = 84 + range.start * 50;
    let end = 84 + range.end * 50;
    let first_diff = first_diff_byte(reference, generated, start, end);
    if first_diff == end {
        let mut details = vec![
            ("kind".to_string(), "bytes".to_string()),
            ("reference_len".to_string(), reference.len().to_string()),
            ("generated_len".to_string(), generated.len().to_string()),
            (
                "stl_facet_range".to_string(),
                format!("{}..{}", range.start, range.end),
            ),
            (
                "reference_facet_count".to_string(),
                reference_count.to_string(),
            ),
            (
                "generated_facet_count".to_string(),
                generated_count.to_string(),
            ),
        ];
        if reference_count != generated_count {
            details.push((
                "stl_count_mismatch_outside_range".to_string(),
                "true".to_string(),
            ));
        }
        return Some(DiffReport {
            passed: true,
            message: format!(
                "PASS {label}: STL facet range {}..{} byte-for-byte match ({} facets)",
                range.start,
                range.end,
                range.end - range.start
            ),
            details,
        });
    }

    let reference_byte = reference
        .get(first_diff)
        .map(|byte| format!("0x{byte:02x}"))
        .unwrap_or_else(|| "<eof>".to_string());
    let generated_byte = generated
        .get(first_diff)
        .map(|byte| format!("0x{byte:02x}"))
        .unwrap_or_else(|| "<eof>".to_string());
    let stl_context = stl_diff_context(reference, generated, first_diff, options);
    let stl_context_message = stl_context
        .as_ref()
        .map(|context| format!("; {context}"))
        .unwrap_or_default();
    let mut details = vec![
        ("kind".to_string(), "bytes".to_string()),
        ("first_byte".to_string(), first_diff.to_string()),
        ("reference_byte".to_string(), reference_byte.clone()),
        ("generated_byte".to_string(), generated_byte.clone()),
        ("reference_len".to_string(), reference.len().to_string()),
        ("generated_len".to_string(), generated.len().to_string()),
        (
            "stl_facet_range".to_string(),
            format!("{}..{}", range.start, range.end),
        ),
        (
            "reference_facet_count".to_string(),
            reference_count.to_string(),
        ),
        (
            "generated_facet_count".to_string(),
            generated_count.to_string(),
        ),
    ];
    if let Some(context) = stl_context {
        if let Some((_, delta_scan)) = context.split_once("stl_delta_scan=") {
            details.push(("stl_delta_scan".to_string(), delta_scan.to_string()));
        }
        details.push(("stl_context".to_string(), context));
    }

    Some(DiffReport {
        passed: false,
        message: format!(
            "FAIL {label}: first difference in STL facet range {}..{} at byte {first_diff}; reference={reference_byte}, generated={generated_byte}; reference_len={}, generated_len={}{stl_context_message}",
            range.start,
            range.end,
            reference.len(),
            generated.len()
        ),
        details,
    })
}

pub fn diff_text(reference: &str, generated: &str, label: &str) -> DiffReport {
    if reference == generated {
        return DiffReport {
            passed: true,
            message: format!("PASS {label}: text match ({} bytes)", reference.len()),
            details: vec![
                ("kind".to_string(), "text".to_string()),
                ("reference_len".to_string(), reference.len().to_string()),
                ("generated_len".to_string(), generated.len().to_string()),
            ],
        };
    }

    let reference_lines = reference.lines().collect::<Vec<_>>();
    let generated_lines = generated.lines().collect::<Vec<_>>();
    let first_diff = reference_lines
        .iter()
        .zip(generated_lines.iter())
        .position(|(a, b)| a != b)
        .unwrap_or_else(|| reference_lines.len().min(generated_lines.len()));
    let reference_line = reference_lines.get(first_diff).copied().unwrap_or("<eof>");
    let generated_line = generated_lines.get(first_diff).copied().unwrap_or("<eof>");

    DiffReport {
        passed: false,
        message: format!(
            "FAIL {label}: first difference at line {}; reference=\"{}\", generated=\"{}\"; reference_lines={}, generated_lines={}",
            first_diff + 1,
            truncate(reference_line, 160),
            truncate(generated_line, 160),
            reference_lines.len(),
            generated_lines.len()
        ),
        details: vec![
            ("kind".to_string(), "text".to_string()),
            ("first_line".to_string(), (first_diff + 1).to_string()),
            ("reference_line".to_string(), truncate(reference_line, 160)),
            ("generated_line".to_string(), truncate(generated_line, 160)),
            (
                "reference_lines".to_string(),
                reference_lines.len().to_string(),
            ),
            (
                "generated_lines".to_string(),
                generated_lines.len().to_string(),
            ),
        ],
    }
}

fn truncate(value: &str, limit: usize) -> String {
    let mut out = value.chars().take(limit).collect::<String>();
    if value.chars().count() > limit {
        out.push_str("...");
    }
    out
}

fn stl_diff_context(
    reference: &[u8],
    generated: &[u8],
    first_diff: usize,
    options: &DiffBytesOptions,
) -> Option<String> {
    if !looks_like_binary_stl(reference) || !looks_like_binary_stl(generated) || first_diff < 84 {
        return None;
    }
    let delta_scan = options
        .stl_delta_scan
        .then(|| stl_delta_scan_summary(reference, generated))
        .flatten()
        .map(|summary| format!("; stl_delta_scan={summary}"))
        .unwrap_or_default();
    let facet = (first_diff - 84) / 50;
    let in_facet = (first_diff - 84) % 50;
    let component_byte = in_facet % 4;

    if in_facet < 48 {
        let component_offset = 84 + facet * 50 + (in_facet / 4) * 4;
        let reference_value = f32::from_le_bytes(
            reference[component_offset..component_offset + 4]
                .try_into()
                .ok()?,
        );
        let generated_value = f32::from_le_bytes(
            generated[component_offset..component_offset + 4]
                .try_into()
                .ok()?,
        );
        let reference_facet = stl_facet_summary(reference, facet)?;
        let generated_facet = stl_facet_summary(generated, facet)?;
        let facet_delta = stl_facet_delta_summary(reference, generated, facet)?;
        return Some(format!(
            "stl_context=facet {facet} {} byte {component_byte}; reference_f32={reference_value}, generated_f32={generated_value}; reference_facet={reference_facet}; generated_facet={generated_facet}; generated_minus_reference={facet_delta}{delta_scan}",
            stl_component_name(in_facet / 4)
        ));
    }

    let attribute_offset = 84 + facet * 50 + 48;
    let reference_value = u16::from_le_bytes(
        reference[attribute_offset..attribute_offset + 2]
            .try_into()
            .ok()?,
    );
    let generated_value = u16::from_le_bytes(
        generated[attribute_offset..attribute_offset + 2]
            .try_into()
            .ok()?,
    );
    Some(format!(
        "stl_context=facet {facet} attribute_byte_count byte {}; reference_u16={reference_value}, generated_u16={generated_value}{delta_scan}",
        in_facet - 48
    ))
}

fn looks_like_binary_stl(bytes: &[u8]) -> bool {
    if bytes.len() < 84 {
        return false;
    }
    let count = stl_facet_count(bytes).expect("slice length is 4");
    bytes.len() == 84 + count.saturating_mul(50)
}

fn stl_facet_count(bytes: &[u8]) -> Option<usize> {
    Some(u32::from_le_bytes(bytes.get(80..84)?.try_into().ok()?) as usize)
}

fn stl_component_name(component: usize) -> &'static str {
    const NAMES: [&str; 12] = [
        "normal.x",
        "normal.y",
        "normal.z",
        "vertex0.x",
        "vertex0.y",
        "vertex0.z",
        "vertex1.x",
        "vertex1.y",
        "vertex1.z",
        "vertex2.x",
        "vertex2.y",
        "vertex2.z",
    ];
    NAMES.get(component).copied().unwrap_or("unknown")
}

fn stl_facet_summary(bytes: &[u8], facet: usize) -> Option<String> {
    let base = 84 + facet * 50;
    Some(format!(
        "{{normal:{},vertices:[{},{},{}]}}",
        stl_vec3_summary(bytes, base)?,
        stl_vec3_summary(bytes, base + 12)?,
        stl_vec3_summary(bytes, base + 24)?,
        stl_vec3_summary(bytes, base + 36)?
    ))
}

fn stl_facet_delta_summary(reference: &[u8], generated: &[u8], facet: usize) -> Option<String> {
    let base = 84 + facet * 50;
    let vertex_delta0 = stl_vec3_delta(reference, generated, base + 12)?;
    let vertex_delta1 = stl_vec3_delta(reference, generated, base + 24)?;
    let vertex_delta2 = stl_vec3_delta(reference, generated, base + 36)?;
    let centroid_delta = vec3_average([vertex_delta0, vertex_delta1, vertex_delta2]);
    let residual0 = vec3_sub(vertex_delta0, centroid_delta);
    let residual1 = vec3_sub(vertex_delta1, centroid_delta);
    let residual2 = vec3_sub(vertex_delta2, centroid_delta);
    Some(format!(
        "{{normal:{},vertices:[{},{},{}],vertex_centroid:{},vertex_residuals:[{},{},{}]}}",
        stl_vec3_delta_summary(reference, generated, base)?,
        vec3_summary(vertex_delta0),
        vec3_summary(vertex_delta1),
        vec3_summary(vertex_delta2),
        vec3_summary(centroid_delta),
        vec3_summary(residual0),
        vec3_summary(residual1),
        vec3_summary(residual2)
    ))
}

fn stl_delta_scan_summary(reference: &[u8], generated: &[u8]) -> Option<String> {
    const CENTER_LIKE_RESIDUAL_RATIO_THRESHOLD: f32 = 0.1;
    const CENTER_LIKE_ABSOLUTE_RESIDUAL_EPSILON: f32 = 0.000_001;

    let reference_count = stl_facet_count(reference)?;
    let generated_count = stl_facet_count(generated)?;
    if reference_count != generated_count {
        return Some(format!(
            "{{reference_facet_count:{reference_count},generated_facet_count:{generated_count}}}"
        ));
    }

    let mut nonzero_vertex_delta_facets = 0usize;
    let mut nonzero_vertex_residual_facets = 0usize;
    let mut nonzero_normal_delta_facets = 0usize;
    let mut center_like_vertex_delta_facets = 0usize;
    let mut max_vertex_delta = None;
    let mut max_vertex_centroid = None;
    let mut max_vertex_residual = None;
    let mut max_normal = None;
    let mut max_vertex_residual_to_centroid_ratio = None;
    let mut vertex_centroid_deltas = Vec::new();

    for facet in 0..reference_count {
        let base = 84 + facet * 50;
        let normal_delta = stl_vec3_delta(reference, generated, base)?;
        let vertex_delta0 = stl_vec3_delta(reference, generated, base + 12)?;
        let vertex_delta1 = stl_vec3_delta(reference, generated, base + 24)?;
        let vertex_delta2 = stl_vec3_delta(reference, generated, base + 36)?;
        let vertex_deltas = [vertex_delta0, vertex_delta1, vertex_delta2];
        let centroid_delta = vec3_average(vertex_deltas);
        let residuals = [
            vec3_sub(vertex_delta0, centroid_delta),
            vec3_sub(vertex_delta1, centroid_delta),
            vec3_sub(vertex_delta2, centroid_delta),
        ];
        let has_vertex_delta = vertex_deltas.iter().any(|delta| !vec3_is_zero(*delta));
        let max_residual = residuals
            .iter()
            .copied()
            .max_by(|a, b| vec3_max_component_abs(*a).total_cmp(&vec3_max_component_abs(*b)))
            .unwrap_or([0.0, 0.0, 0.0]);
        let centroid_abs = vec3_max_component_abs(centroid_delta);
        let max_residual_abs = vec3_max_component_abs(max_residual);

        if has_vertex_delta {
            nonzero_vertex_delta_facets += 1;
            vertex_centroid_deltas.push(centroid_delta);
            let residual_to_centroid_ratio = if centroid_abs > 0.0 {
                max_residual_abs / centroid_abs
            } else if max_residual_abs == 0.0 {
                0.0
            } else {
                f32::INFINITY
            };
            update_residual_ratio_max(
                &mut max_vertex_residual_to_centroid_ratio,
                facet,
                residual_to_centroid_ratio,
                centroid_delta,
                max_residual,
            );
            if centroid_abs > 0.0
                && max_residual_abs
                    <= (centroid_abs * CENTER_LIKE_RESIDUAL_RATIO_THRESHOLD)
                        .max(CENTER_LIKE_ABSOLUTE_RESIDUAL_EPSILON)
            {
                center_like_vertex_delta_facets += 1;
            }
        }
        if residuals.iter().any(|delta| !vec3_is_zero(*delta)) {
            nonzero_vertex_residual_facets += 1;
        }
        if !vec3_is_zero(normal_delta) {
            nonzero_normal_delta_facets += 1;
        }

        update_delta_max(&mut max_normal, facet, normal_delta);
        update_delta_max(&mut max_vertex_centroid, facet, centroid_delta);
        for (vertex, delta) in vertex_deltas.iter().copied().enumerate() {
            update_vertex_delta_max(&mut max_vertex_delta, facet, vertex, delta);
        }
        for (vertex, delta) in residuals.iter().copied().enumerate() {
            update_vertex_delta_max(&mut max_vertex_residual, facet, vertex, delta);
        }
    }

    let center_fit = stl_center_fit_summary(reference, generated, &vertex_centroid_deltas)?;
    Some(format!(
        "{{facet_count:{reference_count},nonzero_vertex_delta_facets:{nonzero_vertex_delta_facets},nonzero_vertex_residual_facets:{nonzero_vertex_residual_facets},nonzero_normal_delta_facets:{nonzero_normal_delta_facets},center_like_vertex_delta_facets:{center_like_vertex_delta_facets},shape_like_vertex_delta_facets:{},center_like_residual_ratio_threshold:{CENTER_LIKE_RESIDUAL_RATIO_THRESHOLD},max_vertex_delta_abs:{},max_vertex_centroid_abs:{},max_vertex_residual_abs:{},max_vertex_residual_to_centroid_ratio:{},max_normal_abs:{},center_fit:{center_fit}}}",
        nonzero_vertex_delta_facets.saturating_sub(center_like_vertex_delta_facets),
        vertex_delta_max_summary(max_vertex_delta),
        delta_max_summary(max_vertex_centroid),
        vertex_delta_max_summary(max_vertex_residual),
        residual_ratio_max_summary(max_vertex_residual_to_centroid_ratio),
        delta_max_summary(max_normal)
    ))
}

fn stl_center_fit_summary(
    reference: &[u8],
    generated: &[u8],
    vertex_centroid_deltas: &[[f32; 3]],
) -> Option<String> {
    if vertex_centroid_deltas.is_empty() {
        return Some(
            "{real_vertex_delta_facets:0,median_vertex_centroid_delta:[0,0,0],mean_vertex_centroid_delta:[0,0,0],max_vertex_delta_after_median_center_abs:none,max_vertex_delta_after_mean_center_abs:none}"
                .to_string(),
        );
    }

    let median_center = vec3_component_median(vertex_centroid_deltas);
    let mean_center = vec3_component_mean(vertex_centroid_deltas);
    let reference_count = stl_facet_count(reference)?;
    let mut max_after_median_center = None;
    let mut max_after_mean_center = None;

    for facet in 0..reference_count {
        let base = 84 + facet * 50;
        let vertex_deltas = [
            stl_vec3_delta(reference, generated, base + 12)?,
            stl_vec3_delta(reference, generated, base + 24)?,
            stl_vec3_delta(reference, generated, base + 36)?,
        ];
        if !vertex_deltas.iter().any(|delta| !vec3_is_zero(*delta)) {
            continue;
        }
        for (vertex, delta) in vertex_deltas.iter().copied().enumerate() {
            update_vertex_delta_max(
                &mut max_after_median_center,
                facet,
                vertex,
                vec3_sub(delta, median_center),
            );
            update_vertex_delta_max(
                &mut max_after_mean_center,
                facet,
                vertex,
                vec3_sub(delta, mean_center),
            );
        }
    }

    Some(format!(
        "{{real_vertex_delta_facets:{},median_vertex_centroid_delta:{},mean_vertex_centroid_delta:{},max_vertex_delta_after_median_center_abs:{},max_vertex_delta_after_mean_center_abs:{}}}",
        vertex_centroid_deltas.len(),
        vec3_summary(median_center),
        vec3_summary(mean_center),
        vertex_delta_max_summary(max_after_median_center),
        vertex_delta_max_summary(max_after_mean_center)
    ))
}

fn stl_vec3_summary(bytes: &[u8], offset: usize) -> Option<String> {
    Some(vec3_summary(stl_vec3(bytes, offset)?))
}

fn stl_vec3_delta_summary(reference: &[u8], generated: &[u8], offset: usize) -> Option<String> {
    Some(vec3_summary(stl_vec3_delta(reference, generated, offset)?))
}

fn stl_vec3_delta(reference: &[u8], generated: &[u8], offset: usize) -> Option<[f32; 3]> {
    let reference = stl_vec3(reference, offset)?;
    let generated = stl_vec3(generated, offset)?;
    Some([
        generated[0] - reference[0],
        generated[1] - reference[1],
        generated[2] - reference[2],
    ])
}

fn stl_vec3(bytes: &[u8], offset: usize) -> Option<[f32; 3]> {
    Some([
        stl_f32(bytes, offset)?,
        stl_f32(bytes, offset + 4)?,
        stl_f32(bytes, offset + 8)?,
    ])
}

fn vec3_summary(values: [f32; 3]) -> String {
    format!("[{},{},{}]", values[0], values[1], values[2])
}

fn vec3_average(values: [[f32; 3]; 3]) -> [f32; 3] {
    [
        (values[0][0] + values[1][0] + values[2][0]) / 3.0,
        (values[0][1] + values[1][1] + values[2][1]) / 3.0,
        (values[0][2] + values[1][2] + values[2][2]) / 3.0,
    ]
}

fn vec3_sub(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vec3_component_median(values: &[[f32; 3]]) -> [f32; 3] {
    [
        component_median(values, 0),
        component_median(values, 1),
        component_median(values, 2),
    ]
}

fn component_median(values: &[[f32; 3]], axis: usize) -> f32 {
    let mut sorted = values.iter().map(|value| value[axis]).collect::<Vec<_>>();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

fn vec3_component_mean(values: &[[f32; 3]]) -> [f32; 3] {
    let mut sum = [0.0f64; 3];
    for value in values {
        sum[0] += value[0] as f64;
        sum[1] += value[1] as f64;
        sum[2] += value[2] as f64;
    }
    [
        (sum[0] / values.len() as f64) as f32,
        (sum[1] / values.len() as f64) as f32,
        (sum[2] / values.len() as f64) as f32,
    ]
}

#[derive(Clone, Copy)]
struct DeltaMax {
    facet: usize,
    max_component_abs: f32,
    delta: [f32; 3],
}

#[derive(Clone, Copy)]
struct VertexDeltaMax {
    facet: usize,
    vertex: usize,
    max_component_abs: f32,
    delta: [f32; 3],
}

#[derive(Clone, Copy)]
struct ResidualRatioMax {
    facet: usize,
    ratio: f32,
    centroid_delta: [f32; 3],
    max_residual_delta: [f32; 3],
}

fn update_delta_max(max: &mut Option<DeltaMax>, facet: usize, delta: [f32; 3]) {
    let max_component_abs = vec3_max_component_abs(delta);
    if max
        .map(|current| max_component_abs > current.max_component_abs)
        .unwrap_or(true)
    {
        *max = Some(DeltaMax {
            facet,
            max_component_abs,
            delta,
        });
    }
}

fn update_vertex_delta_max(
    max: &mut Option<VertexDeltaMax>,
    facet: usize,
    vertex: usize,
    delta: [f32; 3],
) {
    let max_component_abs = vec3_max_component_abs(delta);
    if max
        .map(|current| max_component_abs > current.max_component_abs)
        .unwrap_or(true)
    {
        *max = Some(VertexDeltaMax {
            facet,
            vertex,
            max_component_abs,
            delta,
        });
    }
}

fn update_residual_ratio_max(
    max: &mut Option<ResidualRatioMax>,
    facet: usize,
    ratio: f32,
    centroid_delta: [f32; 3],
    max_residual_delta: [f32; 3],
) {
    if max.map(|current| ratio > current.ratio).unwrap_or(true) {
        *max = Some(ResidualRatioMax {
            facet,
            ratio,
            centroid_delta,
            max_residual_delta,
        });
    }
}

fn delta_max_summary(max: Option<DeltaMax>) -> String {
    max.map(|max| {
        format!(
            "{{facet:{},max_component_abs:{},delta:{}}}",
            max.facet,
            max.max_component_abs,
            vec3_summary(max.delta)
        )
    })
    .unwrap_or_else(|| "none".to_string())
}

fn vertex_delta_max_summary(max: Option<VertexDeltaMax>) -> String {
    max.map(|max| {
        format!(
            "{{facet:{},vertex:{},max_component_abs:{},delta:{}}}",
            max.facet,
            max.vertex,
            max.max_component_abs,
            vec3_summary(max.delta)
        )
    })
    .unwrap_or_else(|| "none".to_string())
}

fn residual_ratio_max_summary(max: Option<ResidualRatioMax>) -> String {
    max.map(|max| {
        format!(
            "{{facet:{},ratio:{},centroid_delta:{},max_residual_delta:{}}}",
            max.facet,
            max.ratio,
            vec3_summary(max.centroid_delta),
            vec3_summary(max.max_residual_delta)
        )
    })
    .unwrap_or_else(|| "none".to_string())
}

fn vec3_is_zero(values: [f32; 3]) -> bool {
    values[0] == 0.0 && values[1] == 0.0 && values[2] == 0.0
}

fn vec3_max_component_abs(values: [f32; 3]) -> f32 {
    values[0].abs().max(values[1].abs()).max(values[2].abs())
}

fn stl_f32(bytes: &[u8], offset: usize) -> Option<f32> {
    Some(f32::from_le_bytes(
        bytes.get(offset..offset + 4)?.try_into().ok()?,
    ))
}
