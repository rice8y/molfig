use std::env;
use std::fs;
use std::ops::Range;
use std::process::ExitCode;
use std::time::Instant;

use molfig::{
    convert_to_obj, convert_to_ply, convert_to_stl, diff_bytes_with_options, diff_text,
    stl_export_facet_context, stl_export_facet_context_timed, stl_facet_semantic_context,
    stl_facet_semantic_context_with_vertex_offset,
    stl_facet_semantic_context_with_vertex_offset_timed, DiffBytesOptions,
};

struct CliArgs {
    json: bool,
    timings: bool,
    generated_in: Option<String>,
    generated_out: Option<String>,
    stl_delta_scan: bool,
    stl_semantic_context: bool,
    stl_facet_range: Option<Range<usize>>,
    stl_facet_context: Option<usize>,
    stl_export_facet_context: Option<usize>,
    stl_vertex_offset: Option<[f64; 3]>,
    format: String,
    input_path: String,
    options_path: String,
    reference_path: Option<String>,
}

fn main() -> ExitCode {
    match run() {
        Ok(message) => {
            println!("{message}");
            ExitCode::SUCCESS
        }
        Err(message) => {
            eprintln!("{message}");
            ExitCode::from(1)
        }
    }
}

fn run() -> Result<String, String> {
    let started = Instant::now();
    let mut last_timing = started;
    let args = parse_args(env::args().skip(1))?;

    let format = args.format.as_str();
    if args.generated_in.is_some() && args.generated_out.is_some() {
        return Err("--generated-in and --generated-out are mutually exclusive".to_string());
    }
    if format != "stl" && args.stl_facet_range.is_some() {
        return Err("--stl-facet-range is only valid for stl".to_string());
    }
    let input = fs::read(&args.input_path).map_err(|err| format!("failed to read input: {err}"))?;
    let options = fs::read(&args.options_path)
        .map_err(|err| format!("failed to read options JSON: {err}"))?;
    timing_checkpoint(
        args.timings,
        "read-input-options",
        started,
        &mut last_timing,
    );
    if let Some(facet_index) = args.stl_facet_context {
        if format != "stl" {
            return Err("--stl-facet-context is only valid for stl".to_string());
        }
        if args.stl_export_facet_context.is_some() {
            return Err(
                "--stl-facet-context and --stl-export-facet-context are mutually exclusive"
                    .to_string(),
            );
        }
        if let Some(vertex_offset) = args.stl_vertex_offset {
            let context = if args.timings {
                stl_facet_semantic_context_with_vertex_offset_timed(
                    &input,
                    &options,
                    facet_index,
                    vertex_offset,
                    |label| timing_checkpoint(true, label, started, &mut last_timing),
                )?
            } else {
                stl_facet_semantic_context_with_vertex_offset(
                    &input,
                    &options,
                    facet_index,
                    vertex_offset,
                )?
            };
            timing_checkpoint(args.timings, "stl-facet-context", started, &mut last_timing);
            return Ok(context);
        }
        let context = stl_facet_semantic_context(&input, &options, facet_index)?;
        timing_checkpoint(args.timings, "stl-facet-context", started, &mut last_timing);
        return Ok(context);
    }
    if let Some(facet_index) = args.stl_export_facet_context {
        if format != "stl" {
            return Err("--stl-export-facet-context is only valid for stl".to_string());
        }
        if args.stl_vertex_offset.is_some() {
            return Err("--stl-vertex-offset is only valid with --stl-facet-context".to_string());
        }
        let context = if args.timings {
            stl_export_facet_context_timed(&input, &options, facet_index, |label| {
                timing_checkpoint(true, label, started, &mut last_timing)
            })?
        } else {
            stl_export_facet_context(&input, &options, facet_index)?
        };
        timing_checkpoint(
            args.timings,
            "stl-export-facet-context",
            started,
            &mut last_timing,
        );
        return Ok(context);
    }
    let reference_path = args.reference_path.as_deref().ok_or_else(usage)?;
    let reference = fs::read(reference_path)
        .map_err(|err| format!("failed to read reference export: {err}"))?;
    timing_checkpoint(args.timings, "read-reference", started, &mut last_timing);

    let generated = if let Some(path) = &args.generated_in {
        let generated =
            fs::read(path).map_err(|err| format!("failed to read generated export: {err}"))?;
        timing_checkpoint(args.timings, "read-generated", started, &mut last_timing);
        generated
    } else {
        let generated = match format {
            "obj" => convert_to_obj(&input, &options)?,
            "ply" => convert_to_ply(&input, &options)?,
            "stl" => convert_to_stl(&input, &options)?,
            _ => return Err(usage()),
        };
        timing_checkpoint(args.timings, "generate-export", started, &mut last_timing);
        generated
    };
    if let Some(path) = &args.generated_out {
        fs::write(path, &generated)
            .map_err(|err| format!("failed to write generated export: {err}"))?;
    }

    let report = if matches!(format, "obj" | "ply") {
        let reference = std::str::from_utf8(&reference)
            .map_err(|_| "reference export is not valid UTF-8 text".to_string())?;
        let generated = std::str::from_utf8(&generated)
            .map_err(|_| "generated export is not valid UTF-8 text".to_string())?;
        diff_text(reference, generated, format)
    } else {
        diff_bytes_with_options(
            &reference,
            &generated,
            format,
            DiffBytesOptions {
                stl_delta_scan: args.stl_delta_scan,
                stl_facet_range: args.stl_facet_range.clone(),
            },
        )
    };

    let stl_semantic_context = if format == "stl" && !report.passed && args.stl_semantic_context {
        semantic_context_for_stl_diff(&input, &options, &report)
    } else {
        None
    };

    let message = if args.json {
        diff_report_json(format, &report, stl_semantic_context.as_deref())
    } else {
        append_stl_semantic_context(report.message.clone(), stl_semantic_context.as_deref())
    };

    if report.passed {
        Ok(message)
    } else {
        Err(message)
    }
}

fn parse_args(args: impl Iterator<Item = String>) -> Result<CliArgs, String> {
    let mut json = false;
    let mut timings = false;
    let mut generated_in = None;
    let mut generated_out = None;
    let mut stl_delta_scan = true;
    let mut stl_semantic_context = true;
    let mut stl_facet_range = None;
    let mut stl_facet_context = None;
    let mut stl_export_facet_context = None;
    let mut stl_vertex_offset = None;
    let mut positional = Vec::new();
    let mut args = args.peekable();

    while let Some(arg) = args.next() {
        if arg == "--json" {
            json = true;
        } else if arg == "--timings" {
            timings = true;
        } else if arg == "--generated-in" {
            let value = args
                .next()
                .ok_or_else(|| "--generated-in requires a path".to_string())?;
            generated_in = Some(value);
        } else if arg == "--generated-out" {
            let value = args
                .next()
                .ok_or_else(|| "--generated-out requires a path".to_string())?;
            generated_out = Some(value);
        } else if arg == "--stl-first-diff-only" {
            stl_delta_scan = false;
            stl_semantic_context = false;
        } else if arg == "--skip-stl-delta-scan" {
            stl_delta_scan = false;
        } else if arg == "--no-stl-semantic-context" {
            stl_semantic_context = false;
        } else if arg == "--stl-facet-range" {
            let value = args
                .next()
                .ok_or_else(|| "--stl-facet-range requires start..end".to_string())?;
            stl_facet_range = Some(parse_usize_range(&value, "--stl-facet-range")?);
        } else if arg == "--stl-facet-context" {
            let value = args
                .next()
                .ok_or_else(|| "--stl-facet-context requires a facet index".to_string())?;
            stl_facet_context =
                Some(value.parse::<usize>().map_err(|_| {
                    "--stl-facet-context requires a non-negative integer".to_string()
                })?);
        } else if arg == "--stl-export-facet-context" {
            let value = args
                .next()
                .ok_or_else(|| "--stl-export-facet-context requires a facet index".to_string())?;
            stl_export_facet_context = Some(value.parse::<usize>().map_err(|_| {
                "--stl-export-facet-context requires a non-negative integer".to_string()
            })?);
        } else if arg == "--stl-vertex-offset" {
            let value = args
                .next()
                .ok_or_else(|| "--stl-vertex-offset requires x,y,z".to_string())?;
            stl_vertex_offset = Some(parse_f64_triplet(&value, "--stl-vertex-offset")?);
        } else if arg == "--help" || arg == "-h" {
            return Err(usage());
        } else {
            positional.push(arg);
        }
    }

    let (format, input_path, options_path, reference_path) = match positional.len() {
        3 if stl_facet_context.is_some() || stl_export_facet_context.is_some() => (
            positional.remove(0),
            positional.remove(0),
            positional.remove(0),
            None,
        ),
        4 => (
            positional.remove(0),
            positional.remove(0),
            positional.remove(0),
            Some(positional.remove(0)),
        ),
        _ => return Err(usage()),
    };
    Ok(CliArgs {
        json,
        timings,
        generated_in,
        generated_out,
        stl_delta_scan,
        stl_semantic_context,
        stl_facet_range,
        stl_facet_context,
        stl_export_facet_context,
        stl_vertex_offset,
        format,
        input_path,
        options_path,
        reference_path,
    })
}

fn parse_usize_range(value: &str, flag: &str) -> Result<Range<usize>, String> {
    let (start, end) = value
        .split_once("..")
        .ok_or_else(|| format!("{flag} requires start..end"))?;
    let start = start
        .trim()
        .parse::<usize>()
        .map_err(|_| format!("{flag} requires non-negative integer bounds"))?;
    let end = end
        .trim()
        .parse::<usize>()
        .map_err(|_| format!("{flag} requires non-negative integer bounds"))?;
    if start >= end {
        return Err(format!("{flag} requires start < end"));
    }
    Ok(start..end)
}

fn timing_checkpoint(enabled: bool, label: &str, started: Instant, last: &mut Instant) {
    if !enabled {
        return;
    }
    let now = Instant::now();
    eprintln!(
        "timing {label}: +{:.3}s total={:.3}s",
        now.duration_since(*last).as_secs_f64(),
        now.duration_since(started).as_secs_f64()
    );
    *last = now;
}

fn parse_f64_triplet(value: &str, flag: &str) -> Result<[f64; 3], String> {
    let parts = value
        .split(',')
        .map(str::trim)
        .map(|part| {
            part.parse::<f64>()
                .map_err(|_| format!("{flag} requires comma-separated finite numbers"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let [x, y, z]: [f64; 3] = parts
        .try_into()
        .map_err(|_| format!("{flag} requires exactly three comma-separated numbers"))?;
    if !x.is_finite() || !y.is_finite() || !z.is_finite() {
        return Err(format!("{flag} requires comma-separated finite numbers"));
    }
    Ok([x, y, z])
}

fn semantic_context_for_stl_diff(
    input: &[u8],
    options: &[u8],
    report: &molfig::DiffReport,
) -> Option<String> {
    let first_byte = report
        .details
        .iter()
        .find_map(|(key, value)| (key == "first_byte").then_some(value))?
        .parse::<usize>()
        .ok()?;
    if first_byte < 84 {
        return None;
    }
    let facet_index = (first_byte - 84) / 50;
    Some(
        match stl_facet_semantic_context(input, options, facet_index) {
            Ok(context) => context,
            Err(error) => format!(
                "{{\"found\":false,\"stl_facet\":{},\"error\":{}}}",
                facet_index,
                json_string(&error)
            ),
        },
    )
}

fn append_stl_semantic_context(message: String, context: Option<&str>) -> String {
    match context {
        Some(context) => format!("{message}; stl_semantic_context={context}"),
        None => message,
    }
}

fn diff_report_json(
    format: &str,
    report: &molfig::DiffReport,
    stl_semantic_context: Option<&str>,
) -> String {
    let mut details = report
        .details
        .iter()
        .map(|(key, value)| format!("{}:{}", json_string(key), json_string(value)))
        .collect::<Vec<_>>();
    if let Some(context) = stl_semantic_context {
        details.push(format!(
            "{}:{}",
            json_string("stl_semantic_context"),
            json_string(context)
        ));
    }
    let details = details.join(",");
    format!(
        "{{\"format\":{},\"passed\":{},\"message\":{},\"details\":{{{}}}}}",
        json_string(format),
        report.passed,
        json_string(&report.message),
        details
    )
}

fn json_string(value: &str) -> String {
    let mut out = String::from("\"");
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            ch if ch.is_control() => {
                out.push_str(&format!("\\u{:04x}", ch as u32));
            }
            ch => out.push(ch),
        }
    }
    out.push('"');
    out
}

fn usage() -> String {
    "usage: molfig-diff [--json] [--timings] [--generated-in <path>|--generated-out <path>] [--stl-first-diff-only] [--skip-stl-delta-scan] [--no-stl-semantic-context] [--stl-facet-range start..end] <obj|ply|stl> <input.pdb|cif|bcif> <options.json> <reference-export>\n       molfig-diff [--timings] --stl-facet-context <facet> [--stl-vertex-offset x,y,z] stl <input.pdb|cif|bcif> <options.json>\n       molfig-diff [--timings] --stl-export-facet-context <facet> stl <input.pdb|cif|bcif> <options.json>"
        .to_string()
}
