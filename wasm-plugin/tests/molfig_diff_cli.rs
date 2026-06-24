use std::fs;
use std::process::Command;

#[test]
fn molfig_diff_stl_cli_emits_aggregate_delta_scan() {
    let options = br#"{"format":"pdb","representation":"cartoon","assembly":null}"#;
    let input = include_bytes!("fixtures/pdb/tiny-peptide.pdb");
    let generated = molfig::convert_to_stl(input, options).expect("fixture STL");
    assert!(
        generated.len() > 96,
        "fixture STL should contain facet floats"
    );

    let mut reference = generated.clone();
    let normal_x = f32::from_le_bytes(reference[84..88].try_into().unwrap());
    reference[84..88].copy_from_slice(&(normal_x + 1.0).to_le_bytes());

    let temp_dir = std::env::temp_dir().join(format!(
        "molfig-diff-cli-{}-{}",
        std::process::id(),
        unique_test_suffix()
    ));
    fs::create_dir_all(&temp_dir).expect("create temp dir");
    let input_path = temp_dir.join("input.pdb");
    let options_path = temp_dir.join("options.json");
    let reference_path = temp_dir.join("reference.stl");
    let generated_path = temp_dir.join("generated.stl");
    let generated_in_path = temp_dir.join("generated-in.stl");
    fs::write(&input_path, input).expect("write input");
    fs::write(&options_path, options).expect("write options");
    fs::write(&reference_path, &reference).expect("write reference");
    fs::write(&generated_in_path, &generated).expect("write generated-in");

    let output = Command::new(env!("CARGO_BIN_EXE_molfig-diff"))
        .arg("stl")
        .arg(&input_path)
        .arg(&options_path)
        .arg(&reference_path)
        .output()
        .expect("run molfig-diff");
    let first_diff_only_output = Command::new(env!("CARGO_BIN_EXE_molfig-diff"))
        .arg("--stl-first-diff-only")
        .arg("--generated-in")
        .arg(&generated_in_path)
        .arg("stl")
        .arg(&input_path)
        .arg(&options_path)
        .arg(&reference_path)
        .output()
        .expect("run molfig-diff --stl-first-diff-only --generated-in");
    let facet_range_fail_output = Command::new(env!("CARGO_BIN_EXE_molfig-diff"))
        .arg("--stl-first-diff-only")
        .arg("--stl-facet-range")
        .arg("0..1")
        .arg("--generated-in")
        .arg(&generated_in_path)
        .arg("stl")
        .arg(&input_path)
        .arg(&options_path)
        .arg(&reference_path)
        .output()
        .expect("run molfig-diff --stl-facet-range failing range");
    let facet_range_pass_output = Command::new(env!("CARGO_BIN_EXE_molfig-diff"))
        .arg("--stl-facet-range")
        .arg("1..2")
        .arg("--generated-in")
        .arg(&generated_in_path)
        .arg("stl")
        .arg(&input_path)
        .arg(&options_path)
        .arg(&reference_path)
        .output()
        .expect("run molfig-diff --stl-facet-range passing range");
    let json_output = Command::new(env!("CARGO_BIN_EXE_molfig-diff"))
        .arg("--json")
        .arg("--generated-out")
        .arg(&generated_path)
        .arg("stl")
        .arg(&input_path)
        .arg(&options_path)
        .arg(&reference_path)
        .output()
        .expect("run molfig-diff --json");
    let facet_context_output = Command::new(env!("CARGO_BIN_EXE_molfig-diff"))
        .arg("--stl-facet-context")
        .arg("0")
        .arg("stl")
        .arg(&input_path)
        .arg(&options_path)
        .output()
        .expect("run molfig-diff --stl-facet-context");
    let offset_facet_context_output = Command::new(env!("CARGO_BIN_EXE_molfig-diff"))
        .arg("--stl-facet-context")
        .arg("0")
        .arg("--stl-vertex-offset")
        .arg("1,2,3")
        .arg("stl")
        .arg(&input_path)
        .arg(&options_path)
        .output()
        .expect("run molfig-diff --stl-facet-context --stl-vertex-offset");
    let export_facet_context_output = Command::new(env!("CARGO_BIN_EXE_molfig-diff"))
        .arg("--stl-export-facet-context")
        .arg("0")
        .arg("stl")
        .arg(&input_path)
        .arg(&options_path)
        .output()
        .expect("run molfig-diff --stl-export-facet-context");
    let written_generated = fs::read(&generated_path).expect("generated export was written");

    let _ = fs::remove_dir_all(&temp_dir);
    assert_eq!(
        written_generated, generated,
        "--generated-out should save the generated export even when the diff fails"
    );
    assert!(!output.status.success(), "mutated reference should fail");
    let stderr = String::from_utf8(output.stderr).expect("stderr utf8");
    assert!(stderr.contains("FAIL stl: first difference"));
    assert!(stderr.contains("stl_delta_scan={facet_count:"));
    assert!(stderr.contains("nonzero_normal_delta_facets:1"));
    assert!(stderr.contains("max_normal_abs:{facet:0"));
    assert!(stderr.contains("center_fit:{real_vertex_delta_facets:"));
    assert!(stderr.contains("stl_semantic_context={\"found\":true"));
    assert!(stderr.contains("\"stl_facet\":0"));
    assert!(stderr.contains("\"vertex_offset\":["));
    assert!(stderr.contains("\"stl_vertices\":["));
    assert!(stderr.contains("\"stl_vertex_bits\":["));
    assert!(stderr.contains("\"target_face\":{\"indices\""));

    assert!(
        !first_diff_only_output.status.success(),
        "first-diff-only should still fail on the changed byte"
    );
    let first_diff_only_stderr =
        String::from_utf8(first_diff_only_output.stderr).expect("first-diff stderr utf8");
    assert!(first_diff_only_stderr.contains("FAIL stl: first difference"));
    assert!(first_diff_only_stderr.contains("stl_context=facet 0"));
    assert!(
        !first_diff_only_stderr.contains("stl_delta_scan="),
        "--stl-first-diff-only should skip the aggregate STL scan"
    );
    assert!(
        !first_diff_only_stderr.contains("stl_semantic_context="),
        "--stl-first-diff-only should skip semantic context reconstruction"
    );

    assert!(
        !facet_range_fail_output.status.success(),
        "range containing the changed facet should fail"
    );
    let facet_range_fail_stderr =
        String::from_utf8(facet_range_fail_output.stderr).expect("range fail stderr utf8");
    assert!(facet_range_fail_stderr.contains("first difference in STL facet range 0..1"));
    assert!(facet_range_fail_stderr.contains("stl_context=facet 0"));
    assert!(
        facet_range_pass_output.status.success(),
        "range excluding the changed facet should pass"
    );
    let facet_range_pass_stdout =
        String::from_utf8(facet_range_pass_output.stdout).expect("range pass stdout utf8");
    assert!(facet_range_pass_stdout.contains("PASS stl: STL facet range 1..2 byte-for-byte match"));

    assert!(
        !json_output.status.success(),
        "mutated reference should fail"
    );
    let json_stderr = String::from_utf8(json_output.stderr).expect("json stderr utf8");
    assert!(json_stderr.contains(r#""format":"stl""#));
    assert!(json_stderr.contains(r#""passed":false"#));
    assert!(json_stderr.contains(r#""details":{"kind":"bytes""#));
    assert!(json_stderr.contains(r#""first_byte":"#));
    assert!(json_stderr.contains(r#""stl_delta_scan":"{facet_count:"#));
    assert!(json_stderr.contains(r#""stl_context":"stl_context=facet 0"#));
    assert!(json_stderr.contains(r#""stl_semantic_context":"{\"found\":true"#));
    assert!(json_stderr.contains(r#"\"stl_facet\":0"#));
    assert!(json_stderr.contains(r#"\"stl_vertex_bits\":["#));
    assert!(json_stderr.contains("stl_delta_scan={facet_count:"));
    assert!(json_stderr.contains("max_normal_abs:{facet:0"));
    assert!(json_stderr.contains("center_fit:{real_vertex_delta_facets:"));

    assert!(
        facet_context_output.status.success(),
        "facet context should not require a reference export"
    );
    let facet_context_stdout =
        String::from_utf8(facet_context_output.stdout).expect("facet context stdout utf8");
    assert!(facet_context_stdout.contains(r#""stl_facet":0"#));
    assert!(facet_context_stdout.contains(r#""vertex_offset":["#));
    assert!(facet_context_stdout.contains(r#""stl_vertex_bits":["#));
    assert!(facet_context_stdout.contains(r#""target_face":{"indices""#));
    assert!(
        offset_facet_context_output.status.success(),
        "offset facet context should bypass export-center calculation"
    );
    let offset_facet_context_stdout =
        String::from_utf8(offset_facet_context_output.stdout).expect("offset context stdout utf8");
    assert!(offset_facet_context_stdout.contains(
        r#""vertex_offset":[1.00000000000000000,2.00000000000000000,3.00000000000000000]"#
    ));
    assert!(
        export_facet_context_output.status.success(),
        "export facet context should not require a reference export"
    );
    let export_facet_context_stdout =
        String::from_utf8(export_facet_context_output.stdout).expect("export context stdout utf8");
    assert!(export_facet_context_stdout.contains(r#""stl_facet":0"#));
    assert!(export_facet_context_stdout.contains(r#""export_center":["#));
    assert!(export_facet_context_stdout.contains(r#""export_box_min":["#));
    assert!(export_facet_context_stdout.contains(r#""export_box_max":["#));
    assert!(export_facet_context_stdout.contains(r#""export_box_min_indices":["#));
    assert!(export_facet_context_stdout.contains(r#""export_box_max_indices":["#));
    assert!(export_facet_context_stdout.contains(r#""visible_sphere_report":{"component_count":"#));
    assert!(export_facet_context_stdout.contains(r#""components":["#));
    assert!(export_facet_context_stdout.contains(r#""scene":"#));
    assert!(export_facet_context_stdout.contains(r#""vertex_offset":["#));
    assert!(export_facet_context_stdout.contains(r#""sparse_slot_has_face":true"#));
    assert!(export_facet_context_stdout.contains(r#""stl_normal_bits":["#));
    assert!(export_facet_context_stdout.contains(r#""stl_vertex_bits":[["#));
    assert!(export_facet_context_stdout.contains(r#""target_face":{"indices""#));
}

fn unique_test_suffix() -> u128 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("system time")
        .as_nanos()
}
