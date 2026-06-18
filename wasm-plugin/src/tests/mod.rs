use super::*;
use crate::api::validate_mesh_for_export;
use crate::export::{export_obj, export_ply, export_stl};
use crate::mesh::{
    add_oriented_ribbon, add_profile_tube_for_test, add_ribbon_for_test, add_sheet_for_test,
    add_tube_path_for_test, build_mesh, build_mesh_with_visible_bounding_sphere,
    build_render_objects, coarse_polymer_trace_iterator_reference_json, interpolate_curve_segment,
    interpolate_sizes, polymer_trace_iterator_reference_json,
    polymer_trace_iterator_reference_json_with_helix_orientation, render_object_span_summary_json,
    render_object_summary_json, representation_summary_json, CurveSegmentControls,
    CurveSegmentState, DVec3, PolymerTraceSegmentKind, RenderObject, TestTubeProfile,
};
use crate::model::{
    AtomSiteColumnPresence, AtomicUnitKind, Axes3D, Bond, CoarseElementKind, Face,
    InterUnitBondEdge, InterUnitBondInfo, InterUnitBondProps, Mesh, MoleculeType, PolymerType,
    SecondaryStructureElement, SecondaryStructureType,
};
use crate::options::{ColorTheme, PolymerProfile, Representation, VisualQuality};
use crate::parser::{
    compose_operator_transforms, expand_oper_expression, parse_molecule_with_options, parse_pdb,
    ColumnData,
};

mod molstar_model_parity;

const PDB: &[u8] = b"ATOM      1  N   GLY A   1      11.104  13.207   2.100  1.00 10.00           N\nATOM      2  CA  GLY A   1      12.560  13.207   2.100  1.00 10.00           C\nATOM      3  C   GLY A   1      13.010  14.640   2.100  1.00 10.00           C\nEND\n";

const CIF: &[u8] = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 N N GLY A 1 11.104 13.207 2.100\nATOM 2 C CA GLY A 1 12.560 13.207 2.100\n#\n";

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("wasm-plugin should live under the repository root")
        .to_path_buf()
}

fn read_repo_file_if_present(path: &str) -> Option<String> {
    std::fs::read_to_string(repo_root().join(path)).ok()
}

fn read_repo_bytes_if_present(path: &str) -> Option<Vec<u8>> {
    std::fs::read(repo_root().join(path)).ok()
}

fn read_internal_doc(name: &str) -> Option<String> {
    let path = repo_root().join("dev-docs").join(name);
    std::fs::read_to_string(path).ok()
}

fn read_molstar_source(path: &str) -> Option<String> {
    read_repo_file_if_present(&format!("wasm-plugin/artifacts/molstar/src/{path}"))
}

#[test]
fn parses_pdb_atoms() {
    let mol = parse_molecule(PDB, InputFormat::Pdb).unwrap();
    assert_eq!(mol.atoms.len(), 3);
    assert_eq!(mol.atoms[0].element, "N");
    assert_eq!(mol.source_data.kind, "pdb");
    assert!(!mol.bonds.is_empty());
}

#[test]
fn parses_cif_atoms() {
    let mol = parse_molecule(CIF, InputFormat::Cif).unwrap();
    assert_eq!(mol.atoms.len(), 2);
    assert_eq!(mol.atoms[1].name, "CA");
    assert_eq!(mol.source_data.kind, "mmCIF");
    assert_eq!(mol.source_data.name, "demo");
    assert!(mol
        .source_data
        .categories
        .iter()
        .any(|category| category.name == "atom_site"
            && category.row_count == 2
            && category.column_count == 10));
}

#[test]
fn parses_cif_single_row_categories_like_loop_tables() {
    let cif = b"data_demo\n_entry.id DEMO\n_exptl.method 'ELECTRON MICROSCOPY'\n_entity.id 1\n_entity.type polymer\n_entity.pdbx_description 'single-row entity'\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 0.0 0.0 0.0\n#\n";
    let mol = parse_molecule(cif, InputFormat::Cif).unwrap();

    assert_eq!(mol.entries, vec![Entry { id: "DEMO".into() }]);
    assert_eq!(
        mol.experiments,
        vec![Experiment {
            method: "ELECTRON MICROSCOPY".into()
        }]
    );
    assert_eq!(mol.entities.len(), 1);
    assert_eq!(mol.entities[0].id, "1");
    assert_eq!(mol.entities[0].type_name, "polymer");
    assert_eq!(mol.entities[0].description, "single-row entity");
    assert_eq!(mol.atoms.len(), 1);
}

#[test]
fn cif_single_row_missing_values_and_duplicates_match_molstar_rules() {
    let missing = b"data_demo\n_entry.id\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 0.0 0.0 0.0\n#\n";
    assert_eq!(
        parse_molecule(missing, InputFormat::Cif).unwrap_err(),
        "Expected value."
    );

    let duplicate = b"data_demo\n_entry.id A\n_entry.id B\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 0.0 0.0 0.0\n#\n";
    let mol = parse_molecule(duplicate, InputFormat::Cif).unwrap();
    assert_eq!(mol.entries, vec![Entry { id: "B".into() }]);
}

#[test]
fn cif_single_row_parser_keeps_first_data_block_separate() {
    let cif = b"data_first\n_entry.id FIRST\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 0.0 0.0 0.0\n#\ndata_second\n_entry.id SECOND\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 2 C CA ALA B 1 1.0 0.0 0.0\n#\n";
    let mol = parse_molecule(cif, InputFormat::Cif).unwrap();

    assert_eq!(mol.entries, vec![Entry { id: "FIRST".into() }]);
    assert_eq!(mol.atoms.len(), 1);
    assert_eq!(mol.atoms[0].chain, "A");
}

#[test]
fn molstar_reference_commit_is_pinned() {
    assert_eq!(
        MOLSTAR_REFERENCE_COMMIT,
        "1b8117d3f10f7c978aabb5a0d3d47370635aefe4"
    );
}

#[test]
fn molstar_parity_checklist_progress_matches_checkbox_counts() {
    let Some(checklist) = read_internal_doc("molstar-parity-checklist.md") else {
        eprintln!("skipping internal Mol* parity checklist audit; dev-docs is absent");
        return;
    };
    let progress_line = checklist
        .lines()
        .find(|line| line.starts_with("Progress: "))
        .expect("missing Mol* parity checklist progress line");
    let progress = progress_line
        .strip_prefix("Progress: ")
        .expect("missing progress prefix");
    let documented_percent = progress
        .split_once('%')
        .expect("missing progress percent suffix")
        .0
        .parse::<f64>()
        .expect("invalid progress percent");
    let documented_counts = progress
        .split_once('(')
        .expect("missing progress count prefix")
        .1
        .split_once(' ')
        .expect("missing progress count suffix")
        .0
        .split_once('/')
        .expect("missing checked/total progress separator");
    let documented_checked = documented_counts
        .0
        .parse::<usize>()
        .expect("invalid checked progress count");
    let documented_total = documented_counts
        .1
        .parse::<usize>()
        .expect("invalid total progress count");

    let checked = checklist
        .lines()
        .filter(|line| line.trim_start().starts_with("- [x]"))
        .count();
    let unchecked = checklist
        .lines()
        .filter(|line| line.trim_start().starts_with("- [ ]"))
        .count();
    let total = checked + unchecked;
    let computed_percent = (checked as f64 * 1000.0 / total as f64).round() / 10.0;

    assert_eq!(documented_checked, checked);
    assert_eq!(documented_total, total);
    assert_eq!(
        format!("{documented_percent:.1}"),
        format!("{computed_percent:.1}")
    );
}

#[test]
fn molstar_parity_evidence_manifest_is_well_formed() {
    let entries = parity_evidence_entries();
    let test_sources = parity_evidence_test_sources();
    let repo = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut domains = std::collections::BTreeSet::new();

    assert!(
        entries.len() >= 20,
        "parity evidence should cover the main checked domains"
    );

    for entry in &entries {
        assert!(
            domains.insert(entry.domain.clone()),
            "duplicate parity evidence domain: {}",
            entry.domain
        );
        assert!(
            !entry.unit_tests.is_empty(),
            "missing unit tests for {}",
            entry.domain
        );
        assert!(
            !entry.fixtures.is_empty(),
            "missing fixtures for {}",
            entry.domain
        );
        assert!(
            !entry.notes.trim().is_empty(),
            "missing notes for {}",
            entry.domain
        );

        for test in &entry.unit_tests {
            assert!(
                test_sources.contains(&format!("fn {test}(")),
                "unit test listed for {} was not found in source index: {test}",
                entry.domain
            );
        }
        for fixture in &entry.fixtures {
            assert_manifest_token_exists(repo, &test_sources, &entry.domain, "fixture", fixture);
        }
        for reference in &entry.reference_comparisons {
            assert_manifest_token_exists(
                repo,
                &test_sources,
                &entry.domain,
                "reference comparison",
                reference,
            );
        }
    }
}

#[test]
fn molstar_parity_evidence_manifest_covers_export_and_reference_domains() {
    let entries = parity_evidence_entries();
    let by_domain = entries
        .iter()
        .map(|entry| (entry.domain.as_str(), entry))
        .collect::<std::collections::BTreeMap<_, _>>();

    for domain in [
        "parser-input-paths",
        "binary-cif-direct-path",
        "assembly-units-and-operators",
        "chem-comp-struct-conn-index-pair-bonds",
        "coarse-ihm-model-units",
        "polymer-trace-cyclic-ranges",
        "cartoon-ribbon-mesh",
        "nucleotide-direction-visuals",
        "obj-export-parity",
        "stl-export-parity-diagnostics",
        "ply-export-contract",
        "reference-converter-tooling",
    ] {
        assert!(
            by_domain.contains_key(domain),
            "missing parity evidence domain: {domain}"
        );
    }

    for domain in [
        "cartoon-ribbon-mesh",
        "nucleotide-direction-visuals",
        "carbohydrate-branched-visuals",
        "ball-and-stick-component-visuals",
        "obj-export-parity",
        "stl-export-parity-diagnostics",
        "ply-export-contract",
    ] {
        let entry = by_domain
            .get(domain)
            .unwrap_or_else(|| panic!("missing mesh evidence domain: {domain}"));
        assert!(
            !entry.reference_comparisons.is_empty(),
            "mesh evidence domain must list a reference comparison: {domain}"
        );
    }
}

#[test]
fn molstar_parity_checked_items_are_mapped_to_evidence_domains() {
    let Some(checklist) = read_internal_doc("molstar-parity-checklist.md") else {
        eprintln!("skipping internal Mol* parity checklist evidence audit; dev-docs is absent");
        return;
    };
    let entries = parity_evidence_entries();
    let by_domain = entries
        .iter()
        .map(|entry| (entry.domain.as_str(), entry))
        .collect::<std::collections::BTreeMap<_, _>>();
    let section_domains = parity_section_evidence_domains();
    let mut observed_section_checked_items = std::collections::BTreeMap::<String, usize>::new();
    let mesh_sections = [
        "Representation Object Layer",
        "Mesh Builder Exactness",
        "Cartoon And Ribbon Geometry",
        "Nucleotide Geometry",
        "Surface And Volume Representations",
        "Exporter Exactness",
        "Molfig Static PLY Contract",
    ];
    let mut current_section = "";
    let mut checked_items = 0usize;

    for line in checklist.lines() {
        if let Some(section) = line.strip_prefix("## ") {
            current_section = section.trim();
            continue;
        }
        if !line.trim_start().starts_with("- [x]") {
            continue;
        }
        checked_items += 1;
        *observed_section_checked_items
            .entry(current_section.to_string())
            .or_default() += 1;
        let section_evidence = section_domains
            .get(current_section)
            .unwrap_or_else(|| panic!("missing evidence mapping for section: {current_section}"));
        let domains = &section_evidence.domains;
        assert!(
            !domains.is_empty(),
            "empty evidence domain mapping for section: {current_section}"
        );
        for domain in domains {
            assert!(
                by_domain.contains_key(domain.as_str()),
                "section {current_section} maps to missing evidence domain: {domain}"
            );
        }
        assert!(
            domains
                .iter()
                .filter_map(|domain| by_domain.get(domain.as_str()))
                .any(|entry| !entry.unit_tests.is_empty()),
            "checked item lacks unit-test evidence via section {current_section}: {line}"
        );
        assert!(
            domains
                .iter()
                .filter_map(|domain| by_domain.get(domain.as_str()))
                .any(|entry| !entry.fixtures.is_empty()),
            "checked item lacks fixture evidence via section {current_section}: {line}"
        );
        if mesh_sections.contains(&current_section) {
            assert!(
                domains
                    .iter()
                    .filter_map(|domain| by_domain.get(domain.as_str()))
                    .any(|entry| !entry.reference_comparisons.is_empty()),
                "checked mesh item lacks reference-comparison evidence via section {current_section}: {line}"
            );
        }
    }

    for (section, evidence) in &section_domains {
        let observed = observed_section_checked_items
            .get(section)
            .copied()
            .unwrap_or(0);
        assert_eq!(
            observed, evidence.checked_items,
            "checked-item count drifted for section {section}; update tests/expected/molstar-parity-section-evidence.tsv"
        );
    }
    assert!(
        checked_items >= 400,
        "checklist evidence audit should scan the full checked item set"
    );
}

#[test]
fn molstar_parity_checked_mesh_items_have_molstar_reference_export_evidence() {
    let Some(checklist) = read_internal_doc("molstar-parity-checklist.md") else {
        eprintln!("skipping internal Mol* parity checklist mesh audit; dev-docs is absent");
        return;
    };
    let entries = parity_evidence_entries();
    let by_domain = entries
        .iter()
        .map(|entry| (entry.domain.as_str(), entry))
        .collect::<std::collections::BTreeMap<_, _>>();
    let mesh_references = parity_mesh_reference_evidence();
    let section_checked_items = checklist_checked_items_by_section(&checklist);

    assert_eq!(
        mesh_references
            .keys()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec![
            "Cartoon And Ribbon Geometry",
            "Exporter Exactness",
            "Mesh Builder Exactness",
            "Nucleotide Geometry",
            "Representation Object Layer",
            "Surface And Volume Representations",
        ],
        "mesh reference evidence should enumerate the Mol* mesh/export sections explicitly"
    );

    for (section, evidence) in &mesh_references {
        let observed = section_checked_items.get(section).copied().unwrap_or(0);
        assert_eq!(
            observed, evidence.checked_items,
            "checked mesh item count drifted for {section}; update tests/expected/molstar-parity-mesh-reference-evidence.tsv"
        );
        assert!(
            !evidence.notes.trim().is_empty(),
            "missing mesh reference evidence notes for {section}"
        );
        for domain in &evidence.domains {
            let entry = by_domain.get(domain.as_str()).unwrap_or_else(|| {
                panic!("mesh reference section {section} maps to missing domain {domain}")
            });
            assert!(
                !entry.reference_comparisons.is_empty(),
                "mesh reference domain {domain} for {section} has no reference comparisons"
            );
        }
        assert!(
            evidence.domains.iter().any(|domain| {
                by_domain
                    .get(domain.as_str())
                    .map(|entry| {
                        entry.reference_comparisons.iter().any(|reference| {
                            reference.contains("molstar")
                                || reference.contains("9R1O.obj")
                                || reference.contains("9R1O.stl")
                        })
                    })
                    .unwrap_or(false)
            }),
            "mesh reference section {section} must include an explicit Mol* reference artifact domain"
        );
    }
}

#[test]
fn molstar_reference_converter_supports_external_runtime_module_dirs() {
    let script = include_str!("../../scripts/molstar-reference-convert.mjs");
    for snippet in [
        "--runtime-module-dir <path>",
        "--render-object-report",
        "--scene-source <source>",
        "--force-cylinder-impostors",
        "runtimeModuleDirs",
        "runtimeDependencyResolvers",
        "validateSceneSource",
        "configureSceneRuntime",
        "GL_EXT_frag_depth",
        "forceCylinderImpostorSupport",
        "restoreCylinderImpostorSupport",
        "structureVisualCommon",
        "exportableDrawCount",
        "mol-plugin-state/actions/file.js",
        "mol-repr/structure/visual/util/common.js",
        "mol-util/zip/zip.js",
        "installHeadlessDomShim",
        "patchHeadlessCanvas3D",
        "installMolstarVersionShim",
        "configureExporter",
        "export-primitives-quality",
        "scene-source",
        "buildSceneFromDataFormat",
        "buildSceneFromOpenFilesAction",
        "openFilesFormatParam",
        "molstar.OpenFiles",
        "dataFormatProviderId",
        "hierarchyPresetParams",
        "representation-preset",
        "representationPresetParams",
        "TrajectoryFromPDB);",
        "displayPath",
        "inspectToolchain",
        "Runtime toolchain:",
        "pass --runtime-module-dir",
        "validateObjStlSparseSlotDriftReference",
        "objStlSparseSlotDrift",
        "obj_stl_sparse_slot_rounding_mismatch_count",
    ] {
        assert!(
            script.contains(snippet),
            "reference converter should support external runtime module dirs: {snippet}"
        );
    }
}

#[test]
fn molstar_browser_reference_converter_uses_real_chrome_webgl_export_path() {
    let script = include_str!("../../scripts/molstar-browser-reference-convert.mjs");
    let harness = include_str!("../../scripts/molstar-browser-reference-harness.ts");

    for snippet in [
        "molstar-browser-reference-harness.ts",
        "--chrome <path>",
        "--render-object-report",
        "--compare-references",
        "--compare-molfig",
        "--molfig-diff <path|cargo>",
        "--debug-stl-facet <n[,n]>",
        "Chrome DevTools",
        "remote-debugging-port",
        "GL_EXT_frag_depth",
        "headless Chrome",
        "exportableDrawCount",
        "/upload/",
        "compareBrowserOutputsToReferences",
        "compareBrowserOutputsToMolfig",
        "compareBrowserStlFacetDebugToMolfig",
        "--stl-export-facet-context",
        "center_offset_delta",
        "box_min_delta",
        "box_max_delta",
        "box_min_indices_match",
        "box_max_indices_match",
        "box_min_point_deltas",
        "box_max_point_deltas",
        "sphereCenter=",
        "sphereRadius=",
        "sphereExtrema=",
        "molfig-vs-browser",
        "printBrowserStlFacetDebug",
        "stl_context=facet",
    ] {
        assert!(
            script.contains(snippet),
            "browser reference converter should expose real Chrome/WebGL export support: {snippet}"
        );
    }

    for snippet in [
        "window.molfigBrowserReferenceExport",
        "OpenFiles",
        "ObjExporter",
        "StlExporter",
        "getRenderObjects",
        "canvas3d!.webgl",
        "fragDepth",
        "debugStlFacet",
        "centeredVertexBits",
        "boxMin",
        "boxMax",
        "boxMinIndices",
        "boxMaxIndices",
        "boundingSphere",
        "sphereSummary",
        "extremaCount",
        "exporter.add(renderObjects[i], plugin.canvas3d!.webgl, ctx)",
    ] {
        assert!(
            harness.contains(snippet),
            "browser reference harness should run Mol* exporters against browser render objects: {snippet}"
        );
    }
}

#[test]
fn molstar_reference_converter_accepts_format_scoped_artifact_metadata() {
    let repo = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = match std::process::Command::new("node")
        .current_dir(repo)
        .arg("scripts/molstar-reference-convert.mjs")
        .arg("--dry-run")
        .arg("--no-build-from-source")
        .arg("--manifest")
        .arg("tests/expected/molstar-reference/format-scoped-reference-fixtures.txt")
        .arg("--artifact-manifest")
        .arg("tests/expected/molstar-reference/format-scoped-reference-artifacts.json")
        .output()
    {
        Ok(output) => output,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return,
        Err(error) => panic!("failed to run node: {error}"),
    };

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "format-scoped artifact dry-run failed\nstdout:\n{stdout}\nstderr:\n{stderr}"
    );
    assert!(stdout.contains("targets=json"));
    assert!(stdout.contains("- PASS 9r1o-molstar-assembly-1: json;"));
    assert!(stdout.contains("Existing Mol* reference artifact validation:"));
}

#[test]
fn molstar_reference_converter_accepts_scene_source_override() {
    let repo = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = match std::process::Command::new("node")
        .current_dir(repo)
        .arg("scripts/molstar-reference-convert.mjs")
        .arg("--dry-run")
        .arg("--no-build-from-source")
        .arg("--scene-source")
        .arg("open-files")
        .output()
    {
        Ok(output) => output,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return,
        Err(error) => panic!("failed to run node: {error}"),
    };

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "scene-source override dry-run failed\nstdout:\n{stdout}\nstderr:\n{stderr}"
    );
    assert!(stdout.contains("scene=open-files"));
    assert!(stdout.contains("Existing Mol* reference artifact validation:"));
    assert!(stdout.contains("- PASS 9r1o-molstar-assembly-1: json/obj/stl;"));
}

#[test]
fn pinned_molstar_geo_export_has_no_ply_exporter() {
    let Some(controls) = read_molstar_source("extensions/geo-export/controls.ts") else {
        eprintln!("skipping pinned Mol* geo-export source audit; artifacts is absent");
        return;
    };

    for format in ["glb", "stl", "obj", "usdz"] {
        assert!(controls.contains(&format!("['{format}',")));
        assert!(controls.contains(&format!("case '{format}'")));
    }
    assert!(!controls.contains("['ply',"));
    assert!(!controls.contains("case 'ply'"));
    assert!(!controls.contains("PlyExporter"));
    assert!(!controls.contains("./ply-exporter"));
    assert!(!std::path::Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/artifacts/molstar/src/extensions/geo-export/ply-exporter.ts"
    ))
    .exists());
}

#[test]
fn pinned_molstar_obj_face_loop_uses_draw_count_index_order() {
    let Some(exporter) = read_molstar_source("extensions/geo-export/obj-exporter.ts") else {
        eprintln!("skipping pinned Mol* OBJ exporter source audit; artifacts is absent");
        return;
    };

    assert!(exporter.contains("for (let i = 0; i < drawCount; i += 3)"));
    assert!(
        exporter.contains("const v1 = this.vertexOffset + (isGeoTexture ? i : indices![i]) + 1;")
    );
    assert!(exporter
        .contains("const v2 = this.vertexOffset + (isGeoTexture ? i + 1 : indices![i + 1]) + 1;"));
    assert!(exporter
        .contains("const v3 = this.vertexOffset + (isGeoTexture ? i + 2 : indices![i + 2]) + 1;"));
    assert!(exporter.contains("StringBuilder.writeSafe(obj, 'f ');"));
}

#[test]
fn parses_entity_poly_sequence_and_struct_asym_from_cif() {
    let cif = b"data_demo\nloop_\n_entry.id\nDEMO\n#\nloop_\n_exptl.method\n'X-RAY DIFFRACTION'\n#\nloop_\n_entity.id\n_entity.type\n_entity.pdbx_description\n1 polymer 'test protein'\n#\nloop_\n_entity_poly.entity_id\n_entity_poly.type\n_entity_poly.pdbx_seq_one_letter_code\n_entity_poly.nstd_linkage\n_entity_poly.nstd_monomer\n1 'polypeptide(L)' AG no no\n#\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n_entity_poly_seq.hetero\n1 1 ALA n\n1 2 GLY n\n#\nloop_\n_struct_asym.id\n_struct_asym.entity_id\n_struct_asym.details\nA 1 'chain A'\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.B_iso_or_equiv\n_atom_site.pdbx_formal_charge\nATOM 1 C CA ALA A 1 0.000 0.000 0.000 12.5 1\nATOM 2 C CA GLY A 2 1.000 0.000 0.000 13.5 0\n#\nloop_\n_atom_site_anisotrop.id\n_atom_site_anisotrop.U[1][1]\n_atom_site_anisotrop.U[1][2]\n_atom_site_anisotrop.U[1][3]\n_atom_site_anisotrop.U[2][1]\n_atom_site_anisotrop.U[2][2]\n_atom_site_anisotrop.U[2][3]\n_atom_site_anisotrop.U[3][1]\n_atom_site_anisotrop.U[3][2]\n_atom_site_anisotrop.U[3][3]\n1 0.10 0.01 0.02 0.01 0.11 0.03 0.02 0.03 0.12\n#\n";
    let mol = parse_molecule(cif, InputFormat::Cif).unwrap();

    assert_eq!(mol.entries.len(), 1);
    assert_eq!(mol.entries[0].id, "DEMO");
    assert_eq!(mol.experiments.len(), 1);
    assert_eq!(mol.experiments[0].method, "X-RAY DIFFRACTION");
    assert_eq!(mol.entities.len(), 1);
    assert_eq!(mol.entities[0].id, "1");
    assert_eq!(mol.entities[0].type_name, "polymer");
    assert_eq!(mol.entity_polymers.len(), 1);
    assert_eq!(mol.entity_polymers[0].entity_id, "1");
    assert_eq!(mol.entity_polymers[0].sequence, "AG");
    assert_eq!(mol.entity_poly_seq.len(), 2);
    assert_eq!(mol.entity_poly_seq[1].mon_id, "GLY");
    assert_eq!(mol.struct_asym.len(), 1);
    assert_eq!(mol.atoms[0].entity_id, "1");
    assert_eq!(mol.atoms[0].b_iso, 12.5);
    assert_eq!(mol.atoms[0].formal_charge, 1);
    assert_eq!(mol.atom_site_anisotrop.len(), 1);
    assert_eq!(mol.atom_site_anisotrop[0].atom_id, 1);
    assert_eq!(mol.atom_site_anisotrop[0].u[2][2], 0.12);
    let structure = mol.atomic_structure();
    assert_eq!(structure.model.conformation.atom_ids, vec![1, 2]);
    assert_eq!(structure.model.conformation.occupancies, vec![1.0, 1.0]);
    assert_eq!(structure.model.conformation.b_iso, vec![12.5, 13.5]);
    assert_eq!(structure.model.conformation.formal_charges, vec![1, 0]);
    assert!(structure.model.conformation.occupancy_defined);
    assert!(structure.model.conformation.b_iso_defined);
    assert!(structure.model.conformation.xyz_defined);
    assert_eq!(structure.properties.pdbx_formal_charge, vec![1, 0]);
    assert_eq!(structure.model.hierarchy.atoms[0].formal_charge, 1);
    assert_eq!(
        structure.model.conformation.element_to_anisotrop,
        vec![0, -1]
    );
    assert_eq!(
        structure.model.conformation.anisotropic_displacement,
        vec![
            Some([[0.10, 0.01, 0.02], [0.01, 0.11, 0.03], [0.02, 0.03, 0.12]]),
            None
        ]
    );
    assert_eq!(structure.model.hierarchy.chains[0].entity_id, "1");
    assert_eq!(structure.model.sequence.sequences.len(), 1);
    assert_eq!(structure.model.sequence.by_entity_key.get(&0), Some(&0));
    assert_eq!(
        structure.model.sequence.sequences[0]
            .residues
            .iter()
            .map(|residue| (residue.comp_id.as_str(), residue.seq_id))
            .collect::<Vec<_>>(),
        vec![("ALA", 1), ("GLY", 2)]
    );

    let info = String::from_utf8(molecule_info(cif, br#"{"format":"cif"}"#).unwrap()).unwrap();
    assert!(info.contains(r#""entry_count":1"#));
    assert!(info.contains(r#""experiment_count":1"#));
    assert!(info.contains(r#""anisotrop_count":1"#));
    assert!(info.contains(r#""entity_count":1"#));
    assert!(info.contains(r#""entity_poly_count":1"#));
    assert!(info.contains(r#""entity_poly_seq_count":2"#));
    assert!(info.contains(r#""struct_asym_count":1"#));
    assert!(info.contains(r#""atom_id_count":2"#));
    assert!(info.contains(r#""b_iso_count":2"#));
    assert!(info.contains(r#""formal_charge_count":2"#));
    assert!(info.contains(r#""occupancy_defined":true"#));
    assert!(info.contains(r#""b_iso_defined":true"#));
    assert!(info.contains(r#""xyz_defined":true"#));
    assert!(info.contains(r#""element_to_anisotrop_count":2"#));
    assert!(info.contains(r#""anisotropic_displacement_count":1"#));
}

#[test]
fn parses_pdbx_branch_and_sequence_scheme_categories_from_cif() {
    let cif = b"data_demo\nloop_\n_entity.id\n_entity.type\n1 polymer\n2 branched\n3 non-polymer\n#\nloop_\n_pdbx_entity_branch.entity_id\n_pdbx_entity_branch.type\n2 oligosaccharide\n#\nloop_\n_pdbx_entity_branch_link.link_id\n_pdbx_entity_branch_link.details\n_pdbx_entity_branch_link.entity_id\n_pdbx_entity_branch_link.entity_branch_list_num_1\n_pdbx_entity_branch_link.entity_branch_list_num_2\n_pdbx_entity_branch_link.comp_id_1\n_pdbx_entity_branch_link.comp_id_2\n_pdbx_entity_branch_link.atom_id_1\n_pdbx_entity_branch_link.leaving_atom_id_1\n_pdbx_entity_branch_link.atom_stereo_config_1\n_pdbx_entity_branch_link.atom_id_2\n_pdbx_entity_branch_link.leaving_atom_id_2\n_pdbx_entity_branch_link.atom_stereo_config_2\n_pdbx_entity_branch_link.value_order\n1 'test glycosidic link' 2 1 2 NAG MAN C1 O1 n O4 HO4 n sing\n#\nloop_\n_pdbx_branch_scheme.entity_id\n_pdbx_branch_scheme.hetero\n_pdbx_branch_scheme.asym_id\n_pdbx_branch_scheme.mon_id\n_pdbx_branch_scheme.num\n_pdbx_branch_scheme.pdb_asym_id\n_pdbx_branch_scheme.pdb_seq_num\n_pdbx_branch_scheme.pdb_mon_id\n_pdbx_branch_scheme.auth_asym_id\n_pdbx_branch_scheme.auth_seq_num\n_pdbx_branch_scheme.auth_mon_id\n2 n B NAG 1 B 101 NAG BA 501 NAG\n#\nloop_\n_pdbx_nonpoly_scheme.asym_id\n_pdbx_nonpoly_scheme.entity_id\n_pdbx_nonpoly_scheme.mon_id\n_pdbx_nonpoly_scheme.pdb_strand_id\n_pdbx_nonpoly_scheme.ndb_seq_num\n_pdbx_nonpoly_scheme.pdb_seq_num\n_pdbx_nonpoly_scheme.auth_seq_num\n_pdbx_nonpoly_scheme.pdb_mon_id\n_pdbx_nonpoly_scheme.auth_mon_id\n_pdbx_nonpoly_scheme.pdb_ins_code\nL 3 HEM L 1 701 9001 HEM HEM .\n#\nloop_\n_pdbx_poly_seq_scheme.asym_id\n_pdbx_poly_seq_scheme.entity_id\n_pdbx_poly_seq_scheme.seq_id\n_pdbx_poly_seq_scheme.mon_id\n_pdbx_poly_seq_scheme.ndb_seq_num\n_pdbx_poly_seq_scheme.pdb_seq_num\n_pdbx_poly_seq_scheme.auth_seq_num\n_pdbx_poly_seq_scheme.pdb_mon_id\n_pdbx_poly_seq_scheme.auth_mon_id\n_pdbx_poly_seq_scheme.pdb_strand_id\n_pdbx_poly_seq_scheme.pdb_ins_code\n_pdbx_poly_seq_scheme.hetero\nA 1 1 ALA 1 1 10 ALA ALA A . n\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A 1 1 0.000 0.000 0.000\n#\n";
    let mol = parse_molecule(cif, InputFormat::Cif).unwrap();

    assert_eq!(mol.pdbx_entity_branch.len(), 1);
    assert_eq!(mol.pdbx_entity_branch[0].entity_id, "2");
    assert_eq!(mol.pdbx_entity_branch[0].type_name, "oligosaccharide");
    assert_eq!(mol.pdbx_entity_branch_links.len(), 1);
    assert_eq!(mol.pdbx_entity_branch_links[0].link_id, 1);
    assert_eq!(
        mol.pdbx_entity_branch_links[0].details,
        "test glycosidic link"
    );
    assert_eq!(mol.pdbx_entity_branch_links[0].entity_branch_list_num_2, 2);
    assert_eq!(mol.pdbx_entity_branch_links[0].atom_id_2, "O4");
    assert_eq!(mol.pdbx_entity_branch_links[0].value_order, "sing");
    assert_eq!(mol.pdbx_branch_scheme.len(), 1);
    assert_eq!(mol.pdbx_branch_scheme[0].asym_id, "B");
    assert_eq!(mol.pdbx_branch_scheme[0].num, 1);
    assert_eq!(mol.pdbx_branch_scheme[0].auth_seq_num, "501");
    assert_eq!(mol.pdbx_nonpoly_scheme.len(), 1);
    assert_eq!(mol.pdbx_nonpoly_scheme[0].pdb_ins_code, "");
    assert_eq!(mol.pdbx_poly_seq_scheme.len(), 1);
    assert_eq!(mol.pdbx_poly_seq_scheme[0].seq_id, 1);
    assert_eq!(mol.pdbx_poly_seq_scheme[0].pdb_strand_id, "A");

    let info = String::from_utf8(molecule_info(cif, br#"{"format":"cif"}"#).unwrap()).unwrap();
    assert!(info.contains(r#""pdbx_entity_branch_count":1"#));
    assert!(info.contains(r#""pdbx_entity_branch_link_count":1"#));
    assert!(info.contains(r#""pdbx_branch_scheme_count":1"#));
    assert!(info.contains(r#""pdbx_nonpoly_scheme_count":1"#));
    assert!(info.contains(r#""pdbx_poly_seq_scheme_count":1"#));
}

#[test]
fn branched_entity_sequence_and_link_maps_follow_entity_branch_numbers() {
    let mol = Molecule {
        pdbx_branch_scheme: vec![
            PdbxBranchScheme {
                entity_id: "2".into(),
                asym_id: "B".into(),
                mon_id: "NAG".into(),
                num: 1,
                pdb_asym_id: "B".into(),
                pdb_seq_num: "101".into(),
                pdb_mon_id: "NAG".into(),
                auth_asym_id: "BA".into(),
                auth_seq_num: "501".into(),
                auth_mon_id: "NAG".into(),
                ..PdbxBranchScheme::default()
            },
            PdbxBranchScheme {
                entity_id: "2".into(),
                asym_id: "B".into(),
                mon_id: "MAN".into(),
                num: 2,
                pdb_asym_id: "B".into(),
                pdb_seq_num: "102".into(),
                pdb_mon_id: "MAN".into(),
                auth_asym_id: "BA".into(),
                auth_seq_num: "502".into(),
                auth_mon_id: "MAN".into(),
                ..PdbxBranchScheme::default()
            },
            PdbxBranchScheme {
                entity_id: "2".into(),
                asym_id: "C".into(),
                mon_id: "NAG".into(),
                num: 1,
                pdb_asym_id: "C".into(),
                pdb_seq_num: "201".into(),
                pdb_mon_id: "NAG".into(),
                auth_asym_id: "CA".into(),
                auth_seq_num: "601".into(),
                auth_mon_id: "NAG".into(),
                ..PdbxBranchScheme::default()
            },
        ],
        pdbx_entity_branch_links: vec![PdbxEntityBranchLink {
            link_id: 7,
            details: "test glycosidic link".into(),
            entity_id: "2".into(),
            entity_branch_list_num_1: 1,
            entity_branch_list_num_2: 2,
            comp_id_1: "NAG".into(),
            comp_id_2: "MAN".into(),
            atom_id_1: "C1".into(),
            atom_id_2: "O4".into(),
            value_order: "sing".into(),
            ..PdbxEntityBranchLink::default()
        }],
        ..Molecule::default()
    };

    let sequence_map = mol.branched_sequence_map();
    assert_eq!(sequence_map.entries_for_entity("2"), &[0, 1, 2]);
    assert_eq!(sequence_map.entries_for_asym("B"), &[0, 1]);
    assert_eq!(
        sequence_map
            .entry_for_asym_num("B", 2)
            .map(|entry| entry.auth_seq_num.as_str()),
        Some("502")
    );
    assert_eq!(sequence_map.entries_for_entity_num("2", 1), &[0, 2]);

    let link_map = mol.branched_entity_link_map();
    assert_eq!(link_map.links_for_entity("2"), &[0]);
    assert_eq!(link_map.links_for_entity_num("2", 1), &[0]);
    assert_eq!(link_map.links_for_entity_num("2", 2), &[0]);
    assert_eq!(
        link_map
            .link_for_entity_link_id("2", 7)
            .map(|link| link.atom_id_2.as_str()),
        Some("O4")
    );
    assert_eq!(
        link_map.placements,
        vec![BranchedEntityLinkPlacement {
            link_index: 0,
            entry_index_1: 0,
            entry_index_2: 1,
        }]
    );
}

#[test]
fn carbohydrate_branched_entity_fixture_parses_and_derives_carbohydrates() {
    let cif = include_bytes!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/cif/carbohydrate-branched.cif"
    ));
    let mol = parse_molecule(cif, InputFormat::Cif).unwrap();

    assert_eq!(mol.entities.len(), 1);
    assert_eq!(mol.entities[0].type_name, "branched");
    assert_eq!(mol.entity_index.subtype, vec!["oligosaccharide"]);
    assert_eq!(mol.pdbx_entity_branch.len(), 1);
    assert_eq!(mol.pdbx_entity_branch[0].entity_id, "1");
    assert_eq!(mol.pdbx_entity_branch[0].type_name, "oligosaccharide");
    assert_eq!(mol.pdbx_entity_branch_links.len(), 1);
    let branch_link = &mol.pdbx_entity_branch_links[0];
    assert_eq!(branch_link.link_id, 1);
    assert_eq!(branch_link.entity_id, "1");
    assert_eq!(branch_link.entity_branch_list_num_1, 1);
    assert_eq!(branch_link.entity_branch_list_num_2, 2);
    assert_eq!(branch_link.comp_id_1, "NAG");
    assert_eq!(branch_link.comp_id_2, "MAN");
    assert_eq!(branch_link.atom_id_1, "C1");
    assert_eq!(branch_link.atom_id_2, "O4");
    assert_eq!(branch_link.value_order, "sing");
    assert_eq!(mol.pdbx_branch_scheme.len(), 2);
    assert_eq!(
        mol.pdbx_branch_scheme
            .iter()
            .map(|entry| (
                entry.mon_id.as_str(),
                entry.num,
                entry.auth_seq_num.as_str()
            ))
            .collect::<Vec<_>>(),
        vec![("NAG", 1, "501"), ("MAN", 2, "502")]
    );

    let sequence_map = mol.branched_sequence_map();
    assert_eq!(sequence_map.entries_for_entity("1"), &[0, 1]);
    assert_eq!(
        sequence_map
            .entry_for_asym_num("B", 2)
            .map(|entry| entry.mon_id.as_str()),
        Some("MAN")
    );
    let link_map = mol.branched_entity_link_map();
    assert_eq!(link_map.links_for_entity_num("1", 1), &[0]);
    assert_eq!(link_map.links_for_entity_num("1", 2), &[0]);
    assert_eq!(
        link_map.placements,
        vec![BranchedEntityLinkPlacement {
            link_index: 0,
            entry_index_1: 0,
            entry_index_2: 1,
        }]
    );

    assert_eq!(mol.atoms.len(), 14);
    assert_eq!(mol.chemical_components.len(), 2);
    assert!(saccharide_component("NAG").is_some());
    assert!(saccharide_component("MAN").is_some());
    let carbohydrates = mol.carbohydrates();
    assert_eq!(carbohydrates.elements.len(), 2);
    assert!(carbohydrates.partial_elements.is_empty());
    assert_eq!(
        carbohydrates
            .elements
            .iter()
            .map(|element| element.component.component_type)
            .collect::<Vec<_>>(),
        vec![SaccharideType::HexNAc, SaccharideType::Hexose]
    );
    assert_eq!(
        carbohydrates
            .links
            .iter()
            .map(|link| (link.carbohydrate_index_a, link.carbohydrate_index_b))
            .collect::<Vec<_>>(),
        vec![(0, 1), (1, 0)]
    );
    let structure = mol.atomic_structure();
    assert_eq!(structure.carbohydrate_element_indices(0, 0), &[0]);
    assert_eq!(structure.carbohydrate_element_indices(0, 7), &[1]);
    assert_eq!(structure.carbohydrate_link_indices(0, 0), &[0]);
    assert_eq!(structure.carbohydrate_link_indices(0, 7), &[1]);

    let info = String::from_utf8(molecule_info(cif, br#"{"format":"cif"}"#).unwrap()).unwrap();
    assert!(info.contains(r#""pdbx_entity_branch_count":1"#));
    assert!(info.contains(r#""pdbx_entity_branch_link_count":1"#));
    assert!(info.contains(r#""pdbx_branch_scheme_count":2"#));
    assert!(info.contains(r#""chem_comp_count":2"#));
    assert!(info.contains(r#""chem_comp_bond_count":14"#));
}

#[test]
fn carbohydrate_branched_fixture_matches_molstar_reference_summary() {
    let cif = include_bytes!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/cif/carbohydrate-branched.cif"
    ));
    let mol = parse_molecule(cif, InputFormat::Cif).unwrap();
    let actual = carbohydrate_reference_summary_json(&mol);
    let expected =
        include_str!("../../tests/expected/carbohydrate-reference-summary.json").trim_end();

    assert_eq!(actual, expected, "\n{actual}");
}

#[test]
fn ligand_metal_aromatic_fixture_covers_dictionary_rings_and_struct_conn() {
    let cif = include_bytes!("../../tests/fixtures/cif/ligand-metal-aromatic.cif");
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.atoms.len(), 12);
    assert_eq!(mol.chemical_components.len(), 2);
    assert_eq!(mol.chemical_components[0].id, "LMR");
    assert_eq!(mol.chemical_components[0].formula, "C10 O1");
    assert_eq!(mol.chemical_components[0].pdbx_formal_charge, Some(-1));
    assert_eq!(mol.chemical_component_atoms.len(), 11);
    assert_eq!(mol.chemical_component_bonds.len(), 12);

    let oxygen = mol
        .chemical_component_atoms
        .iter()
        .find(|atom| atom.atom_id == "O1")
        .unwrap();
    assert_eq!(oxygen.charge, -1);
    assert!(oxygen.leaving_atom);
    assert_eq!(oxygen.stereo_config, "R");
    assert_eq!(oxygen.model_cartn, Some(vec3(-1.2, -0.8, 0.0)));

    let substituent = mol
        .chemical_component_bonds
        .iter()
        .find(|bond| bond.atom_id_1 == "C1" && bond.atom_id_2 == "O1")
        .unwrap();
    assert_eq!(substituent.order, 1);
    assert_eq!(substituent.ordinal, Some(112));
    assert_eq!(substituent.stereo_config, "R");
    assert!(!substituent.flags.contains(BondFlags::AROMATIC));

    assert_eq!(mol.bonds.len(), 13);
    assert_eq!(mol.resonance.ring_count, 2);
    assert_eq!(mol.resonance.aromatic_ring_count, 2);
    assert_eq!(mol.resonance.delocalized_bond_count, 11);
    assert_eq!(mol.rings[0].atom_indices, vec![0, 1, 2, 3, 4, 5]);
    assert_eq!(mol.rings[1].atom_indices, vec![3, 4, 6, 7, 8, 9]);
    assert!(mol.rings.iter().all(|ring| ring.aromatic));

    let struct_conn_index = mol
        .bond_metadata
        .iter()
        .position(|metadata| metadata.source == BondSource::StructConn)
        .unwrap();
    let metadata = &mol.bond_metadata[struct_conn_index];
    assert_eq!(metadata.order, 1);
    assert!(metadata.flags.contains(BondFlags::METALLIC_COORDINATION));
    assert!(!metadata.flags.contains(BondFlags::COVALENT));
    assert_eq!(metadata.distance, Some(2.10));
    let struct_conn = metadata.struct_conn.as_ref().unwrap();
    assert_eq!(struct_conn.id, "zn-o1");
    assert_eq!(struct_conn.partner_a_comp_id, "ZN");
    assert_eq!(struct_conn.partner_b_comp_id, "LMR");
    assert_eq!(mol.atoms[struct_conn.partner_a_atom_index].name, "ZN");
    assert_eq!(mol.atoms[struct_conn.partner_b_atom_index].name, "O1");

    let structure = mol.atomic_structure();
    assert_eq!(structure.units.len(), 2);
    assert_eq!(structure.intra_unit_bond_count, 12);
    assert_eq!(structure.inter_unit_bonds.len(), 1);
    assert!(structure.inter_unit_bonds[0]
        .flags
        .contains(BondFlags::METALLIC_COORDINATION));

    let info =
        String::from_utf8(molecule_info(cif, br#"{"format":"cif","infer-bonds":false}"#).unwrap())
            .unwrap();
    assert!(info.contains(r#""chem_comp_atom_count":11"#));
    assert!(info.contains(r#""chem_comp_bond_count":12"#));
    assert!(info.contains(r#""struct_conn":1"#));
    assert!(info.contains(r#""aromatic_rings":2"#));
}

#[test]
fn parses_ihm_coarse_sites_and_exports_semantic_meshes() {
    let cif = b"data_demo\nloop_\n_ihm_model_list.model_id\n_ihm_model_list.model_name\n_ihm_model_list.assembly_id\n_ihm_model_list.protocol_id\n_ihm_model_list.representation_id\n1 'model one' 1 1 1\n#\nloop_\n_ihm_model_group.id\n_ihm_model_group.name\n_ihm_model_group.details\n1 ensemble 'test group'\n#\nloop_\n_ihm_model_group_link.model_id\n_ihm_model_group_link.group_id\n1 1\n#\nloop_\n_ihm_sphere_obj_site.id\n_ihm_sphere_obj_site.entity_id\n_ihm_sphere_obj_site.asym_id\n_ihm_sphere_obj_site.seq_id_begin\n_ihm_sphere_obj_site.seq_id_end\n_ihm_sphere_obj_site.Cartn_x\n_ihm_sphere_obj_site.Cartn_y\n_ihm_sphere_obj_site.Cartn_z\n_ihm_sphere_obj_site.object_radius\n1 1 A 1 10 0.0 0.0 0.0 2.0\n#\nloop_\n_ihm_gaussian_obj_site.id\n_ihm_gaussian_obj_site.entity_id\n_ihm_gaussian_obj_site.asym_id\n_ihm_gaussian_obj_site.seq_id_begin\n_ihm_gaussian_obj_site.seq_id_end\n_ihm_gaussian_obj_site.mean_Cartn_x\n_ihm_gaussian_obj_site.mean_Cartn_y\n_ihm_gaussian_obj_site.mean_Cartn_z\n_ihm_gaussian_obj_site.weight\n_ihm_gaussian_obj_site.covariance_matrix[1][1]\n_ihm_gaussian_obj_site.covariance_matrix[1][2]\n_ihm_gaussian_obj_site.covariance_matrix[1][3]\n_ihm_gaussian_obj_site.covariance_matrix[2][1]\n_ihm_gaussian_obj_site.covariance_matrix[2][2]\n_ihm_gaussian_obj_site.covariance_matrix[2][3]\n_ihm_gaussian_obj_site.covariance_matrix[3][1]\n_ihm_gaussian_obj_site.covariance_matrix[3][2]\n_ihm_gaussian_obj_site.covariance_matrix[3][3]\n1 1 A 11 20 4.0 0.0 0.0 1.0 4.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0\n#\nloop_\n_ihm_cross_link_restraint.id\n_ihm_cross_link_restraint.group_id\n_ihm_cross_link_restraint.entity_id_1\n_ihm_cross_link_restraint.entity_id_2\n_ihm_cross_link_restraint.asym_id_1\n_ihm_cross_link_restraint.asym_id_2\n_ihm_cross_link_restraint.comp_id_1\n_ihm_cross_link_restraint.comp_id_2\n_ihm_cross_link_restraint.seq_id_1\n_ihm_cross_link_restraint.seq_id_2\n_ihm_cross_link_restraint.atom_id_1\n_ihm_cross_link_restraint.atom_id_2\n_ihm_cross_link_restraint.restraint_type\n_ihm_cross_link_restraint.conditional_crosslink_flag\n_ihm_cross_link_restraint.model_granularity\n_ihm_cross_link_restraint.distance_threshold\n_ihm_cross_link_restraint.psi\n_ihm_cross_link_restraint.sigma_1\n_ihm_cross_link_restraint.sigma_2\n1 1 1 1 A A ALA GLY 1 10 CA CA 'upper bound' ALL by-residue 25.0 0.1 1.5 2.5\n#\n";
    let mol = parse_molecule(cif, InputFormat::Cif).unwrap();
    assert_eq!(mol.atoms.len(), 0);
    assert_eq!(mol.ihm_model_list.len(), 1);
    assert_eq!(mol.ihm_model_list[0].model_id, 1);
    assert_eq!(mol.ihm_model_list[0].model_name, "model one");
    assert_eq!(mol.ihm_model_groups.len(), 1);
    assert_eq!(mol.ihm_model_groups[0].details, "test group");
    assert_eq!(mol.ihm_model_group_links.len(), 1);
    assert_eq!(mol.ihm_model_group_links[0].group_id, 1);
    assert_eq!(mol.ihm_cross_link_restraints.len(), 1);
    assert_eq!(
        mol.ihm_cross_link_restraints[0].distance_threshold,
        Some(25.0)
    );
    assert_eq!(mol.coarse_spheres.len(), 1);
    assert_eq!(mol.coarse_gaussians.len(), 1);
    let structure = mol.atomic_structure();
    assert!(structure.coarse.hierarchy.is_defined);
    assert_eq!(structure.coarse.conformation.id.len(), 22);
    assert_eq!(
        structure.coarse.conformation.spheres,
        CoarseSphereConformation {
            x: vec![0.0],
            y: vec![0.0],
            z: vec![0.0],
            radius: vec![2.0],
            rmsf: vec![0.0],
        }
    );
    assert_eq!(
        structure.coarse.conformation.gaussians,
        CoarseGaussianConformation {
            x: vec![4.0],
            y: vec![0.0],
            z: vec![0.0],
            weight: vec![1.0],
            covariance_matrix: vec![[[4.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],]],
        }
    );
    assert_eq!(structure.coarse.hierarchy.spheres.polymer_ranges.len(), 1);
    assert_eq!(structure.coarse.hierarchy.gaussians.polymer_ranges.len(), 1);
    assert_eq!(structure.units.len(), 2);
    assert_eq!(structure.units[0].kind, UnitKind::Spheres);
    assert_eq!(structure.units[1].kind, UnitKind::Gaussians);
    assert_eq!(structure.units[0].invariant_id, 0);
    assert_eq!(structure.units[1].invariant_id, 1);
    assert_eq!(structure.units[0].traits, UnitTraits::NONE);
    assert_eq!(structure.units[1].traits, UnitTraits::NONE);
    assert_eq!(structure.units[0].props.polymer_elements, vec![0]);
    assert!(structure.units[0].props.gap_elements.is_empty());
    assert_eq!(structure.units[1].props.polymer_elements, vec![0]);
    assert!(structure.units[1].props.gap_elements.is_empty());
    assert!((structure.units[0].props.boundary.sphere.radius - 2.0).abs() < 0.000_1);
    assert_eq!(structure.units[1].props.boundary.sphere.radius, 0.0);
    assert_eq!(structure.position(0, 0).unwrap(), vec3(0.0, 0.0, 0.0));
    assert_eq!(structure.position(1, 0).unwrap(), vec3(4.0, 0.0, 0.0));
    assert!(structure.lookup3d.check(vec3(4.0, 0.0, 0.0), 0.1));

    let options = MeshOptions {
        format: InputFormat::Cif,
        center: false,
        ..MeshOptions::default()
    };
    let summary = render_object_summary_json(&mol, &options);
    assert!(summary.contains(r#""secondary_type":"coarse-sphere""#));
    assert!(summary.contains(r#""geometry_type":"ellipsoid""#));
    let info =
        String::from_utf8(molecule_info(cif, br#"{"format":"cif","center":false}"#).unwrap())
            .unwrap();
    assert!(info.contains(r#""ihm_model_count":1"#));
    assert!(info.contains(r#""ihm_model_group_count":1"#));
    assert!(info.contains(r#""ihm_model_group_link_count":1"#));
    assert!(info.contains(r#""ihm_cross_link_restraint_count":1"#));
    assert!(info.contains(r#""coarse_sphere_count":1"#));
    assert!(info.contains(r#""coarse_gaussian_count":1"#));
    assert!(info.contains(r#""unit_kind_counts":{"atomic":0,"spheres":1,"gaussians":1}"#));
    let obj =
        String::from_utf8(convert_to_obj(cif, br#"{"format":"cif","center":false}"#).unwrap())
            .unwrap();
    assert!(obj.contains("\nv "));
    assert!(obj.contains("\nf "));
}

#[test]
fn assembly_expanded_coarse_units_use_source_conformation_and_unit_operators() {
    let cif = b"data_demo\nloop_\n_ihm_sphere_obj_site.id\n_ihm_sphere_obj_site.entity_id\n_ihm_sphere_obj_site.asym_id\n_ihm_sphere_obj_site.seq_id_begin\n_ihm_sphere_obj_site.seq_id_end\n_ihm_sphere_obj_site.Cartn_x\n_ihm_sphere_obj_site.Cartn_y\n_ihm_sphere_obj_site.Cartn_z\n_ihm_sphere_obj_site.object_radius\n1 1 A 1 10 0.0 0.0 0.0 1.5\n2 1 B 1 10 100.0 0.0 0.0 1.0\n#\nloop_\n_ihm_gaussian_obj_site.id\n_ihm_gaussian_obj_site.entity_id\n_ihm_gaussian_obj_site.asym_id\n_ihm_gaussian_obj_site.seq_id_begin\n_ihm_gaussian_obj_site.seq_id_end\n_ihm_gaussian_obj_site.mean_Cartn_x\n_ihm_gaussian_obj_site.mean_Cartn_y\n_ihm_gaussian_obj_site.mean_Cartn_z\n_ihm_gaussian_obj_site.weight\n_ihm_gaussian_obj_site.covariance_matrix[1][1]\n_ihm_gaussian_obj_site.covariance_matrix[1][2]\n_ihm_gaussian_obj_site.covariance_matrix[1][3]\n_ihm_gaussian_obj_site.covariance_matrix[2][1]\n_ihm_gaussian_obj_site.covariance_matrix[2][2]\n_ihm_gaussian_obj_site.covariance_matrix[2][3]\n_ihm_gaussian_obj_site.covariance_matrix[3][1]\n_ihm_gaussian_obj_site.covariance_matrix[3][2]\n_ihm_gaussian_obj_site.covariance_matrix[3][3]\n1 1 A 11 20 4.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0\n#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n1 1,2 A\n#\nloop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\n1 1 0 0 0 0 1 0 0 0 0 1 0\n2 1 0 0 10 0 1 0 0 0 0 1 0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: Some("1".to_string()),
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = mol.atomic_structure();

    assert_eq!(structure.coarse.conformation.spheres.len(), 2);
    assert_eq!(structure.coarse.conformation.gaussians.len(), 1);
    assert_eq!(
        structure.coarse.conformation.spheres.position(0),
        Some(vec3(0.0, 0.0, 0.0))
    );
    assert_eq!(
        structure.coarse.conformation.spheres.position(1),
        Some(vec3(100.0, 0.0, 0.0))
    );
    assert_eq!(structure.units.len(), 4);
    assert_eq!(
        structure
            .units
            .iter()
            .map(|unit| (
                unit.kind,
                unit.elements.clone(),
                unit.operator.instance_id.clone()
            ))
            .collect::<Vec<_>>(),
        vec![
            (UnitKind::Spheres, vec![0], "ASM-1".to_string()),
            (UnitKind::Spheres, vec![0], "ASM-2".to_string()),
            (UnitKind::Gaussians, vec![0], "ASM-1".to_string()),
            (UnitKind::Gaussians, vec![0], "ASM-2".to_string()),
        ]
    );
    assert_eq!(
        structure
            .units
            .iter()
            .map(|unit| unit.invariant_id)
            .collect::<Vec<_>>(),
        vec![0, 0, 2, 2]
    );
    assert_eq!(structure.position(0, 0).unwrap(), vec3(0.0, 0.0, 0.0));
    assert_eq!(structure.position(1, 0).unwrap(), vec3(10.0, 0.0, 0.0));
    assert_eq!(structure.position(2, 0).unwrap(), vec3(4.0, 0.0, 0.0));
    assert_eq!(structure.position(3, 0).unwrap(), vec3(14.0, 0.0, 0.0));
    assert!(structure
        .units
        .iter()
        .all(|unit| unit.atom_indices == vec![0]));

    let geometry = mol.expanded_for_geometry();
    assert_eq!(
        geometry
            .coarse_spheres
            .iter()
            .map(|sphere| sphere.position.x)
            .collect::<Vec<_>>(),
        vec![0.0, 10.0]
    );
    assert_eq!(
        geometry
            .coarse_gaussians
            .iter()
            .map(|gaussian| gaussian.position.x)
            .collect::<Vec<_>>(),
        vec![4.0, 14.0]
    );
}

#[test]
fn assembly_coarse_units_keep_asymmetric_model_invariant_ids() {
    let cif = b"data_demo\n_entry.id DEMO\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 0.0 0.0 0.0\nATOM 2 C CA GLY B 1 1.0 0.0 0.0\n#\nloop_\n_ihm_sphere_obj_site.id\n_ihm_sphere_obj_site.entity_id\n_ihm_sphere_obj_site.asym_id\n_ihm_sphere_obj_site.seq_id_begin\n_ihm_sphere_obj_site.seq_id_end\n_ihm_sphere_obj_site.Cartn_x\n_ihm_sphere_obj_site.Cartn_y\n_ihm_sphere_obj_site.Cartn_z\n_ihm_sphere_obj_site.object_radius\n1 1 C 1 10 5.0 0.0 0.0 1.0\n#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n1 1 C\n#\nloop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\n1 1 0 0 0 0 1 0 0 0 0 1 0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: Some("1".to_string()),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = mol.atomic_structure();

    assert_eq!(structure.units.len(), 1);
    assert_eq!(structure.units[0].kind, UnitKind::Spheres);
    assert_eq!(structure.units[0].invariant_id, 2);
    assert_eq!(structure.units[0].chain_group_id, 2);
    assert_eq!(structure.symmetry_groups[0].invariant_id, 2);
}

#[test]
fn mixed_atomic_and_coarse_structure_keeps_all_unit_kinds_and_bounds() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 -2.0 0.0 0.0\n#\nloop_\n_ihm_sphere_obj_site.id\n_ihm_sphere_obj_site.entity_id\n_ihm_sphere_obj_site.asym_id\n_ihm_sphere_obj_site.seq_id_begin\n_ihm_sphere_obj_site.seq_id_end\n_ihm_sphere_obj_site.Cartn_x\n_ihm_sphere_obj_site.Cartn_y\n_ihm_sphere_obj_site.Cartn_z\n_ihm_sphere_obj_site.object_radius\n1 1 B 1 10 4.0 0.0 0.0 1.0\n#\nloop_\n_ihm_gaussian_obj_site.id\n_ihm_gaussian_obj_site.entity_id\n_ihm_gaussian_obj_site.asym_id\n_ihm_gaussian_obj_site.seq_id_begin\n_ihm_gaussian_obj_site.seq_id_end\n_ihm_gaussian_obj_site.mean_Cartn_x\n_ihm_gaussian_obj_site.mean_Cartn_y\n_ihm_gaussian_obj_site.mean_Cartn_z\n_ihm_gaussian_obj_site.weight\n_ihm_gaussian_obj_site.covariance_matrix[1][1]\n_ihm_gaussian_obj_site.covariance_matrix[1][2]\n_ihm_gaussian_obj_site.covariance_matrix[1][3]\n_ihm_gaussian_obj_site.covariance_matrix[2][1]\n_ihm_gaussian_obj_site.covariance_matrix[2][2]\n_ihm_gaussian_obj_site.covariance_matrix[2][3]\n_ihm_gaussian_obj_site.covariance_matrix[3][1]\n_ihm_gaussian_obj_site.covariance_matrix[3][2]\n_ihm_gaussian_obj_site.covariance_matrix[3][3]\n1 1 C 11 20 8.0 0.0 0.0 1.0 4.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = mol.atomic_structure();

    assert_eq!(
        structure
            .units
            .iter()
            .map(|unit| unit.kind)
            .collect::<Vec<_>>(),
        vec![UnitKind::Atomic, UnitKind::Spheres, UnitKind::Gaussians]
    );
    assert_eq!(structure.element_count, 3);
    assert_eq!(structure.position(0, 0).unwrap(), vec3(-2.0, 0.0, 0.0));
    assert_eq!(structure.position(1, 0).unwrap(), vec3(4.0, 0.0, 0.0));
    assert_eq!(structure.position(2, 0).unwrap(), vec3(8.0, 0.0, 0.0));
    assert_eq!(structure.coarse.hierarchy.spheres.elements[0].asym_id, "B");
    assert_eq!(
        structure.coarse.hierarchy.gaussians.elements[0].asym_id,
        "C"
    );
    assert!(structure.boundary.box_min.x <= -2.0);
    assert!(structure.boundary.box_max.x >= 8.0);
    assert!(structure.boundary.box_max.x < 10.0);
    assert_eq!(structure.lookup3d.find(vec3(4.0, 0.0, 0.0), 10.0).len(), 3);
    assert_eq!(structure.model.sequence.sequences.len(), 2);
    assert_eq!(
        structure.model.sequence.sequences[0].ranges,
        vec![SequenceRange {
            seq_id_begin: 1,
            seq_id_end: 10
        }]
    );
    assert_eq!(
        structure.model.sequence.sequences[1].ranges,
        vec![SequenceRange {
            seq_id_begin: 11,
            seq_id_end: 20
        }]
    );
}

#[test]
fn asymmetric_unit_order_and_traits_match_molstar_of_model() {
    let mut atoms = vec![
        test_atom(1, "N", "A", 1, vec3(0.0, 0.0, 0.0)),
        test_atom(2, "CA", "A", 1, vec3(1.0, 0.0, 0.0)),
        test_atom(3, "C", "A", 1, vec3(2.0, 0.0, 0.0)),
        test_atom(4, "O", "A", 1, vec3(2.2, 0.0, 0.0)),
        test_atom(5, "O", "W1", 1, vec3(4.0, 0.0, 0.0)),
        test_atom(6, "O", "W2", 1, vec3(5.0, 0.0, 0.0)),
        test_atom(7, "ZN", "I1", 1, vec3(7.0, 0.0, 0.0)),
        test_atom(8, "CL", "I2", 1, vec3(8.0, 0.0, 0.0)),
    ];
    for atom in &mut atoms[0..4] {
        atom.entity_id = "P".to_string();
    }
    for atom in &mut atoms[4..6] {
        atom.entity_id = "W".to_string();
        atom.residue = "HOH".to_string();
    }
    for atom in &mut atoms[6..8] {
        atom.entity_id = "I".to_string();
        atom.residue = "LIG".to_string();
        atom.auth_chain = "I".to_string();
    }
    let entities = vec![
        Entity {
            id: "P".to_string(),
            type_name: "polymer".to_string(),
            description: String::new(),
        },
        Entity {
            id: "W".to_string(),
            type_name: "water".to_string(),
            description: String::new(),
        },
        Entity {
            id: "I".to_string(),
            type_name: "non-polymer".to_string(),
            description: String::new(),
        },
    ];
    let structure = Molecule {
        atoms,
        entity_index: EntityIndexMap::from_entities(&entities, &[], &[]),
        entities,
        coarse_spheres: vec![CoarseSphere {
            id: 1,
            model_num: 1,
            entity_id: "S".to_string(),
            asym_id: "S".to_string(),
            seq_id_begin: 1,
            seq_id_end: 5,
            position: vec3(10.0, 0.0, 0.0),
            radius: 1.0,
            rmsf: 0.0,
        }],
        coarse_gaussians: vec![CoarseGaussian {
            id: 1,
            model_num: 1,
            entity_id: "G".to_string(),
            asym_id: "G".to_string(),
            seq_id_begin: 6,
            seq_id_end: 10,
            position: vec3(12.0, 0.0, 0.0),
            weight: 1.0,
            covariance: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }],
        ..Molecule::default()
    }
    .atomic_structure();

    assert_eq!(
        structure
            .units
            .iter()
            .map(|unit| unit.kind)
            .collect::<Vec<_>>(),
        vec![
            UnitKind::Atomic,
            UnitKind::Atomic,
            UnitKind::Atomic,
            UnitKind::Spheres,
            UnitKind::Gaussians,
        ]
    );
    assert_eq!(
        structure
            .units
            .iter()
            .map(|unit| unit.chain_indices.clone())
            .collect::<Vec<_>>(),
        vec![vec![0], vec![1, 2], vec![3, 4], vec![0], vec![0]]
    );
    assert_eq!(structure.units[0].traits, UnitTraits::NONE);
    assert_eq!(
        structure.units[1].traits,
        UnitTraits::WATER.union(UnitTraits::MULTI_CHAIN)
    );
    assert_eq!(structure.units[2].traits, UnitTraits::MULTI_CHAIN);
    assert_eq!(structure.units[3].traits, UnitTraits::NONE);
    assert_eq!(structure.units[4].traits, UnitTraits::NONE);
}

#[test]
fn cif_asymmetric_unit_order_follows_molstar_chain_segment_order() {
    let cif = b"data_demo\nloop_\n_entity.id\n_entity.type\n1 polymer\n2 non-polymer\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.auth_asym_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 1 A 0.0 0.0 0.0\nHETATM 2 C C1 LIG B 2 1 B 5.0 0.0 0.0\nATOM 3 C CA GLY C 1 1 C 10.0 0.0 0.0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let structure = mol.atomic_structure();
    let units = structure
        .units
        .iter()
        .map(|unit| {
            (
                structure.model.hierarchy.chains[unit.chain_index]
                    .id
                    .as_str(),
                unit.atom_indices.clone(),
            )
        })
        .collect::<Vec<_>>();

    assert_eq!(units, vec![("A", vec![0]), ("B", vec![1]), ("C", vec![2])]);
}

#[test]
fn structure_unit_array_ordering_and_unit_maps_match_molstar_create() {
    let cif = b"data_demo\nloop_\n_entity.id\n_entity.type\n1 polymer\n2 non-polymer\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.auth_asym_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 1 A 0.0 0.0 0.0\nHETATM 2 C C1 LIG B 2 1 B 5.0 0.0 0.0\nATOM 3 C CA GLY C 1 1 C 10.0 0.0 0.0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let mut structure = mol.atomic_structure();
    assert!(structure.units_are_molstar_sorted());
    structure.units[0].id = 20;
    structure.units[1].id = 10;
    structure.units[2].id = 30;
    assert!(!structure.units_are_molstar_sorted());

    crate::model::sort_structure_units_molstar(&mut structure.units);

    assert!(structure.units_are_molstar_sorted());
    assert_eq!(
        structure
            .units
            .iter()
            .map(|unit| unit.id)
            .collect::<Vec<_>>(),
        vec![10, 20, 30]
    );
    assert_eq!(structure.unit_index_by_id(20), Some(1));
    assert_eq!(structure.unit_index_map().get(&10), Some(&0));
    assert_eq!(structure.position(20, 0).unwrap(), vec3(0.0, 0.0, 0.0));
    assert_eq!(structure.position(10, 0).unwrap(), vec3(5.0, 0.0, 0.0));
    assert_eq!(structure.serial_index(20, 0), Some(1));
    assert_eq!(
        structure.serial_mapping(),
        (vec![0, 1, 2], vec![0, 1, 2], vec![1, 0, 2])
    );
}

#[test]
fn structure_element_and_bond_loci_remap_by_unit_id_and_element_identity() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 0.0 0.0 0.0\nATOM 2 C CB GLY A 1 1.0 0.0 0.0\nATOM 3 C C GLY A 1 2.0 0.0 0.0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let source = mol.atomic_structure();
    let mut target = source.clone();
    target.units[0].elements = vec![1, 2];
    target.units[0].atom_indices = vec![1, 2];
    target.element_count = 2;

    assert_eq!(
        source.remap_element_loci_to(&[(0, vec![0, 2])], &target),
        vec![(0, vec![1])]
    );
    assert_eq!(
        source.remap_bond_loci_to(&[(0, 1, 0, 2)], &target),
        vec![(0, 0, 0, 1)]
    );
    assert!(source
        .remap_bond_loci_to(&[(0, 0, 0, 2)], &target)
        .is_empty());
}

#[test]
fn partitioned_loci_extend_to_whole_chains_across_same_chain_operator_group() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 0.0 0.0 0.0\nATOM 2 C CB GLY A 1 1.0 0.0 0.0\nATOM 3 C C GLY A 1 2.0 0.0 0.0\nATOM 4 O O GLY A 1 3.0 0.0 0.0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let mut structure = mol.atomic_structure();
    let mut left = structure.units[0].clone();
    left.id = 0;
    left.chain_group_id = 7;
    left.traits = left.traits.union(UnitTraits::PARTITIONED);
    left.elements = vec![0, 1];
    left.atom_indices = vec![0, 1];
    let mut right = structure.units[0].clone();
    right.id = 1;
    right.chain_group_id = 7;
    right.traits = right.traits.union(UnitTraits::PARTITIONED);
    right.elements = vec![2, 3];
    right.atom_indices = vec![2, 3];
    assert!(left.are_same_chain_operator_group(&right));

    structure.units = vec![left, right];
    structure.element_count = 4;

    assert_eq!(
        structure.extend_element_loci_to_whole_chains(&[(0, vec![0])]),
        vec![(0, vec![0, 1]), (1, vec![0, 1])]
    );
}

#[test]
fn chem_comp_bond_dictionary_creates_bonds_without_distance_inference() {
    let cif = b"data_demo\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n_chem_comp_bond.pdbx_stereo_config\n_chem_comp_bond.pdbx_ordinal\nLIG C1 C2 DOUB N N 7\nLIG C2 C3 delo N E 8\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 LIG A 1 0.0 0.0 0.0\nHETATM 2 C C2 LIG A 1 20.0 0.0 0.0\nHETATM 3 C C3 LIG A 1 40.0 0.0 0.0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            infer_bonds: false,
            assembly: None,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bonds.len(), 2);
    assert_eq!(mol.chemical_component_bonds[1].stereo_config, "E");
    assert_eq!(mol.chemical_component_bonds[1].ordinal, Some(8));
    assert_eq!(mol.bond_metadata[0].source, BondSource::ChemComp);
    assert_eq!(mol.bond_metadata[0].order, 2);
    assert_eq!(mol.bond_metadata[0].key, 7);
    assert_eq!(mol.bond_metadata[1].key, 8);
    assert!(mol.bond_metadata[1].flags.contains(BondFlags::AROMATIC));
    assert!(mol.bond_metadata[1].flags.contains(BondFlags::RESONANCE));
}

#[test]
fn chem_comp_bond_missing_or_not_present_ordinal_uses_molstar_zero_key() {
    let missing_ordinal = b"data_demo\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\nLIG C1 C2 sing\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 LIG A 1 0.0 0.0 0.0\nHETATM 2 C C2 LIG A 1 20.0 0.0 0.0\n#\n";
    let not_present_ordinal = b"data_demo\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_ordinal\nLIG C1 C2 sing .\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 LIG A 1 0.0 0.0 0.0\nHETATM 2 C C2 LIG A 1 20.0 0.0 0.0\n#\n";

    for cif in [missing_ordinal.as_slice(), not_present_ordinal.as_slice()] {
        let mol = parse_molecule_with_options(
            cif,
            &MeshOptions {
                format: InputFormat::Cif,
                infer_bonds: false,
                assembly: None,
                ..MeshOptions::default()
            },
        )
        .unwrap();

        assert_eq!(mol.bonds.len(), 1);
        assert_eq!(mol.bond_metadata[0].source, BondSource::ChemComp);
        assert_eq!(mol.bond_metadata[0].key, 0);
    }
}

#[test]
fn chem_comp_bond_alt_loc_compatibility_matches_molstar() {
    let cif = b"data_demo\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_ordinal\nLIG C1 C2 sing 11\nLIG C0 C2 sing 12\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C0 . LIG A 1 -20.0 0.0 0.0\nHETATM 2 C C1 A LIG A 1 0.0 0.0 0.0\nHETATM 3 C C1 B LIG A 1 0.0 8.0 0.0\nHETATM 4 C C2 B LIG A 1 20.0 8.0 0.0\nHETATM 5 C C2 A LIG A 1 20.0 0.0 0.0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            alt_loc: "all".to_string(),
            infer_bonds: false,
            assembly: None,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let mut actual = mol
        .bonds
        .iter()
        .map(|bond| {
            let a = &mol.atoms[bond.a];
            let b = &mol.atoms[bond.b];
            (
                a.name.clone(),
                a.alt_id.clone(),
                b.name.clone(),
                b.alt_id.clone(),
            )
        })
        .collect::<Vec<_>>();
    actual.sort();

    assert_eq!(
        actual,
        vec![
            ("C0".into(), "".into(), "C2".into(), "A".into()),
            ("C0".into(), "".into(), "C2".into(), "B".into()),
            ("C1".into(), "A".into(), "C2".into(), "A".into()),
            ("C1".into(), "B".into(), "C2".into(), "B".into()),
        ]
    );
    assert_eq!(
        mol.bond_metadata
            .iter()
            .map(|metadata| (metadata.source.clone(), metadata.key))
            .collect::<Vec<_>>(),
        vec![
            (BondSource::ChemComp, 11),
            (BondSource::ChemComp, 11),
            (BondSource::ChemComp, 12),
            (BondSource::ChemComp, 12),
        ]
    );
}

#[test]
fn chem_comp_bond_normalizes_deuterium_atom_names_like_molstar() {
    let cif = b"data_demo\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\nHOH O H1 sing\nDOD O H1 sing\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 O O HOH A 1 0.0 0.0 0.0\nHETATM 2 D D1 HOH A 1 1.0 0.0 0.0\nHETATM 3 O O DOD B 1 0.0 4.0 0.0\nHETATM 4 D D1 DOD B 1 1.0 4.0 0.0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            infer_bonds: false,
            assembly: None,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bonds.len(), 1);
    assert_eq!((mol.bonds[0].a, mol.bonds[0].b), (0, 1));
    assert_eq!(mol.bond_metadata[0].source, BondSource::ChemComp);
}

#[test]
fn molstar_intra_component_bond_order_table_matches_reference_entries() {
    let entries = [
        ("HIS", "CD2", "CG"),
        ("HIS", "CE1", "ND1"),
        ("ARG", "CZ", "NH2"),
        ("PHE", "CE1", "CZ"),
        ("PHE", "CD2", "CE2"),
        ("PHE", "CD1", "CG"),
        ("TRP", "CD1", "CG"),
        ("TRP", "CD2", "CE2"),
        ("TRP", "CE3", "CZ3"),
        ("TRP", "CH2", "CZ2"),
        ("ASN", "CG", "OD1"),
        ("GLN", "CD", "OE1"),
        ("TYR", "CD1", "CG"),
        ("TYR", "CD2", "CE2"),
        ("TYR", "CE1", "CZ"),
        ("ASP", "CG", "OD1"),
        ("GLU", "CD", "OE1"),
        ("G", "C8", "N7"),
        ("G", "C4", "C5"),
        ("G", "C2", "N3"),
        ("G", "C6", "O6"),
        ("C", "C4", "N3"),
        ("C", "C5", "C6"),
        ("C", "C2", "O2"),
        ("A", "C2", "N3"),
        ("A", "C6", "N1"),
        ("A", "C4", "C5"),
        ("A", "C8", "N7"),
        ("U", "C5", "C6"),
        ("U", "C2", "O2"),
        ("U", "C4", "O4"),
        ("DG", "C8", "N7"),
        ("DG", "C4", "C5"),
        ("DG", "C2", "N3"),
        ("DG", "C6", "O6"),
        ("DC", "C4", "N3"),
        ("DC", "C5", "C6"),
        ("DC", "C2", "O2"),
        ("DA", "C2", "N3"),
        ("DA", "C6", "N1"),
        ("DA", "C4", "C5"),
        ("DA", "C8", "N7"),
        ("DT", "C5", "C6"),
        ("DT", "C2", "O2"),
        ("DT", "C4", "O4"),
    ];

    for (comp_id, atom_id_1, atom_id_2) in entries {
        assert_eq!(
            crate::model::intra_bond_order_from_table(comp_id, atom_id_1, atom_id_2),
            2,
            "{comp_id}|{atom_id_1}|{atom_id_2}"
        );
        assert_eq!(
            crate::model::intra_bond_order_from_table(comp_id, atom_id_2, atom_id_1),
            2,
            "{comp_id}|{atom_id_2}|{atom_id_1}"
        );
    }
    assert_eq!(
        crate::model::intra_bond_order_from_table("ALA", "C", "O"),
        2
    );
    assert_eq!(
        crate::model::intra_bond_order_from_table("A", "P", "OP1"),
        2
    );
    assert_eq!(
        crate::model::intra_bond_order_from_table("DG", "P", "OP1"),
        2
    );
    assert_eq!(
        crate::model::intra_bond_order_from_table("ALA", "N", "CA"),
        1
    );
}

#[test]
fn molstar_inter_component_bond_order_table_matches_reference_entries() {
    assert_eq!(
        crate::parser::inter_bond_order_from_table("LYS", "NZ", "RET", "C15"),
        2
    );
    assert_eq!(
        crate::parser::inter_bond_order_from_table("RET", "C15", "LYS", "NZ"),
        2
    );
    assert_eq!(
        crate::parser::inter_bond_order_from_table("LYS", "CA", "RET", "C15"),
        1
    );
}

#[test]
fn parses_chemical_component_dictionary_atoms_and_angles() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.name\n_chem_comp.type\n_chem_comp.formula\n_chem_comp.formula_weight\n_chem_comp.one_letter_code\n_chem_comp.three_letter_code\n_chem_comp.mon_nstd_flag\n_chem_comp.pdbx_synonyms\n_chem_comp.pdbx_formal_charge\n_chem_comp.pdbx_initial_date\n_chem_comp.pdbx_modified_date\n_chem_comp.pdbx_ambiguous_flag\n_chem_comp.pdbx_release_status\nLIG 'TEST LIGAND' non-polymer 'C2 O1' 42.0 ? LIG y 'TEST;LIG' -1 2020-01-01 2024-05-14 N REL\n#\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.alt_atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n_chem_comp_atom.pdbx_leaving_atom_flag\n_chem_comp_atom.pdbx_stereo_config\n_chem_comp_atom.model_Cartn_x\n_chem_comp_atom.model_Cartn_y\n_chem_comp_atom.model_Cartn_z\n_chem_comp_atom.pdbx_model_Cartn_x_ideal\n_chem_comp_atom.pdbx_model_Cartn_y_ideal\n_chem_comp_atom.pdbx_model_Cartn_z_ideal\nLIG C1 C1 C 1 Y N N 0.0 1.0 2.0 0.1 1.1 2.1\nLIG O1 O1 O -1 N Y R 3.0 4.0 5.0 3.1 4.1 5.1\n#\nloop_\n_chem_comp_angle.comp_id\n_chem_comp_angle.atom_id_1\n_chem_comp_angle.atom_id_2\n_chem_comp_angle.atom_id_3\n_chem_comp_angle.value_angle\n_chem_comp_angle.value_angle_esd\nLIG C1 O1 C1 109.5 0.8\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 LIG A 1 0.0 0.0 0.0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            infer_bonds: false,
            assembly: None,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.chemical_components.len(), 1);
    assert_eq!(mol.chemical_components[0].name, "TEST LIGAND");
    assert_eq!(mol.chemical_components[0].formula, "C2 O1");
    assert_eq!(mol.chemical_components[0].formula_weight, Some(42.0));
    assert_eq!(mol.chemical_components[0].one_letter_code, "");
    assert_eq!(mol.chemical_components[0].three_letter_code, "LIG");
    assert_eq!(mol.chemical_components[0].mon_nstd_flag, "y");
    assert_eq!(mol.chemical_components[0].pdbx_formal_charge, Some(-1));
    assert_eq!(mol.chemical_components[0].pdbx_release_status, "REL");
    assert_eq!(mol.chemical_component_atoms.len(), 2);
    assert_eq!(mol.chemical_component_atoms[0].type_symbol, "C");
    assert_eq!(mol.chemical_component_atoms[0].charge, 1);
    assert!(mol.chemical_component_atoms[0].aromatic);
    assert!(!mol.chemical_component_atoms[0].leaving_atom);
    assert_eq!(
        mol.chemical_component_atoms[0].model_cartn,
        Some(vec3(0.0, 1.0, 2.0))
    );
    assert_eq!(
        mol.chemical_component_atoms[1].ideal_cartn,
        Some(vec3(3.1, 4.1, 5.1))
    );
    assert_eq!(mol.chemical_component_angles.len(), 1);
    assert_eq!(mol.chemical_component_angles[0].atom_id_2, "O1");
    assert_eq!(mol.chemical_component_angles[0].value_angle, Some(109.5));

    let info =
        String::from_utf8(molecule_info(cif, br#"{"format":"cif","infer-bonds":false}"#).unwrap())
            .unwrap();
    assert!(info.contains("\"chem_comp_count\":1"));
    assert!(info.contains("\"chem_comp_atom_count\":2"));
    assert!(info.contains("\"chem_comp_angle_count\":1"));
}

#[test]
fn info_exposes_bond_flag_and_resonance_counts() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CG PHE A 1 1.0 0.0 0.0\nATOM 2 C CD1 PHE A 1 0.5 0.8 0.0\nATOM 3 C CE1 PHE A 1 -0.5 0.8 0.0\nATOM 4 C CZ PHE A 1 -1.0 0.0 0.0\nATOM 5 C CE2 PHE A 1 -0.5 -0.8 0.0\nATOM 6 C CD2 PHE A 1 0.5 -0.8 0.0\n#\nloop_\n_molstar_bond_site.atom_id_1\n_molstar_bond_site.atom_id_2\n_molstar_bond_site.value_order\n_molstar_bond_site.type_id\n1 2 sing covale\n2 3 sing covale\n3 4 sing covale\n4 5 sing covale\n5 6 sing covale\n1 6 sing covale\n#\n";
    let info =
        String::from_utf8(molecule_info(cif, br#"{"format":"cif","infer-bonds":false}"#).unwrap())
            .unwrap();

    assert!(info.contains("\"covalent\":6"));
    assert!(info.contains("\"aromatic\":6"));
    assert!(info.contains("\"resonance\":6"));
    assert!(info.contains("\"rings\":1"));
    assert!(info.contains("\"aromatic_rings\":1"));
    assert!(info.contains("\"delocalized_bonds\":6"));
}

#[test]
fn detects_aromatic_ring_resonance_from_index_pair_bonds() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CG PHE A 1 1.0 0.0 0.0\nATOM 2 C CD1 PHE A 1 0.5 0.8 0.0\nATOM 3 C CE1 PHE A 1 -0.5 0.8 0.0\nATOM 4 C CZ PHE A 1 -1.0 0.0 0.0\nATOM 5 C CE2 PHE A 1 -0.5 -0.8 0.0\nATOM 6 C CD2 PHE A 1 0.5 -0.8 0.0\n#\nloop_\n_molstar_bond_site.atom_id_1\n_molstar_bond_site.atom_id_2\n_molstar_bond_site.value_order\n_molstar_bond_site.type_id\n1 2 sing covale\n2 3 sing covale\n3 4 sing covale\n4 5 sing covale\n5 6 sing covale\n1 6 sing covale\n#\n";
    let mol = parse_molecule(cif, InputFormat::Cif).unwrap();
    assert_eq!(mol.resonance.ring_count, 1);
    assert_eq!(mol.resonance.aromatic_ring_count, 1);
    assert_eq!(mol.resonance.delocalized_bond_count, 6);
    assert!(mol
        .bond_metadata
        .iter()
        .all(|metadata| metadata.flags.contains(BondFlags::RESONANCE)));
    assert!(mol
        .bond_metadata
        .iter()
        .all(|metadata| metadata.flags.contains(BondFlags::AROMATIC)));
}

#[test]
fn amino_acid_aromatic_ring_reference_tests_cover_his_phe_trp_tyr() {
    type AromaticRingCase<'a> = (&'a str, &'a [&'a str], &'a [(usize, usize)], usize, usize);
    let cases: &[AromaticRingCase<'_>] = &[
        (
            "HIS",
            &["CG", "ND1", "CE1", "NE2", "CD2"],
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
            1,
            5,
        ),
        (
            "PHE",
            &["CG", "CD1", "CE1", "CZ", "CE2", "CD2"],
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
            1,
            6,
        ),
        (
            "TRP",
            &["CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ3", "CH2", "CZ2"],
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 0),
                (4, 5),
                (5, 6),
                (6, 7),
                (7, 8),
                (8, 3),
            ],
            2,
            10,
        ),
        (
            "TYR",
            &["CG", "CD1", "CE1", "CZ", "CE2", "CD2"],
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
            1,
            6,
        ),
    ];

    for (residue, atom_names, bond_pairs, expected_rings, expected_delocalized) in cases {
        let molecule = topology_test_molecule(residue, atom_names, bond_pairs);

        assert_eq!(molecule.resonance.ring_count, *expected_rings, "{residue}");
        assert_eq!(
            molecule.resonance.aromatic_ring_count, *expected_rings,
            "{residue}"
        );
        assert_eq!(
            molecule.resonance.delocalized_bond_count, *expected_delocalized,
            "{residue}"
        );
        assert!(molecule.rings.iter().all(|ring| ring.aromatic), "{residue}");
        assert!(molecule.bond_metadata.iter().all(|metadata| {
            metadata.flags.contains(BondFlags::AROMATIC)
                && metadata.flags.contains(BondFlags::RESONANCE)
        }));
    }
}

#[test]
fn unflagged_non_aromatic_five_and_six_member_rings_stay_non_aromatic() {
    type NonAromaticRingCase<'a> = (&'a str, &'a [&'a str], &'a [(usize, usize)]);
    let cases: &[NonAromaticRingCase<'_>] = &[
        (
            "L5",
            &["C1", "C2", "N1", "O1", "S1"],
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
        ),
        (
            "L6",
            &["C1", "C2", "C3", "N1", "O1", "S1"],
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
        ),
    ];

    for (residue, atom_names, bond_pairs) in cases {
        let mut molecule = topology_test_molecule(residue, atom_names, bond_pairs);
        if let Some(atom) = molecule.atoms.get_mut(1) {
            atom.position.y = 10.0;
        }
        if let Some(atom) = molecule.atoms.last_mut() {
            atom.position.z = 10.0;
        }
        molecule.refresh_topology_metadata();

        assert_eq!(molecule.resonance.ring_count, 1, "{residue}");
        assert_eq!(molecule.resonance.aromatic_ring_count, 0, "{residue}");
        assert!(!molecule.rings[0].aromatic, "{residue}");
        assert_eq!(molecule.resonance.delocalized_bond_count, 0, "{residue}");
        assert!(molecule.bond_metadata.iter().all(|metadata| {
            !metadata.flags.contains(BondFlags::AROMATIC)
                && !metadata.flags.contains(BondFlags::RESONANCE)
        }));
    }
}

#[test]
fn chem_comp_bond_aromatic_component_graph_marks_ring_aromatic() {
    let cif = b"data_demo\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\nLIG C1 C2 sing Y\nLIG C2 C3 sing Y\nLIG C3 C4 sing Y\nLIG C4 C5 sing Y\nLIG C5 C6 sing Y\nLIG C6 C1 sing Y\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 LIG A 1 0.0 0.0 0.0\nHETATM 2 C C2 LIG A 1 20.0 0.0 0.0\nHETATM 3 C C3 LIG A 1 40.0 0.0 0.0\nHETATM 4 C C4 LIG A 1 60.0 0.0 0.0\nHETATM 5 C C5 LIG A 1 80.0 0.0 0.0\nHETATM 6 C C6 LIG A 1 100.0 0.0 0.0\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(molecule.bonds.len(), 6);
    assert!(molecule
        .bond_metadata
        .iter()
        .all(|metadata| metadata.source == BondSource::ChemComp));
    assert_eq!(molecule.rings.len(), 1);
    assert!(molecule.rings[0].aromatic);
    assert_eq!(molecule.resonance.aromatic_ring_count, 1);
    assert_eq!(molecule.resonance.delocalized_bond_count, 6);
    assert!(molecule.bond_metadata.iter().all(|metadata| {
        metadata.flags.contains(BondFlags::AROMATIC)
            && metadata.flags.contains(BondFlags::RESONANCE)
    }));
}

#[test]
fn ligand_fused_ring_reference_tracks_membership_fingerprint_and_aromaticity() {
    let cif = b"data_demo\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\nLIG C1 C2 arom Y\nLIG C2 C3 arom Y\nLIG C3 C4 arom Y\nLIG C4 C5 arom Y\nLIG C5 C6 arom Y\nLIG C6 C1 arom Y\nLIG C4 C7 arom Y\nLIG C7 C8 arom Y\nLIG C8 C9 arom Y\nLIG C9 C10 arom Y\nLIG C10 C5 arom Y\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 LIG A 1 0.000 0.000 0.000\nHETATM 2 C C2 LIG A 1 1.000 0.000 0.000\nHETATM 3 C C3 LIG A 1 1.500 0.866 0.000\nHETATM 4 C C4 LIG A 1 1.000 1.732 0.000\nHETATM 5 C C5 LIG A 1 0.000 1.732 0.000\nHETATM 6 C C6 LIG A 1 -0.500 0.866 0.000\nHETATM 7 C C7 LIG A 1 1.500 2.598 0.000\nHETATM 8 C C8 LIG A 1 1.000 3.464 0.000\nHETATM 9 C C9 LIG A 1 0.000 3.464 0.000\nHETATM 10 C C10 LIG A 1 -0.500 2.598 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(molecule.bonds.len(), 11);
    assert!(molecule
        .bond_metadata
        .iter()
        .all(|metadata| metadata.source == BondSource::ChemComp));
    assert_eq!(molecule.resonance.ring_count, 2);
    assert_eq!(molecule.resonance.aromatic_ring_count, 2);
    assert_eq!(molecule.resonance.delocalized_bond_count, 11);
    assert_eq!(molecule.rings.len(), 2);
    assert_eq!(molecule.rings[0].atom_indices, vec![0, 1, 2, 3, 4, 5]);
    assert_eq!(molecule.rings[1].atom_indices, vec![3, 4, 6, 7, 8, 9]);
    assert!(molecule.rings.iter().all(|ring| ring.aromatic));
    assert!(molecule
        .rings
        .iter()
        .all(|ring| ring.fingerprint == "C-C-C-C-C-C"));
    assert_eq!(
        molecule.resonance.element_ring_indices,
        vec![
            vec![0],
            vec![0],
            vec![0],
            vec![0, 1],
            vec![0, 1],
            vec![0],
            vec![1],
            vec![1],
            vec![1],
            vec![1],
        ]
    );
    assert_eq!(
        molecule.resonance.element_aromatic_ring_indices,
        molecule.resonance.element_ring_indices
    );
    assert_eq!(molecule.resonance.ring_component_index, vec![0, 0]);
    assert_eq!(molecule.resonance.ring_components, vec![vec![0, 1]]);
    assert!(molecule.bond_metadata.iter().all(|metadata| {
        metadata.flags.contains(BondFlags::AROMATIC)
            && metadata.flags.contains(BondFlags::RESONANCE)
    }));
}

#[test]
fn aromatic_ring_bond_flags_take_precedence_over_size_and_fallback_rules() {
    let all_aromatic_four_member = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 LIG A 1 0.000 0.000 0.000\nHETATM 2 C C2 LIG A 1 1.000 0.000 0.000\nHETATM 3 C C3 LIG A 1 1.000 1.000 0.000\nHETATM 4 C C4 LIG A 1 0.000 1.000 0.000\n#\nloop_\n_molstar_bond_site.atom_id_1\n_molstar_bond_site.atom_id_2\n_molstar_bond_site.value_order\n_molstar_bond_site.type_id\n1 2 arom covale\n2 3 arom covale\n3 4 arom covale\n4 1 arom covale\n#\n";
    let molecule = parse_molecule_with_options(
        all_aromatic_four_member,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(molecule.rings.len(), 1);
    assert!(molecule.rings[0].aromatic);
    assert_eq!(molecule.resonance.aromatic_ring_count, 1);
    assert_eq!(molecule.resonance.delocalized_bond_count, 4);
    assert!(molecule.resonance.delocalized_triplets.is_empty());
    assert!(molecule
        .bond_metadata
        .iter()
        .all(|metadata| metadata.flags.contains(BondFlags::RESONANCE)));

    let partial_aromatic_phe = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CG PHE A 1 1.0 0.0 0.0\nATOM 2 C CD1 PHE A 1 0.5 0.8 0.0\nATOM 3 C CE1 PHE A 1 -0.5 0.8 0.0\nATOM 4 C CZ PHE A 1 -1.0 0.0 0.0\nATOM 5 C CE2 PHE A 1 -0.5 -0.8 0.0\nATOM 6 C CD2 PHE A 1 0.5 -0.8 0.0\n#\nloop_\n_molstar_bond_site.atom_id_1\n_molstar_bond_site.atom_id_2\n_molstar_bond_site.value_order\n_molstar_bond_site.type_id\n1 2 arom covale\n2 3 sing covale\n3 4 sing covale\n4 5 sing covale\n5 6 sing covale\n1 6 sing covale\n#\n";
    let molecule = parse_molecule_with_options(
        partial_aromatic_phe,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(molecule.rings.len(), 1);
    assert!(!molecule.rings[0].aromatic);
    assert_eq!(molecule.resonance.aromatic_ring_count, 0);
    assert_eq!(molecule.resonance.delocalized_bond_count, 1);
    assert!(molecule.resonance.delocalized_triplets.is_empty());
    assert!(molecule.bond_metadata[0]
        .flags
        .contains(BondFlags::AROMATIC));
    assert!(!molecule.bond_metadata[0]
        .flags
        .contains(BondFlags::RESONANCE));
    assert!(molecule
        .bond_metadata
        .iter()
        .skip(1)
        .all(|metadata| !metadata.flags.contains(BondFlags::AROMATIC)));
}

#[test]
fn non_ring_aromatic_centers_build_molstar_delocalized_triplets() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 LIG A 1 -1.000 0.000 0.000\nHETATM 2 C C2 LIG A 1 0.000 0.000 0.000\nHETATM 3 C C3 LIG A 1 1.000 0.000 0.000\nHETATM 4 N N1 LIG A 1 0.000 1.000 0.000\n#\nloop_\n_molstar_bond_site.atom_id_1\n_molstar_bond_site.atom_id_2\n_molstar_bond_site.value_order\n_molstar_bond_site.type_id\n1 2 arom covale\n2 3 arom covale\n2 4 arom covale\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert!(molecule.rings.is_empty());
    assert_eq!(
        molecule.resonance.delocalized_triplets,
        vec![[0, 1, 2], [0, 1, 2], [0, 1, 3]]
    );
    assert_eq!(
        molecule.resonance.delocalized_triplet_lookup.triplets,
        molecule.resonance.delocalized_triplets
    );
    assert_eq!(
        molecule
            .resonance
            .delocalized_triplet_lookup
            .get_third_element(0, 1),
        Some(2)
    );
    assert_eq!(
        molecule
            .resonance
            .delocalized_triplet_lookup
            .get_third_element(1, 0),
        Some(2)
    );
    assert_eq!(
        molecule
            .resonance
            .delocalized_triplet_lookup
            .get_third_element(1, 3),
        Some(0)
    );
    assert_eq!(
        molecule
            .resonance
            .delocalized_triplet_lookup
            .get_triplet_indices(1),
        Some([0, 1, 2].as_slice())
    );
    assert_eq!(
        molecule
            .resonance
            .delocalized_triplet_lookup
            .get_triplet_indices(0),
        None
    );
    assert_eq!(molecule.resonance.delocalized_bond_count, 3);
}

#[test]
fn molstar_ring_detection_keeps_small_and_proline_rings_non_aromatic() {
    let triangle = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 LIG A 1 0.000 0.000 0.000\nHETATM 2 C C2 LIG A 1 1.000 0.000 0.000\nHETATM 3 C C3 LIG A 1 0.500 0.800 0.000\n#\nloop_\n_molstar_bond_site.atom_id_1\n_molstar_bond_site.atom_id_2\n_molstar_bond_site.value_order\n_molstar_bond_site.type_id\n1 2 sing covale\n2 3 sing covale\n3 1 sing covale\n#\n";
    let triangle = parse_molecule_with_options(
        triangle,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(triangle.resonance.ring_count, 1);
    assert_eq!(triangle.rings[0].atom_indices, vec![0, 1, 2]);
    assert!(!triangle.rings[0].aromatic);
    assert_eq!(triangle.resonance.aromatic_ring_count, 0);
    assert_eq!(triangle.resonance.delocalized_bond_count, 0);
    assert!(triangle
        .bond_metadata
        .iter()
        .all(|metadata| !metadata.flags.contains(BondFlags::RESONANCE)));

    let proline = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 N N PRO A 1 0.000 0.000 0.000\nATOM 2 C CA PRO A 1 1.000 0.000 0.000\nATOM 3 C CB PRO A 1 1.500 0.800 0.000\nATOM 4 C CG PRO A 1 0.700 1.300 0.000\nATOM 5 C CD PRO A 1 -0.200 0.800 0.000\n#\nloop_\n_molstar_bond_site.atom_id_1\n_molstar_bond_site.atom_id_2\n_molstar_bond_site.value_order\n_molstar_bond_site.type_id\n1 2 sing covale\n2 3 sing covale\n3 4 sing covale\n4 5 sing covale\n5 1 sing covale\n#\n";
    let proline = parse_molecule_with_options(
        proline,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(proline.resonance.ring_count, 1);
    assert!(!proline.rings[0].aromatic);
    assert_eq!(proline.resonance.aromatic_ring_count, 0);
    assert_eq!(proline.resonance.delocalized_bond_count, 0);
}

#[test]
fn molstar_ring_fingerprint_and_membership_lookup_match_unit_rings_contract() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 N N1 LIG A 1 0.000 0.000 0.000\nHETATM 2 O O1 LIG A 1 1.000 0.000 0.000\nHETATM 3 C C1 LIG A 1 0.500 0.800 0.000\n#\nloop_\n_molstar_bond_site.atom_id_1\n_molstar_bond_site.atom_id_2\n_molstar_bond_site.value_order\n_molstar_bond_site.type_id\n1 2 sing covale\n2 3 sing covale\n3 1 sing covale\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(molecule.rings.len(), 1);
    assert_eq!(molecule.rings[0].fingerprint, "C-N-O");
    assert_eq!(
        molecule.resonance.element_ring_indices,
        vec![vec![0], vec![0], vec![0]]
    );
    assert_eq!(
        molecule.resonance.element_aromatic_ring_indices,
        vec![Vec::<usize>::new(), Vec::new(), Vec::new()]
    );
    assert_eq!(molecule.resonance.ring_component_index, vec![0]);
    assert_eq!(molecule.resonance.ring_components, vec![vec![0]]);
}

#[test]
fn molstar_ring_fingerprint_uses_sorted_ring_indices_not_discovery_path() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 LIG A 1 0.000 0.000 0.000\nHETATM 2 N N1 LIG A 1 1.000 0.000 0.000\nHETATM 3 C C2 LIG A 1 1.000 1.000 0.000\nHETATM 4 O O1 LIG A 1 0.000 1.000 0.000\n#\nloop_\n_molstar_bond_site.atom_id_1\n_molstar_bond_site.atom_id_2\n_molstar_bond_site.value_order\n_molstar_bond_site.type_id\n1 3 sing covale\n3 2 sing covale\n2 4 sing covale\n4 1 sing covale\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(molecule.rings.len(), 1);
    assert_eq!(molecule.rings[0].atom_indices, vec![0, 1, 2, 3]);
    assert_eq!(molecule.rings[0].fingerprint, "C-N-C-O");
}

#[test]
fn ring_order_is_stable_for_different_bond_input_orders() {
    let atoms = ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]
        .iter()
        .enumerate()
        .map(|(index, name)| {
            let mut atom = test_atom(index + 1, name, "A", 1, vec3(index as f32, 0.0, 0.0));
            atom.residue = "PHE".to_string();
            atom.auth_residue = "PHE".to_string();
            atom
        })
        .collect::<Vec<_>>();
    let metadata = vec![
        BondMetadata {
            source: BondSource::IndexPair,
            order: 1,
            flags: BondFlags::COVALENT,
            key: -1,
            distance: None,
            operator_a: -1,
            operator_b: -1,
            struct_conn: None,
        };
        6
    ];
    let bonds_a = vec![
        Bond { a: 0, b: 1 },
        Bond { a: 1, b: 2 },
        Bond { a: 2, b: 3 },
        Bond { a: 3, b: 4 },
        Bond { a: 4, b: 5 },
        Bond { a: 0, b: 5 },
    ];
    let bonds_b = vec![
        Bond { a: 3, b: 4 },
        Bond { a: 0, b: 5 },
        Bond { a: 1, b: 2 },
        Bond { a: 4, b: 5 },
        Bond { a: 0, b: 1 },
        Bond { a: 2, b: 3 },
    ];
    let mut molecule_a = Molecule {
        atoms: atoms.clone(),
        bonds: bonds_a,
        bond_metadata: metadata.clone(),
        ..Molecule::default()
    };
    let mut molecule_b = Molecule {
        atoms,
        bonds: bonds_b,
        bond_metadata: metadata,
        ..Molecule::default()
    };
    molecule_a.refresh_topology_metadata();
    molecule_b.refresh_topology_metadata();

    assert_eq!(molecule_a.rings, molecule_b.rings);
    assert_eq!(molecule_a.rings[0].atom_indices, vec![0, 1, 2, 3, 4, 5]);
    assert!(molecule_a.rings[0].aromatic);
}

#[test]
fn molstar_ring_detection_suppresses_superset_cycles_with_chords() {
    let molecule = topology_test_molecule(
        "LIG",
        &["C1", "C2", "C3", "C4"],
        &[(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)],
    );

    assert_eq!(
        molecule
            .rings
            .iter()
            .map(|ring| ring.atom_indices.clone())
            .collect::<Vec<_>>(),
        vec![vec![0, 1, 2], vec![0, 2, 3]]
    );
    assert!(molecule.rings.iter().all(|ring| !ring.aromatic));
}

#[test]
fn molstar_ring_detection_uses_only_covalent_bonds() {
    let atoms = ["C1", "C2", "C3"]
        .iter()
        .enumerate()
        .map(|(index, name)| test_atom(index + 1, name, "A", 1, vec3(index as f32, 0.0, 0.0)))
        .collect::<Vec<_>>();
    let mut molecule = Molecule {
        atoms,
        bonds: vec![
            Bond { a: 0, b: 1 },
            Bond { a: 1, b: 2 },
            Bond { a: 2, b: 0 },
        ],
        bond_metadata: vec![
            BondMetadata {
                source: BondSource::IndexPair,
                flags: BondFlags::COVALENT,
                ..BondMetadata::computed()
            },
            BondMetadata {
                source: BondSource::IndexPair,
                flags: BondFlags::HYDROGEN_BOND,
                ..BondMetadata::computed()
            },
            BondMetadata {
                source: BondSource::IndexPair,
                flags: BondFlags::COVALENT,
                ..BondMetadata::computed()
            },
        ],
        ..Molecule::default()
    };
    molecule.refresh_topology_metadata();

    assert!(molecule.rings.is_empty());
    assert_eq!(molecule.resonance.ring_count, 0);
}

#[test]
fn molstar_ring_detection_splits_alt_loc_ring_searches() {
    let mut atoms = ["C1", "C2", "C3", "C2"]
        .iter()
        .enumerate()
        .map(|(index, name)| test_atom(index + 1, name, "A", 1, vec3(index as f32, 0.0, 0.0)))
        .collect::<Vec<_>>();
    atoms[1].alt_id = "A".to_string();
    atoms[3].alt_id = "B".to_string();
    let mut molecule = Molecule {
        atoms,
        bonds: vec![
            Bond { a: 0, b: 1 },
            Bond { a: 1, b: 2 },
            Bond { a: 2, b: 0 },
            Bond { a: 0, b: 3 },
            Bond { a: 3, b: 2 },
        ],
        bond_metadata: vec![
            BondMetadata {
                source: BondSource::IndexPair,
                flags: BondFlags::COVALENT,
                ..BondMetadata::computed()
            };
            5
        ],
        ..Molecule::default()
    };
    molecule.refresh_topology_metadata();

    assert_eq!(
        molecule
            .rings
            .iter()
            .map(|ring| ring.atom_indices.clone())
            .collect::<Vec<_>>(),
        vec![vec![0, 1, 2], vec![0, 2, 3]]
    );
}

#[test]
fn bond_ring_fixture_graphs_match_molstar_reference_summary() {
    let expected =
        include_str!("../../tests/expected/bond-ring-graph-reference-summary.json").trim();
    let actual = bond_ring_reference_summary_json();
    assert_eq!(actual, expected);
}

#[test]
fn molstar_helix_uses_oriented_elliptical_tube() {
    let pdb = include_bytes!("../../tests/fixtures/pdb/assembly-altloc-helix.pdb");
    let molecule = parse_molecule_with_options(
        pdb,
        &MeshOptions {
            format: InputFormat::Pdb,
            assembly: None,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };
    let objects = build_render_objects(&molecule, &options);

    let helix_segments = objects
        .iter()
        .filter_map(|object| {
            if let RenderObject::PolymerTraceSegment {
                widths,
                heights,
                kind: PolymerTraceSegmentKind::Tube { profile, .. },
                ..
            } = object
            {
                Some((*profile, widths[1], heights[1]))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    assert_eq!(helix_segments.len(), 3);
    assert!(helix_segments.iter().all(|(profile, width, height)| {
        *profile == PolymerProfile::Elliptical
            && (*width - 0.20).abs() <= 0.000_001
            && (*height - 1.00).abs() <= 0.000_001
    }));
}

#[test]
fn tubular_helices_use_molstar_helix_orientation_centers() {
    let mut atoms = Vec::new();
    for (seq, ca) in [
        (1, vec3(1.0, 0.0, 0.0)),
        (2, vec3(0.0, 1.0, 1.5)),
        (3, vec3(-1.0, 0.0, 3.0)),
        (4, vec3(0.0, -1.0, 4.5)),
        (5, vec3(1.0, 0.0, 6.0)),
    ] {
        atoms.push(test_atom(
            atoms.len() + 1,
            "N",
            "H",
            seq,
            ca - vec3(0.25, 0.0, 0.0),
        ));
        atoms.push(test_atom(atoms.len() + 1, "CA", "H", seq, ca));
        atoms.push(test_atom(
            atoms.len() + 1,
            "C",
            "H",
            seq,
            ca + vec3(0.25, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "O",
            "H",
            seq,
            ca + vec3(0.25, 0.6, 0.0),
        ));
    }
    let molecule = Molecule {
        atoms,
        helices: vec![SecondaryRange {
            chain: "H".to_string(),
            start: 1,
            start_insertion_code: String::new(),
            end: 5,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        tubular_helices: true,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let objects = build_render_objects(&molecule, &options);
    let helix_segments = objects
        .iter()
        .filter_map(|object| {
            if let RenderObject::PolymerTraceSegment {
                controls,
                widths,
                heights,
                ..
            } = object
            {
                Some((
                    controls.p2.to_vec3(),
                    controls.d12.to_vec3(),
                    widths[1],
                    heights[1],
                ))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    assert_eq!(helix_segments.len(), 5);
    assert!(helix_segments
        .iter()
        .all(|(_, normal, _, _)| *normal == vec3(1.0, 0.0, 0.0)));
    assert!(helix_segments
        .iter()
        .all(|(_, _, width, height)| (*width - 1.50).abs() <= 0.000_001
            && (*height - 1.50).abs() <= 0.000_001));
    assert_vec3_close(helix_segments[2].0, vec3(0.0, 0.0, 3.0), 0.000_1);
    assert!(helix_segments
        .iter()
        .all(|(center, _, _, _)| vec3_is_finite(*center)));
}

#[test]
fn molstar_trace_quality_options_are_parsed() {
    let Some(molstar_polymer_trace) =
        read_molstar_source("mol-repr/structure/visual/polymer-trace-mesh.ts")
    else {
        eprintln!("skipping pinned Mol* polymer trace source audit; artifacts is absent");
        return;
    };
    assert!(molstar_polymer_trace
        .contains("linearSegments: PD.Numeric(8, { min: 1, max: 48, step: 1 }"));
    assert!(molstar_polymer_trace
        .contains("radialSegments: PD.Numeric(16, { min: 2, max: 56, step: 2 }"));
    assert!(molstar_polymer_trace.contains("sizeFactor: PD.Numeric(0.2"));

    let defaults = MeshOptions::default();
    assert_eq!(defaults.linear_segments, 8);
    assert_eq!(defaults.radial_segments, 16);
    assert_eq!(defaults.ribbon_radius, 0.2);
    assert_eq!(defaults.sheet_arrow_factor, 1.5);
    assert!(!defaults.tubular_helices);
    assert!(!defaults.round_cap);
    assert_eq!(defaults.block_index, None);
    assert_eq!(defaults.block_header, None);
    assert_eq!(defaults.color_theme, ColorTheme::ChainId);

    let options = MeshOptions::from_json(
            br#"{"color-theme":"chain-id","tubular-helices":true,"linear-segments":6,"radial-segments":12,"sheet-arrow-factor":0.75,"block-index":1,"block-header":"second"}"#,
        )
        .unwrap();

    assert_eq!(options.color_theme, ColorTheme::ChainId);
    assert!(options.tubular_helices);
    assert_eq!(options.quality, Some(VisualQuality::Auto));
    assert_eq!(options.linear_segments, 6);
    assert_eq!(options.radial_segments, 12);
    assert_eq!(options.sheet_arrow_factor, 0.75);
    assert_eq!(options.block_index, Some(1));
    assert_eq!(options.block_header.as_deref(), Some("second"));
    let resolved = options.resolved_for_quality(2_870, false);
    assert_eq!(resolved.sphere_detail, 2);
    assert_eq!(resolved.linear_segments, 10);
    assert_eq!(resolved.radial_segments, 20);

    let clamped = MeshOptions::from_json(
        br#"{"linear-segments":0,"radial-segments":1,"sheet-arrow-factor":9}"#,
    )
    .unwrap();
    assert_eq!(clamped.linear_segments, 1);
    assert_eq!(clamped.radial_segments, 2);
    assert_eq!(clamped.ribbon_radius, 0.2);
    assert_eq!(clamped.sheet_arrow_factor, 3.0);

    let clamped =
        MeshOptions::from_json(br#"{"linear-segments":99,"radial-segments":99,"ribbon-radius":9}"#)
            .unwrap();
    assert_eq!(clamped.linear_segments, 48);
    assert_eq!(clamped.radial_segments, 56);
    assert_eq!(clamped.ribbon_radius, 2.0);

    assert_eq!(
        MeshOptions::from_json(br#"{"color-theme":"residue-name"}"#).unwrap_err(),
        "unsupported color-theme: residue-name; expected \"chain-id\""
    );
}

#[test]
fn molstar_helix_normals_keep_carbonyl_orientation_continuous() {
    let mut atoms = Vec::new();
    for (i, carbonyl_y) in [1.0, -1.0, 1.0].into_iter().enumerate() {
        let seq = i as i32 + 1;
        let x = i as f32 * 1.5;
        atoms.push(test_atom(
            atoms.len() + 1,
            "CA",
            "A",
            seq,
            vec3(x, 0.0, 0.0),
        ));
        atoms.push(test_atom(atoms.len() + 1, "C", "A", seq, vec3(x, 0.1, 0.0)));
        atoms.push(test_atom(
            atoms.len() + 1,
            "O",
            "A",
            seq,
            vec3(x, 0.1 + carbonyl_y, 0.0),
        ));
    }
    let molecule = Molecule {
        atoms,
        helices: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 1,
            start_insertion_code: String::new(),
            end: 3,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let normals = build_render_objects(&molecule, &options)
        .into_iter()
        .filter_map(|object| {
            if let RenderObject::PolymerTraceSegment { controls, .. } = object {
                Some(controls.d12.to_vec3())
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    assert_eq!(normals.len(), 3);
    assert!(normals.iter().all(|normal| vec3_is_finite(*normal)));
    assert!(normals
        .iter()
        .all(|normal| (normal.length() - 1.0).abs() < 0.000_1));
    assert!(normals.iter().all(|normal| normal.y.abs() > 0.9));
}

#[test]
fn oriented_ribbon_mesh_has_finite_molstar_tube_counts() {
    let mut mesh = Mesh::default();
    add_oriented_ribbon(
        &mut mesh,
        &[
            Vec3 {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            Vec3 {
                x: 1.0,
                y: 0.2,
                z: 0.0,
            },
            Vec3 {
                x: 2.0,
                y: 0.0,
                z: 0.1,
            },
        ],
        &[
            Vec3 {
                x: 0.0,
                y: 1.0,
                z: 0.0,
            },
            Vec3 {
                x: 0.0,
                y: 1.0,
                z: 0.1,
            },
            Vec3 {
                x: 0.0,
                y: 1.0,
                z: 0.0,
            },
        ],
        0.20,
        1.00,
    );

    assert_eq!(mesh.vertices.len(), 866);
    assert_eq!(mesh.normals.len(), 866);
    assert_eq!(mesh.faces.len(), 1600);
    assert!(mesh.vertices.iter().all(|v| vec3_is_finite(*v)));
    assert!(mesh.normals.iter().all(|v| vec3_is_finite(*v)));
    assert!(faces_have_valid_indices(&mesh));

    let mut rounded = Mesh::default();
    add_profile_tube_for_test(
        &mut rounded,
        vec![vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 1.0)],
        vec![vec3(0.0, 1.0, 0.0), vec3(0.0, 1.0, 0.0)],
        vec![vec3(-1.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0)],
        vec![0.2, 0.2],
        vec![1.0, 1.0],
        8,
        TestTubeProfile::Rounded,
        false,
        false,
        false,
    );
    assert_eq!(rounded.vertices.len(), 16);
    assert_eq!(rounded.faces.len(), 16);
    assert_vec3_close(
        rounded.vertices[0],
        vec3(-0.076_536_69, 0.984_775_9, 0.0),
        0.000_01,
    );
    assert_vec3_close(
        rounded.normals[0],
        vec3(-0.382_683_43, 0.923_879_5, 0.0),
        0.000_01,
    );
}

#[test]
fn oriented_ribbon_mesh_preserves_elliptical_profile_axes() {
    let mut mesh = Mesh::default();
    add_oriented_ribbon(
        &mut mesh,
        &[vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0)],
        &[vec3(0.0, 1.0, 0.0), vec3(0.0, 1.0, 0.0)],
        0.20,
        1.00,
    );

    let (min, max) = mesh_bounds(&mesh);
    let y_span = max.y - min.y;
    let z_span = max.z - min.z;
    assert!(y_span > 1.8, "y_span={y_span}");
    assert!(z_span < 0.5, "z_span={z_span}");
    assert!(y_span > z_span * 4.0, "y_span={y_span} z_span={z_span}");
    assert!(mesh.vertices.iter().all(|v| vec3_is_finite(*v)));
    assert!(mesh.normals.iter().all(|v| vec3_is_finite(*v)));
    assert!(faces_have_valid_indices(&mesh));
}

#[test]
fn curve_segment_interpolation_matches_molstar_reference_vectors() {
    let controls = CurveSegmentControls {
        sec_struc_first: false,
        sec_struc_last: true,
        p0: DVec3::from_vec3(vec3(-1.2, 0.1, 0.3)),
        p1: DVec3::from_vec3(vec3(0.0, 0.0, 0.0)),
        p2: DVec3::from_vec3(vec3(1.0, 0.25, -0.1)),
        p3: DVec3::from_vec3(vec3(2.1, -0.15, 0.4)),
        p4: DVec3::from_vec3(vec3(3.0, 0.2, 0.0)),
        d12: DVec3::from_vec3(vec3(0.0, 1.0, 0.2)),
        d23: DVec3::from_vec3(vec3(0.3, 0.8, -0.1)),
    };
    let mut state = CurveSegmentState::new(4);

    interpolate_curve_segment(&mut state, &controls, 0.9, 0.5);
    interpolate_sizes(&mut state, 0.3, 0.5, 0.9, 0.7, 1.1, 1.4, 0.5);

    assert_vec3_close(
        state.curve_points[0],
        vec3(0.511_25, 0.158_75, -0.14),
        0.000_01,
    );
    assert_vec3_close(state.curve_points[2], vec3(1.0, 0.25, -0.1), 0.000_01);
    assert_vec3_close(
        state.curve_points[4],
        vec3(1.556_25, 0.043_75, 0.168_75),
        0.000_01,
    );
    assert_vec3_close(
        state.tangent_vectors[2],
        vec3(0.979_941_96, -0.069_948_67, 0.186_603_77),
        0.000_01,
    );
    assert_vec3_close(
        state.normal_vectors[1],
        vec3(-0.228_962_1, 0.903_766_04, 0.281_519_56),
        0.000_01,
    );
    assert_vec3_close(
        state.normal_vectors[3],
        vec3(0.291_595_52, 0.917_723_54, 0.069_312_7),
        0.000_01,
    );
    assert_vec3_close(
        state.binormal_vectors[2],
        vec3(-0.196_282_95, -0.176_913_1, 0.964_455_66),
        0.000_01,
    );
    assert!((state.width_values[0] - 0.4).abs() <= 0.000_001);
    assert!((state.width_values[4] - 0.7).abs() <= 0.000_001);
    assert!((state.height_values[0] - 0.9).abs() <= 0.000_001);
    assert!((state.height_values[4] - 1.25).abs() <= 0.000_001);
}

#[test]
fn profile_tube_extrusion_preserves_molstar_frame_scale() {
    let mut mesh = Mesh::default();
    add_profile_tube_for_test(
        &mut mesh,
        vec![vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 1.0)],
        vec![vec3(0.0, 2.0, 0.0), vec3(0.0, 2.0, 0.0)],
        vec![vec3(0.0, 0.0, 0.5), vec3(0.0, 0.0, 0.5)],
        vec![0.3, 0.3],
        vec![0.7, 0.7],
        4,
        TestTubeProfile::Elliptical,
        false,
        false,
        false,
    );

    assert_eq!(mesh.vertices.len(), 8);
    assert_eq!(mesh.faces.len(), 8);
    assert_vec3_close(mesh.vertices[0], vec3(0.0, 1.4, 0.0), 0.000_01);
    assert_vec3_close(mesh.vertices[1], vec3(0.0, 0.0, 0.15), 0.000_01);
    assert_vec3_close(mesh.normals[0], vec3(0.0, 1.0, 0.0), 0.000_01);
    assert_vec3_close(mesh.normals[1], vec3(0.0, 0.0, 1.0), 0.000_01);
    assert_eq!(
        (mesh.faces[0].a, mesh.faces[0].b, mesh.faces[0].c),
        (1, 5, 0)
    );
    assert_eq!(
        (mesh.faces[1].a, mesh.faces[1].b, mesh.faces[1].c),
        (5, 4, 0)
    );
    assert!(mesh.vertices.iter().all(|v| vec3_is_finite(*v)));
    assert!(mesh.normals.iter().all(|v| vec3_is_finite(*v)));
    assert!(faces_have_valid_indices(&mesh));
}

#[test]
fn profile_tube_round_caps_match_molstar_single_segment_double_cap() {
    let mut mesh = Mesh::default();
    add_profile_tube_for_test(
        &mut mesh,
        vec![vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 1.0)],
        vec![vec3(0.0, 1.0, 0.0), vec3(0.0, 1.0, 0.0)],
        vec![vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0)],
        vec![0.2, 0.2],
        vec![0.2, 0.2],
        4,
        TestTubeProfile::Elliptical,
        true,
        true,
        true,
    );

    let molstar_epsilon_radius = 0.2 * f64::EPSILON as f32;
    assert_eq!(mesh.vertices.len(), 18);
    assert_eq!(mesh.normals.len(), 18);
    assert_eq!(mesh.faces.len(), 16);
    assert_vec3_close(
        mesh.vertices[0],
        vec3(0.0, molstar_epsilon_radius, 0.0),
        0.000_000_000_000_001,
    );
    assert_vec3_close(mesh.normals[0], vec3(0.0, 0.0, 1.0), 0.000_001);
    assert_vec3_close(
        mesh.vertices[4],
        vec3(0.0, molstar_epsilon_radius, 1.0),
        0.000_000_000_000_001,
    );
    assert_vec3_close(mesh.normals[4], vec3(0.0, 0.0, -1.0), 0.000_001);
    assert_vec3_close(mesh.vertices[8], vec3(0.0, 0.0, 0.0), 0.000_001);
    assert_vec3_close(mesh.vertices[9], vec3(0.0, 0.0, 0.0), 0.000_001);
    assert_vec3_close(mesh.vertices[13], vec3(0.0, 0.0, 1.0), 0.000_001);
    assert_eq!(
        (mesh.faces[8].a, mesh.faces[8].b, mesh.faces[8].c),
        (10, 9, 8)
    );
    assert_eq!(
        (mesh.faces[12].a, mesh.faces[12].b, mesh.faces[12].c),
        (14, 15, 13)
    );
    assert!(faces_have_valid_indices(&mesh));
}

#[test]
fn tube_path_uses_molstar_tube_builder_vertex_and_cap_order() {
    let mut mesh = Mesh::default();
    add_tube_path_for_test(
        &mut mesh,
        &[vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 1.0)],
        0.2,
        4,
    );

    assert_eq!(mesh.vertices.len(), 30);
    assert_eq!(mesh.normals.len(), 30);
    assert_eq!(mesh.faces.len(), 40);
    assert_vec3_close(mesh.vertices[0], vec3(0.0, 0.2, 0.0), 0.000_01);
    assert_vec3_close(mesh.vertices[1], vec3(-0.2, 0.0, 0.0), 0.000_01);
    assert_vec3_close(mesh.normals[0], vec3(0.0, 1.0, 0.0), 0.000_01);
    assert_vec3_close(mesh.normals[1], vec3(-1.0, 0.0, 0.0), 0.000_01);
    assert_eq!(
        (mesh.faces[0].a, mesh.faces[0].b, mesh.faces[0].c),
        (1, 5, 0)
    );
    assert_eq!(
        (mesh.faces[1].a, mesh.faces[1].b, mesh.faces[1].c),
        (5, 4, 0)
    );
    assert_vec3_close(mesh.vertices[20], vec3(0.0, 0.0, 0.0), 0.000_01);
    assert_eq!(
        (mesh.faces[32].a, mesh.faces[32].b, mesh.faces[32].c),
        (22, 21, 20)
    );
    assert_vec3_close(mesh.vertices[25], vec3(0.0, 0.0, 1.0), 0.000_01);
    assert_eq!(
        (mesh.faces[36].a, mesh.faces[36].b, mesh.faces[36].c),
        (26, 27, 25)
    );
    assert!(faces_have_valid_indices(&mesh));
}

#[test]
fn ribbon_builder_matches_molstar_add_ribbon_vertices() {
    let mut mesh = Mesh::default();
    add_ribbon_for_test(
        &mut mesh,
        &[vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 1.0)],
        0.2,
        0.7,
        1,
    );

    assert_eq!(mesh.vertices.len(), 8);
    assert_eq!(mesh.normals.len(), 8);
    assert_eq!(mesh.faces.len(), 4);
    assert_vec3_close(mesh.vertices[0], vec3(0.0, 0.7, 0.0), 0.000_01);
    assert_vec3_close(mesh.vertices[1], vec3(0.0, -0.7, 0.0), 0.000_01);
    assert_vec3_close(mesh.vertices[4], vec3(0.0, 0.7, 1.0), 0.000_01);
    assert_vec3_close(mesh.normals[0], vec3(1.0, 0.0, 0.0), 0.000_01);
    assert_vec3_close(mesh.normals[2], vec3(-1.0, 0.0, 0.0), 0.000_01);
    assert_eq!(
        (mesh.faces[0].a, mesh.faces[0].b, mesh.faces[0].c),
        (0, 5, 1)
    );
    assert_eq!(
        (mesh.faces[1].a, mesh.faces[1].b, mesh.faces[1].c),
        (0, 4, 5)
    );
    assert_eq!(
        (mesh.faces[2].a, mesh.faces[2].b, mesh.faces[2].c),
        (3, 7, 2)
    );
    assert_eq!(
        (mesh.faces[3].a, mesh.faces[3].b, mesh.faces[3].c),
        (2, 7, 6)
    );
}

#[test]
fn sheet_builder_matches_molstar_add_sheet_vertices_and_caps() {
    let mut mesh = Mesh::default();
    add_sheet_for_test(
        &mut mesh,
        &[vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 1.0)],
        0.2,
        0.7,
        0.0,
        true,
        true,
        1,
    );

    assert_eq!(mesh.vertices.len(), 24);
    assert_eq!(mesh.normals.len(), 24);
    assert_eq!(mesh.faces.len(), 12);
    assert_vec3_close(mesh.vertices[0], vec3(-0.2, 0.7, 0.0), 0.000_01);
    assert_vec3_close(mesh.vertices[1], vec3(0.2, 0.7, 0.0), 0.000_01);
    assert_vec3_close(mesh.vertices[4], vec3(0.2, -0.7, 0.0), 0.000_01);
    assert_vec3_close(mesh.normals[0], vec3(0.0, 1.0, 0.0), 0.000_01);
    assert_vec3_close(mesh.normals[2], vec3(1.0, 0.0, 0.0), 0.000_01);
    assert_vec3_close(mesh.normals[6], vec3(-1.0, 0.0, 0.0), 0.000_01);
    assert_eq!(
        (mesh.faces[0].a, mesh.faces[0].b, mesh.faces[0].c),
        (0, 9, 1)
    );
    assert_eq!(
        (mesh.faces[1].a, mesh.faces[1].b, mesh.faces[1].c),
        (0, 8, 9)
    );
    assert_eq!(
        (mesh.faces[8].a, mesh.faces[8].b, mesh.faces[8].c),
        (18, 17, 16)
    );
    assert_eq!(
        (mesh.faces[10].a, mesh.faces[10].b, mesh.faces[10].c),
        (20, 21, 22)
    );
    assert!(faces_have_valid_indices(&mesh));
}

#[test]
fn sheet_builder_keeps_molstar_arrow_normal_offset_unnormalized() {
    let mut mesh = Mesh::default();
    add_sheet_for_test(
        &mut mesh,
        &[vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 1.0)],
        0.2,
        0.7,
        0.5,
        false,
        false,
        1,
    );

    assert_eq!(mesh.vertices.len(), 24);
    assert_eq!(mesh.normals.len(), 24);
    assert_eq!(mesh.faces.len(), 12);
    assert_vec3_close(mesh.vertices[0], vec3(-0.2, 0.5, 0.0), 0.000_01);
    assert_vec3_close(mesh.vertices[4], vec3(0.2, -0.5, 0.0), 0.000_01);
    assert_vec3_close(mesh.normals[0], vec3(0.0, 1.0, 0.5), 0.000_01);
    assert_vec3_close(mesh.normals[4], vec3(0.0, -1.0, 0.5), 0.000_01);
    assert!(mesh.normals[0].length() > 1.0);
    assert_eq!(
        (mesh.faces[8].a, mesh.faces[8].b, mesh.faces[8].c),
        (18, 17, 16)
    );
    assert!(faces_have_valid_indices(&mesh));
}

#[test]
fn molstar_sheet_and_tube_defaults_use_size_factor_aspect_and_caps() {
    let molecule = sheet_test_molecule();
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };
    let objects = build_render_objects(&molecule, &options);
    let sheet_segments = objects
        .iter()
        .filter_map(|object| {
            if let RenderObject::PolymerTraceSegment {
                widths,
                heights,
                kind: PolymerTraceSegmentKind::Sheet { arrow_height },
                start_cap,
                end_cap,
                ..
            } = object
            {
                Some((widths[1], heights[1], *arrow_height, *start_cap, *end_cap))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    assert!(!sheet_segments.is_empty());
    assert!(sheet_segments
        .iter()
        .all(
            |(width, height, _, _, _)| (*width - 0.20).abs() <= 0.000_001
                && (*height - 1.00).abs() <= 0.000_001
        ));
    assert!(sheet_segments
        .iter()
        .all(|(_, _, arrow_height, _, _)| arrow_height.abs() <= 0.000_001));
    assert!(sheet_segments.first().unwrap().3);
    assert!(sheet_segments.last().unwrap().4);

    let pdb = b"ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C\nATOM      2  CA  GLY A   2       1.000   0.000   0.000  1.00 10.00           C\nEND\n";
    let molecule = parse_molecule(pdb, InputFormat::Pdb).unwrap();
    let objects = build_render_objects(
        &molecule,
        &MeshOptions {
            representation: Representation::Cartoon,
            center: false,
            assembly: None,
            ..MeshOptions::default()
        },
    );
    let RenderObject::PolymerTraceSegment {
        widths, heights, ..
    } = objects
        .iter()
        .find(|object| matches!(object, RenderObject::PolymerTraceSegment { .. }))
        .expect("cartoon trace segment")
    else {
        unreachable!()
    };
    assert!((widths[1] - 0.20).abs() <= 0.000_001);
    assert!((heights[1] - 0.20).abs() <= 0.000_001);

    let cartoon_objects = build_render_objects(
        &molecule,
        &MeshOptions {
            representation: Representation::Cartoon,
            ribbon_radius: 1.0,
            center: false,
            assembly: None,
            ..MeshOptions::default()
        },
    );
    let RenderObject::PolymerTraceSegment {
        widths, heights, ..
    } = cartoon_objects
        .iter()
        .find(|object| matches!(object, RenderObject::PolymerTraceSegment { .. }))
        .expect("cartoon trace segment")
    else {
        unreachable!()
    };
    assert!((widths[1] - 0.20).abs() <= 0.000_001);
    assert!((heights[1] - 0.20).abs() <= 0.000_001);

    let ribbon_objects = build_render_objects(
        &molecule,
        &MeshOptions {
            representation: Representation::Ribbon,
            ribbon_radius: 1.0,
            center: false,
            assembly: None,
            ..MeshOptions::default()
        },
    );
    let RenderObject::Tube { radius, .. } = ribbon_objects
        .iter()
        .find(|object| matches!(object, RenderObject::Tube { .. }))
        .expect("ribbon polymer-tube")
    else {
        unreachable!()
    };
    assert!((*radius - 1.0).abs() <= 0.000_001);
}

#[test]
fn cartoon_trace_widths_use_molstar_uniform_size_factor() {
    let pdb = b"ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N\nATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00 10.00           C\nATOM      3  C   ALA A   1       1.500   1.000   0.000  1.00 10.00           C\nATOM      4  N   GLY A   2       2.500   1.000   0.000  1.00 10.00           N\nATOM      5  CA  GLY A   2       3.000   2.000   0.000  1.00 10.00           C\nATOM      6  C   GLY A   2       4.000   2.000   0.000  1.00 10.00           C\nEND\n";
    let molecule = parse_molecule(pdb, InputFormat::Pdb).unwrap();
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let objects = build_render_objects(&molecule, &options);
    let RenderObject::PolymerTraceSegment {
        widths, heights, ..
    } = objects
        .iter()
        .find(|object| matches!(object, RenderObject::PolymerTraceSegment { .. }))
        .expect("expected cartoon coil trace segment")
    else {
        unreachable!()
    };

    assert!((widths[1] - 0.20).abs() <= 0.000_001);
    assert!((heights[1] - 0.20).abs() <= 0.000_001);
}

#[test]
fn polymer_trace_widths_sample_uniform_previous_current_next_sizes() {
    let mut atoms = [
        test_atom(1, "CA", "A", 1, vec3(0.0, 0.0, 0.0)),
        test_atom(2, "CA", "A", 2, vec3(1.5, 0.0, 0.0)),
        test_atom(3, "CA", "A", 3, vec3(3.0, 0.0, 0.0)),
    ];
    atoms[1].element = "N".to_string();
    atoms[1].type_symbol = "N".to_string();
    atoms[2].element = "O".to_string();
    atoms[2].type_symbol = "O".to_string();
    let molecule = Molecule {
        atoms: atoms.into_iter().collect(),
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let segments = build_render_objects(&molecule, &options)
        .into_iter()
        .filter_map(|object| {
            if let RenderObject::PolymerTraceSegment {
                controls,
                widths,
                heights,
                overhang_width,
                ..
            } = object
            {
                Some((controls.p2.to_vec3(), widths, heights, overhang_width))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    assert_eq!(segments.len(), 3);
    let (_, widths, heights, overhang_width) = segments
        .iter()
        .find(|(center, _, _, _)| center.distance(vec3(1.5, 0.0, 0.0)) < 0.000_001)
        .expect("middle trace segment");

    for (actual, expected) in widths.iter().zip([0.20, 0.20, 0.20]) {
        assert!((*actual - expected).abs() <= 0.000_001);
    }
    assert_eq!(heights, widths);
    assert!((*overhang_width - 0.20).abs() <= 0.000_001);
}

#[test]
fn polymer_trace_widths_ignore_type_symbol_for_uniform_size_theme() {
    let mut atoms = [
        test_atom(1, "CA", "A", 1, vec3(0.0, 0.0, 0.0)),
        test_atom(2, "CA", "A", 2, vec3(1.5, 0.0, 0.0)),
        test_atom(3, "CA", "A", 3, vec3(3.0, 0.0, 0.0)),
    ];
    atoms[1].element = "C".to_string();
    atoms[1].type_symbol = "H".to_string();
    let molecule = Molecule {
        atoms: atoms.into_iter().collect(),
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let RenderObject::PolymerTraceSegment {
        controls, widths, ..
    } = build_render_objects(&molecule, &options)
        .into_iter()
        .find(|object| {
            matches!(
                object,
                RenderObject::PolymerTraceSegment { controls, .. }
                    if controls.p2.to_vec3().distance(vec3(1.5, 0.0, 0.0)) < 0.000_001
            )
        })
        .expect("middle trace segment")
    else {
        unreachable!()
    };

    assert!(controls.p2.to_vec3().distance(vec3(1.5, 0.0, 0.0)) < 0.000_001);
    assert!((widths[1] - 0.20).abs() <= 0.000_001);
}

#[test]
fn spacefill_spheres_use_type_symbol_for_physical_size_theme() {
    let mut atom = test_atom(1, "CA", "A", 1, vec3(0.0, 0.0, 0.0));
    atom.element = "C".to_string();
    atom.type_symbol = "H".to_string();
    let molecule = Molecule {
        atoms: vec![atom],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Spacefill,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let RenderObject::Sphere { radius, .. } = build_render_objects(&molecule, &options)
        .into_iter()
        .find(|object| matches!(object, RenderObject::Sphere { .. }))
        .expect("spacefill sphere")
    else {
        unreachable!()
    };

    assert!((radius - 1.1).abs() <= 0.000_001);
}

#[test]
fn nine_r1o_molstar_ply_regression_is_finite() {
    let pdb = include_bytes!("../../../package/examples/data/9R1O.pdb");
    let molecule = parse_pdb(std::str::from_utf8(pdb).unwrap()).unwrap();
    assert_eq!(molecule.atoms.len(), 2870);
    assert_eq!(molecule.helices.len(), 12);
    assert_eq!(molecule.sheets.len(), 0);
    assert_eq!(molecule.assemblies.len(), 2);
    assert_eq!(
        molecule
            .atoms
            .iter()
            .filter(|atom| atom.het && atom.residue == "1PE")
            .count(),
        13
    );

    let ply = String::from_utf8(
            convert_to_ply(
                pdb,
                br#"{"format":"pdb","representation":"molstar","assembly":"asymmetric-unit","sphere-detail":1}"#,
            )
            .unwrap(),
        )
        .unwrap();
    let vertex_count = ply_header_count(&ply, "element vertex ");
    let face_count = ply_header_count(&ply, "element face ");
    if std::env::var("MOLFIG_PRINT_9R1O_PLY_CONTRACT").as_deref() == Ok("1") {
        let (bounds_min, bounds_max) = ply_vertex_bounds(&ply);
        eprintln!("vertex_count={vertex_count}");
        eprintln!("face_count={face_count}");
        eprintln!(
            "group_count={}",
            ply_header_comment_usize(&ply, "comment molfig_group_count ")
        );
        eprintln!("byte_len={}", ply.len());
        eprintln!("fnv1a64={:016x}", stable_test_hash64(ply.as_bytes()));
        eprintln!("bounds_min={}", contract_triplet(bounds_min));
        eprintln!("bounds_max={}", contract_triplet(bounds_max));
        eprintln!("vertex_samples={}", ply_vertex_samples_json(&ply));
        eprintln!("face_samples={}", ply_face_samples_json(&ply));
    }
    assert!(vertex_count > 10_000);
    assert!(face_count > 10_000);
    assert!(ply_vertices_are_finite(&ply));
    assert!(ply_faces_have_valid_indices(&ply));

    let contract =
        include_str!("../../tests/expected/ply/9r1o-molstar-asymmetric-unit.ply.contract");
    assert_eq!(
        contract_value(contract, "fixture"),
        "package/examples/data/9R1O.pdb"
    );
    assert_eq!(contract_value(contract, "format"), "ascii-ply");
    assert_eq!(vertex_count, contract_usize(contract, "vertex_count"));
    assert_eq!(face_count, contract_usize(contract, "face_count"));
    assert_eq!(
        ply_header_comment_usize(&ply, "comment molfig_group_count "),
        contract_usize(contract, "group_count")
    );
    assert_eq!(ply.len(), contract_usize(contract, "byte_len"));
    assert_eq!(
        format!("{:016x}", stable_test_hash64(ply.as_bytes())),
        contract_value(contract, "fnv1a64")
    );
    assert_eq!(
        ply_vertex_samples_json(&ply),
        contract_value(contract, "vertex_samples")
    );
    assert_eq!(
        ply_face_samples_json(&ply),
        contract_value(contract, "face_samples")
    );
}

#[test]
fn nine_r1o_default_assembly_molstar_mesh_stays_near_reference_density() {
    let pdb = include_bytes!("../../../package/examples/data/9R1O.pdb");
    let Some(reference_obj) = read_repo_file_if_present("package/examples/data/9R1O.obj") else {
        eprintln!(
            "skipping 9R1O reference density audit; package/examples/data/9R1O.obj is absent"
        );
        return;
    };
    let ply = String::from_utf8(
        convert_to_ply(
            pdb,
            br#"{"format":"pdb","representation":"molstar","sphere-detail":1}"#,
        )
        .unwrap(),
    )
    .unwrap();

    let vertex_count = ply_header_count(&ply, "element vertex ");
    let face_count = ply_header_count(&ply, "element face ");
    let reference_obj_stats = obj_stats(&reference_obj);
    let reference_obj_vertices = obj_vectors(&reference_obj, "v ");
    let reference_obj_normals = obj_vectors(&reference_obj, "vn ");
    let (bounds_min, bounds_max) = ply_vertex_bounds(&ply);
    let vertex_samples = ply_vertex_samples(&ply);
    let obj = String::from_utf8(
        convert_to_obj(
            pdb,
            br#"{"format":"pdb","representation":"molstar","sphere-detail":1,"obj-basename":"9R1O","operator-metadata":false,"obj-groups":false}"#,
        )
        .unwrap(),
    )
    .unwrap();
    let obj_normals = obj_vectors(&obj, "vn ");
    let sampled_normal_deltas = normal_sample_abs_deltas_json(&obj_normals, &reference_obj_normals);
    if std::env::var("MOLFIG_PRINT_9R1O_DEFAULT_ASSEMBLY_PLY").as_deref() == Ok("1") {
        eprintln!("vertex_count={vertex_count}");
        eprintln!("face_count={face_count}");
        eprintln!(
            "group_count={}",
            ply_header_comment_usize(&ply, "comment molfig_group_count ")
        );
        eprintln!("byte_len={}", ply.len());
        eprintln!("fnv1a64={:016x}", stable_test_hash64(ply.as_bytes()));
        eprintln!("bounds_min={}", contract_triplet(bounds_min));
        eprintln!("bounds_max={}", contract_triplet(bounds_max));
        eprintln!("vertex_samples={}", vertex_samples_json(&vertex_samples));
        eprintln!("face_samples={}", ply_face_samples_json(&ply));
        eprintln!(
            "obj_bounds_abs_delta_min={}",
            f32_triplet_json(abs_delta_array(bounds_min, reference_obj_stats.min))
        );
        eprintln!(
            "obj_bounds_abs_delta_max={}",
            f32_triplet_json(abs_delta_array(bounds_max, reference_obj_stats.max))
        );
        eprintln!(
            "obj_vertex_sample_abs_delta_vs_ply=[{}]",
            vertex_sample_abs_deltas_json(&vertex_samples, &reference_obj_vertices)
        );
        eprintln!("obj_normal_sample_abs_delta=[{sampled_normal_deltas}]");
    }
    assert!(
        (90_000..=115_000).contains(&vertex_count),
        "vertex_count={vertex_count}"
    );
    assert!(
        (160_000..=195_000).contains(&face_count),
        "face_count={face_count}"
    );
    assert_eq!(
        vertex_count, reference_obj_stats.vertex_count,
        "default 9R1O Molstar PLY vertex count should match the pinned Mol* OBJ export"
    );
    assert_eq!(
        face_count, reference_obj_stats.face_count,
        "default 9R1O Molstar PLY face count should match the pinned Mol* OBJ export"
    );
    assert_eq!(
        f32_triplet_json(abs_delta_array(bounds_min, reference_obj_stats.min)),
        "[0.0003,0.0001,0.0005]",
        "default 9R1O bounds min delta vs pinned Mol* OBJ changed"
    );
    assert_eq!(
        f32_triplet_json(abs_delta_array(bounds_max, reference_obj_stats.max)),
        "[0.0005,0,0.0002]",
        "default 9R1O bounds max delta vs pinned Mol* OBJ changed"
    );
    assert_eq!(
        vertex_sample_abs_deltas_json(&vertex_samples, &reference_obj_vertices),
        r#"{"label":"start","molfig_index":0,"reference_index":0,"abs_delta":[0.0002,0.0002,0.0004]},{"label":"quarter","molfig_index":25278,"reference_index":25278,"abs_delta":[0.0002,0.0005,0.0002]},{"label":"middle","molfig_index":50557,"reference_index":50557,"abs_delta":[0.0001,0.0001,0.0002]},{"label":"three_quarter","molfig_index":75835,"reference_index":75835,"abs_delta":[0.0004,0.0005,0.0004]},{"label":"end","molfig_index":101113,"reference_index":101113,"abs_delta":[0,0.0003,0.0005]}"#,
        "default 9R1O sampled vertex deltas vs pinned Mol* OBJ changed"
    );
    assert_eq!(
        obj_normals.len(),
        reference_obj_normals.len(),
        "default 9R1O Molstar OBJ normal count should match the pinned Mol* OBJ export"
    );
    assert_eq!(
        sampled_normal_deltas,
        r#"{"label":"start","molfig_index":0,"reference_index":0,"abs_delta":[0,0,0]},{"label":"quarter","molfig_index":25278,"reference_index":25278,"abs_delta":[0,0,0]},{"label":"middle","molfig_index":50557,"reference_index":50557,"abs_delta":[0,0,0]},{"label":"three_quarter","molfig_index":75835,"reference_index":75835,"abs_delta":[0,0,0]},{"label":"end","molfig_index":101113,"reference_index":101113,"abs_delta":[0,0,0]}"#,
        "default 9R1O sampled normal deltas vs pinned Mol* OBJ changed"
    );
    assert!(ply_vertices_are_finite(&ply));
    assert!(ply_faces_have_valid_indices(&ply));
}

#[test]
#[ignore = "generates full 9R1O OBJ/STL diagnostics on demand"]
fn nine_r1o_default_assembly_live_export_diff_report_is_actionable() {
    let Some(actual) = nine_r1o_default_assembly_live_export_diff_json() else {
        eprintln!("skipping 9R1O live export diff report; reference OBJ/STL artifacts are absent");
        return;
    };
    if std::env::var("MOLFIG_PRINT_9R1O_DEFAULT_ASSEMBLY_LIVE_DIFF").as_deref() == Ok("1") {
        eprintln!("{actual}");
    }

    assert!(actual.contains(r#""name": "9r1o-default-assembly-live-export-diff""#));
    assert!(actual.contains(r#""source": "live molfig export""#));
    assert!(actual.contains(r#""source": "pinned Mol* geo-export""#));
    assert!(actual.contains(r#""obj_counts_delta": {"vertex": 0, "normal": 0, "face": 0}"#));
    assert!(actual.contains(r#""stl_counts_delta": {"facet": 0}"#));
    assert!(actual.contains(r#""obj_sample_abs_delta": ["#));
    assert!(actual.contains(r#""normal_abs_delta":"#));
    assert!(
        actual.matches(r#""normal_abs_delta":[0,0,0]"#).count() >= 5,
        "all sampled OBJ normals should match the pinned Mol* OBJ reference exactly"
    );
    assert!(actual.contains(r#""stl_sample_abs_delta": ["#));
    assert!(actual.contains(r#""stl_sample_signed_delta": ["#));
    assert!(actual.contains(r#""stl_first_facet_signed_delta": {"#));
    assert!(actual.contains(r#""vertex_centroid_signed_delta":"#));
    assert!(actual.contains(r#""vertex_residual_signed_delta":"#));
    assert!(actual.contains(r#""vertices_abs_delta":"#));
    assert!(actual.contains(r#""vertices_signed_delta":"#));
}

#[test]
#[ignore = "generates full 9R1O render-object span diagnostics on demand"]
fn nine_r1o_default_assembly_render_object_spans_are_actionable() {
    let pdb = include_bytes!("../../../package/examples/data/9R1O.pdb");
    let options =
        MeshOptions::from_json(br#"{"format":"pdb","representation":"molstar","sphere-detail":1}"#)
            .unwrap();
    let molecule = parse_molecule_with_options(pdb, &options).unwrap();
    let actual = render_object_span_summary_json(&molecule, &options);
    if std::env::var("MOLFIG_PRINT_9R1O_DEFAULT_ASSEMBLY_SPANS").as_deref() == Ok("1") {
        eprintln!("{actual}");
    }

    assert!(actual.contains(r#""vertex_count":101114"#));
    assert!(actual.contains(r#""normal_count":101114"#));
    assert!(actual.contains(r#""face_count":178864"#));
    assert!(actual.contains(r#""vertex_end":101114"#));
    assert!(actual.contains(r#""stl_facet_start":0,"stl_facet_end":1260"#));
    assert!(actual.contains(r#""index":0,"geometry_type":"tube","visual":"polymer-trace","representation":"molstar","secondary_type":"helix","chain":"A","residue_start":3,"residue_end":3,"group_id":0,"polymer_trace":{"initial":true,"final":false,"sec_struc_first":false,"sec_struc_last":false},"vertex_start":0,"vertex_end":253,"face_start":0,"face_end":420,"stl_facet_start":0,"stl_facet_end":1260"#));
}

#[test]
fn molstar_bond_cylinder_export_uses_two_open_half_cylinders() {
    let molecule = Molecule {
        atoms: vec![
            test_atom(1, "C1", "A", 1, vec3(0.0, 0.0, 0.0)),
            test_atom(2, "C2", "A", 1, vec3(1.5, 0.0, 0.0)),
        ],
        bonds: vec![Bond { a: 0, b: 1 }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::BallAndStick,
        center: false,
        assembly: None,
        quality: Some(VisualQuality::Higher),
        ..MeshOptions::default()
    };

    let mesh = build_mesh(&molecule, &options);

    // Mol* geo-export converts the bond cylinder impostor into two uncapped
    // half-cylinders. With primitivesQuality:auto and fewer than 2000
    // cylinders, each half uses 36 radial segments: 74 vertices and 72 faces.
    let sphere_vertices = 2 * 642;
    let sphere_faces = 2 * 1280;
    let bond_vertices = 2 * 74;
    let bond_faces = 2 * 72;
    assert_eq!(mesh.vertices.len(), sphere_vertices + bond_vertices);
    assert_eq!(mesh.faces.len(), sphere_faces + bond_faces);

    let first_cylinder_normal = mesh.normals[sphere_vertices];
    assert!(
        (first_cylinder_normal.length() - 1.0).abs() < 0.000_01,
        "Mol* cylinder exporter keeps radius in primitive props, so side normals are not scaled by 1/radius"
    );
}

#[test]
fn molfig_obj_outputs_are_diffed_against_all_molstar_reference_fixtures() {
    let manifest = include_str!("../../tests/expected/molstar-reference/reference-fixtures.txt");
    let mut checked = 0;

    for contract_path in reference_manifest_entries_for_format(manifest, "obj") {
        let contract = read_manifest_file(contract_path);
        let reference = read_manifest_file(contract_value(&contract, "obj_reference"));
        let reference_stats = obj_stats(&reference);
        assert_eq!(
            reference_stats.vertex_count,
            contract_usize(&contract, "obj_vertex_count"),
            "{contract_path}: reference vertex_count"
        );
        assert_eq!(
            reference_stats.normal_count,
            contract_usize(&contract, "obj_normal_count"),
            "{contract_path}: reference normal_count"
        );
        assert_eq!(
            reference_stats.face_count,
            contract_usize(&contract, "obj_face_count"),
            "{contract_path}: reference face_count"
        );

        let fixture = read_manifest_file_bytes(contract_value(&contract, "fixture"));
        let options = contract_value(&contract, "options");
        let generated = String::from_utf8(convert_to_obj(&fixture, options.as_bytes()).unwrap())
            .unwrap_or_else(|_| panic!("{contract_path}: OBJ output is not valid UTF-8"));
        let actual_stats = obj_stats(&generated);

        assert!(
            actual_stats.vertex_count > 0,
            "{contract_path}: OBJ vertices"
        );
        assert_eq!(
            actual_stats.vertex_count, actual_stats.normal_count,
            "{contract_path}: OBJ normals"
        );
        assert!(actual_stats.face_count > 0, "{contract_path}: OBJ faces");
        assert_ratio_in_range(
            actual_stats.vertex_count,
            reference_stats.vertex_count,
            0.25,
            1.25,
            &format!("{contract_path}: OBJ vertex count"),
        );
        assert_ratio_in_range(
            actual_stats.face_count,
            reference_stats.face_count,
            0.25,
            1.25,
            &format!("{contract_path}: OBJ face count"),
        );
        let report = crate::diff_text(
            &reference,
            &generated,
            &format!("{contract_path}: OBJ exact export"),
        );
        assert!(
            report.passed
                || (report.message.starts_with("FAIL ")
                    && report.message.contains("first difference")),
            "{contract_path}: OBJ exact diff report was not actionable: {}",
            report.message
        );
        assert_ne!(
            stable_test_hash64(generated.as_bytes()),
            0,
            "{contract_path}: generated OBJ hash"
        );
        checked += 1;
    }

    assert!(
        checked > 0,
        "expected at least one Mol* OBJ reference fixture"
    );
}

#[test]
fn molfig_stl_outputs_are_diffed_against_all_molstar_reference_fixtures() {
    let manifest = include_str!("../../tests/expected/molstar-reference/reference-fixtures.txt");
    let mut checked = 0;

    for contract_path in reference_manifest_entries_for_format(manifest, "stl") {
        let contract = read_manifest_file(contract_path);
        let reference = read_manifest_file_bytes(contract_value(&contract, "stl_reference"));
        let reference_stats = stl_stats(&reference);
        assert_eq!(
            reference_stats.byte_len,
            contract_usize(&contract, "stl_byte_len"),
            "{contract_path}: reference byte_len"
        );
        assert_eq!(
            reference_stats.facet_count,
            contract_usize(&contract, "stl_facet_count"),
            "{contract_path}: reference facet_count"
        );

        let fixture = read_manifest_file_bytes(contract_value(&contract, "fixture"));
        let options = contract_value(&contract, "options");
        let generated = convert_to_stl(&fixture, options.as_bytes()).unwrap();
        let actual_stats = stl_stats(&generated);

        assert!(actual_stats.facet_count > 0, "{contract_path}: STL facets");
        assert_ratio_in_range(
            actual_stats.facet_count,
            reference_stats.facet_count,
            0.25,
            1.25,
            &format!("{contract_path}: STL facet count"),
        );
        let report = crate::diff_bytes(
            &reference,
            &generated,
            &format!("{contract_path}: STL exact export"),
        );
        assert!(
            report.passed
                || (report.message.starts_with("FAIL ")
                    && report.message.contains("first difference")),
            "{contract_path}: STL exact diff report was not actionable: {}",
            report.message
        );
        assert_ne!(
            stable_test_hash64(&generated),
            0,
            "{contract_path}: generated STL hash"
        );
        checked += 1;
    }

    assert!(
        checked > 0,
        "expected at least one Mol* STL reference fixture"
    );
}

#[test]
fn molstar_reference_obj_stl_artifacts_are_cross_checked_at_sparse_slots() {
    let manifest = include_str!("../../tests/expected/molstar-reference/reference-fixtures.txt");
    let mut checked = 0;

    for contract_path in reference_manifest_entries_for_formats(manifest, &["obj", "stl"]) {
        let contract = read_manifest_file(contract_path);
        let reference_obj = read_manifest_file(contract_value(&contract, "obj_reference"));
        let reference_stl = read_manifest_file_bytes(contract_value(&contract, "stl_reference"));
        let obj_vertices = obj_vectors(&reference_obj, "v ");
        let obj_faces = obj_face_indices(&reference_obj);
        let obj_stats = obj_stats(&reference_obj);
        let stl_stats = stl_stats(&reference_stl);

        assert_eq!(
            obj_vertices.len(),
            obj_stats.vertex_count,
            "{contract_path}"
        );
        assert_eq!(obj_faces.len(), obj_stats.face_count, "{contract_path}");
        assert_eq!(
            stl_stats.facet_count,
            obj_stats.face_count * 3,
            "{contract_path}: Mol* STL stores drawCount sparse facet slots"
        );

        let drift = obj_stl_sparse_slot_drift(&obj_vertices, &obj_faces, &reference_stl);

        assert!(
            drift.rounded_mismatch_count > 0,
            "{contract_path}: reference STL raw coordinates now round back to the reference OBJ at every sparse slot; update this diagnostic"
        );
        assert!(
            drift.max_delta > 0.0005,
            "{contract_path}: reference OBJ/STL are now raw-coordinate identical across all sparse slots; update this diagnostic"
        );
        assert!(
            drift.max_delta < 0.0025,
            "{contract_path}: reference OBJ/STL sparse slots drifted beyond the pinned full-scan envelope: {drift:?}"
        );
        assert_eq!(
            drift.total_components,
            contract_usize(&contract, "obj_stl_sparse_slot_total_components"),
            "{contract_path}: obj_stl_sparse_slot_total_components"
        );
        assert_eq!(
            drift.total_components,
            obj_faces.len() * 9,
            "{contract_path}: full sparse-slot component scan coverage"
        );
        assert_eq!(
            drift.rounded_mismatch_count,
            contract_usize(&contract, "obj_stl_sparse_slot_rounding_mismatch_count"),
            "{contract_path}: obj_stl_sparse_slot_rounding_mismatch_count"
        );
        assert!(
            (drift.max_delta - contract_f32(&contract, "obj_stl_sparse_slot_max_delta")).abs()
                <= 0.000_000_1,
            "{contract_path}: obj_stl_sparse_slot_max_delta drift changed: {drift:?}"
        );
        assert!(
            drift.rounded_mismatch_count > obj_faces.len(),
            "{contract_path}: expected full-scan reference OBJ/STL raw-coordinate drift to affect more than one axis per face on average: {drift:?}"
        );
        let first = drift
            .first_rounded_mismatch
            .expect("drift count is nonzero");
        assert_eq!(
            (first.face_index, first.stl_facet_index, first.vertex_slot),
            (0, 0, 0),
            "{contract_path}: first reference OBJ/STL sparse-slot drift location changed: {drift:?}"
        );
        assert_eq!(
            [
                first.face_index,
                first.stl_facet_index,
                first.vertex_slot,
                first.axis
            ],
            contract_usize_array4(&contract, "obj_stl_sparse_slot_first_mismatch"),
            "{contract_path}: obj_stl_sparse_slot_first_mismatch"
        );
        assert_eq!(
            first.axis, 2,
            "{contract_path}: first reference OBJ/STL sparse-slot drift axis changed: {drift:?}"
        );
        assert!(
            (first.rounded_stl_value - first.obj_value).abs() > 0.000_01,
            "{contract_path}: first drift component no longer fails OBJ rounding: {drift:?}"
        );
        let max_component = drift
            .max_delta_component
            .expect("max delta should be present for a non-empty OBJ/STL reference");
        assert_eq!(
            max_component.stl_facet_index,
            max_component.face_index * 3,
            "{contract_path}: max drift component should map OBJ face to sparse STL slot"
        );
        assert!(max_component.axis < 3, "{contract_path}: {drift:?}");
        assert!(
            ((max_component.obj_value - max_component.stl_value).abs() - drift.max_delta).abs()
                <= f32::EPSILON,
            "{contract_path}: max drift component must explain max_delta: {drift:?}"
        );
        assert!(
            (max_component.rounded_stl_value - max_component.obj_value).abs() > 0.000_01,
            "{contract_path}: max drift component should cross OBJ rounding: {drift:?}"
        );
        checked += 1;
    }

    assert!(
        checked > 0,
        "expected at least one Mol* OBJ/STL reference fixture"
    );
}

#[test]
fn molfig_generated_obj_stl_exports_share_sparse_slots_and_obj_rounding() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A 1 1 0.1254 0.3333 0.7777\nATOM 2 O O ALA A 1 1 1.9876 0.4321 -0.3333\n#\n";
    let options =
        br#"{"format":"cif","representation":"spacefill","sphere-detail":1,"center":true}"#;
    let obj = String::from_utf8(convert_to_obj(cif, options).unwrap()).unwrap();
    let stl = convert_to_stl(cif, options).unwrap();
    let obj_vertices = obj_vectors(&obj, "v ");
    let obj_faces = obj_face_indices(&obj);
    let stl_stats = stl_stats(&stl);

    assert_eq!(
        stl_stats.facet_count,
        obj_faces.len() * 3,
        "molfig STL must keep Mol* drawCount sparse facet slots"
    );

    for (face_index, face) in obj_faces.iter().enumerate() {
        let stl_facet_index = face_index * 3;
        for vertex_slot in 0..3 {
            let obj_vertex = obj_vertices[face[vertex_slot]];
            let stl_vertex = stl_facet_vertex(&stl, stl_facet_index, vertex_slot);
            for axis in 0..3 {
                assert!(
                    (molstar_obj_rounded_coordinate(stl_vertex[axis]) - obj_vertex[axis]).abs()
                        <= 0.000_01,
                    "face {face_index} vertex {vertex_slot} axis {axis}: STL raw {:?} did not round back to OBJ {:?}",
                    stl_vertex,
                    obj_vertex
                );
            }
        }
    }
}

fn molstar_obj_rounded_coordinate(value: f32) -> f32 {
    let rounded = ((value as f64 * 1000.0) + 0.5).floor() / 1000.0;
    if rounded == 0.0 {
        0.0
    } else {
        rounded as f32
    }
}

#[derive(Clone, Copy, Debug)]
struct ObjStlSparseSlotDrift {
    total_components: usize,
    rounded_mismatch_count: usize,
    first_rounded_mismatch: Option<ObjStlSparseSlotDriftComponent>,
    max_delta: f32,
    max_delta_component: Option<ObjStlSparseSlotDriftComponent>,
}

#[derive(Clone, Copy, Debug)]
struct ObjStlSparseSlotDriftComponent {
    face_index: usize,
    stl_facet_index: usize,
    vertex_slot: usize,
    axis: usize,
    obj_value: f32,
    stl_value: f32,
    rounded_stl_value: f32,
}

fn obj_stl_sparse_slot_drift(
    obj_vertices: &[[f32; 3]],
    obj_faces: &[[usize; 3]],
    stl: &[u8],
) -> ObjStlSparseSlotDrift {
    let mut drift = ObjStlSparseSlotDrift {
        total_components: 0,
        rounded_mismatch_count: 0,
        first_rounded_mismatch: None,
        max_delta: 0.0,
        max_delta_component: None,
    };
    for (face_index, face) in obj_faces.iter().enumerate() {
        let stl_facet_index = face_index * 3;
        for vertex_slot in 0..3 {
            let obj_vertex = obj_vertices[face[vertex_slot]];
            let stl_vertex = stl_facet_vertex(stl, stl_facet_index, vertex_slot);
            for axis in 0..3 {
                drift.total_components += 1;
                let obj_value = obj_vertex[axis];
                let stl_value = stl_vertex[axis];
                let rounded_stl_value = molstar_obj_rounded_coordinate(stl_value);
                let component = ObjStlSparseSlotDriftComponent {
                    face_index,
                    stl_facet_index,
                    vertex_slot,
                    axis,
                    obj_value,
                    stl_value,
                    rounded_stl_value,
                };
                let delta = (obj_value - stl_value).abs();
                if delta > drift.max_delta {
                    drift.max_delta = delta;
                    drift.max_delta_component = Some(component);
                }
                if (rounded_stl_value - obj_value).abs() > 0.000_01 {
                    drift.rounded_mismatch_count += 1;
                    if drift.first_rounded_mismatch.is_none() {
                        drift.first_rounded_mismatch = Some(component);
                    }
                }
            }
        }
    }
    drift
}

#[test]
fn nine_r1o_asymmetric_ply_vs_assembly_one_reference_gap_matches_snapshot() {
    let Some(actual) = nine_r1o_asymmetric_ply_vs_assembly_one_reference_gap_json() else {
        eprintln!("skipping 9R1O reference gap snapshot; reference OBJ/STL artifacts are absent");
        return;
    };
    if std::env::var("MOLFIG_PRINT_9R1O_GEOMETRY_SNAPSHOT").as_deref() == Ok("1") {
        eprintln!("{actual}");
    }

    assert!(actual.contains(&format!(
        r#""molstar_reference_commit": "{}""#,
        MOLSTAR_REFERENCE_COMMIT
    )));
    assert!(
        actual.contains(r#""name": "9r1o-asymmetric-ply-vs-assembly-1-reference-gap""#),
        "{actual}"
    );
    assert!(
        actual.contains(
            r#""comparison_scope": "diagnostic only: package-owned asymmetric-unit PLY contract vs pinned Mol* assembly 1 OBJ/STL reference""#
        ),
        "{actual}"
    );
    assert!(
        actual.contains(r#""obj_vertex_count_delta": -1160"#),
        "{actual}"
    );
    assert!(
        actual.contains(r#""obj_face_count_delta": -12736"#),
        "{actual}"
    );
    assert!(
        actual.contains(r#""stl_facet_count_delta": -38208"#),
        "{actual}"
    );
    assert!(
        actual.contains(
            r#""face_samples":[{"label":"start","index":0,"indices":[1,21,0],"group":0}"#
        ),
        "{actual}"
    );
    assert!(
        actual.contains(r#""obj_vertex_sample_abs_delta_vs_ply": [{"label":"start","molfig_index":0,"reference_index":0,"abs_delta":[13.9143,5.8752,2.087]}"#),
        "{actual}"
    );
    assert!(
        actual.contains(
            r#""ply_reference": "Mol* geo-export has no PLY exporter; molfig PLY is package-owned""#
        ),
        "{actual}"
    );
}

#[test]
fn molfig_ply_contract_outputs_match_all_reference_fixtures() {
    let manifest = include_str!("../../tests/expected/ply/reference-fixtures.txt");
    let mut checked = 0;

    for contract_path in manifest
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
    {
        let contract = read_manifest_file(contract_path);
        let fixture_path = contract_value(&contract, "fixture");
        let options = contract_value(&contract, "options");
        let fixture = read_manifest_file_bytes(fixture_path);

        let ply = String::from_utf8(convert_to_ply(&fixture, options.as_bytes()).unwrap())
            .unwrap_or_else(|_| panic!("{contract_path}: PLY output is not valid UTF-8"));
        assert_eq!(
            contract_value(&contract, "format"),
            "ascii-ply",
            "{contract_path}: format"
        );
        assert_eq!(
            ply_header_count(&ply, "element vertex "),
            contract_usize(&contract, "vertex_count"),
            "{contract_path}: vertex_count"
        );
        assert_eq!(
            ply_header_count(&ply, "element face "),
            contract_usize(&contract, "face_count"),
            "{contract_path}: face_count"
        );
        assert_eq!(
            ply_header_comment_usize(&ply, "comment molfig_group_count "),
            contract_usize(&contract, "group_count"),
            "{contract_path}: group_count"
        );
        assert_eq!(
            ply.len(),
            contract_usize(&contract, "byte_len"),
            "{contract_path}: byte_len"
        );
        assert_eq!(
            format!("{:016x}", stable_test_hash64(ply.as_bytes())),
            contract_value(&contract, "fnv1a64"),
            "{contract_path}: fnv1a64"
        );
        assert!(
            ply_vertices_are_finite(&ply),
            "{contract_path}: non-finite vertex"
        );
        assert!(
            ply_faces_have_valid_indices(&ply),
            "{contract_path}: out-of-range face index"
        );
        if let Some(expected) = contract_optional_value(&contract, "vertex_samples") {
            assert_eq!(
                ply_vertex_samples_json(&ply),
                expected,
                "{contract_path}: vertex_samples"
            );
        }
        if let Some(expected) = contract_optional_value(&contract, "face_samples") {
            assert_eq!(
                ply_face_samples_json(&ply),
                expected,
                "{contract_path}: face_samples"
            );
        }
        checked += 1;
    }

    assert!(checked > 0, "expected at least one PLY reference fixture");
}

#[test]
fn molstar_reference_artifact_contracts_match_all_reference_fixtures() {
    let manifest = include_str!("../../tests/expected/molstar-reference/reference-fixtures.txt");
    let mut checked = 0;

    for entry in reference_manifest_entries_with_formats(manifest) {
        let contract_path = entry.contract_path;
        let contract = read_manifest_file(contract_path);
        assert_eq!(
            contract_value(&contract, "molstar_reference_commit"),
            MOLSTAR_REFERENCE_COMMIT,
            "{contract_path}: Mol* commit"
        );
        assert!(
            !read_manifest_file_bytes(contract_value(&contract, "fixture")).is_empty(),
            "{contract_path}: fixture"
        );
        MeshOptions::from_json(contract_value(&contract, "options").as_bytes())
            .unwrap_or_else(|error| panic!("{contract_path}: invalid options JSON: {error}"));

        if entry.includes_format("obj") {
            let obj_path = contract_value(&contract, "obj_reference");
            let obj = read_manifest_file(obj_path);
            let obj_stats = obj_stats(&obj);
            assert_eq!(
                obj.len(),
                contract_usize(&contract, "obj_byte_len"),
                "{contract_path}: obj_byte_len"
            );
            assert_eq!(
                format!("{:016x}", stable_test_hash64(obj.as_bytes())),
                contract_value(&contract, "obj_fnv1a64"),
                "{contract_path}: obj_fnv1a64"
            );
            assert_eq!(
                obj_stats.vertex_count,
                contract_usize(&contract, "obj_vertex_count"),
                "{contract_path}: obj_vertex_count"
            );
            assert_eq!(
                obj_stats.normal_count,
                contract_usize(&contract, "obj_normal_count"),
                "{contract_path}: obj_normal_count"
            );
            assert_eq!(
                obj_stats.face_count,
                contract_usize(&contract, "obj_face_count"),
                "{contract_path}: obj_face_count"
            );
            assert_f32_array_close(
                obj_stats.min,
                contract_f32_array(&contract, "obj_bounds_min"),
                0.0001,
                &format!("{contract_path}: obj_bounds_min"),
            );
            assert_f32_array_close(
                obj_stats.max,
                contract_f32_array(&contract, "obj_bounds_max"),
                0.0001,
                &format!("{contract_path}: obj_bounds_max"),
            );
        }

        if entry.includes_format("stl") {
            let stl_path = contract_value(&contract, "stl_reference");
            let stl = read_manifest_file_bytes(stl_path);
            let stl_stats = stl_stats(&stl);
            assert_eq!(
                stl_stats.byte_len,
                contract_usize(&contract, "stl_byte_len"),
                "{contract_path}: stl_byte_len"
            );
            assert_eq!(
                format!("{:016x}", stable_test_hash64(&stl)),
                contract_value(&contract, "stl_fnv1a64"),
                "{contract_path}: stl_fnv1a64"
            );
            assert_eq!(
                stl_stats.facet_count,
                contract_usize(&contract, "stl_facet_count"),
                "{contract_path}: stl_facet_count"
            );
            assert_f32_array_close(
                stl_stats.min,
                contract_f32_array(&contract, "stl_bounds_min"),
                0.0001,
                &format!("{contract_path}: stl_bounds_min"),
            );
            assert_f32_array_close(
                stl_stats.max,
                contract_f32_array(&contract, "stl_bounds_max"),
                0.0001,
                &format!("{contract_path}: stl_bounds_max"),
            );
        }

        if entry.includes_format("json") {
            let json_path = contract_value(&contract, "json_reference");
            let json = read_manifest_file(json_path);
            assert_eq!(
                json.len(),
                contract_usize(&contract, "json_byte_len"),
                "{contract_path}: json_byte_len"
            );
            assert_eq!(
                format!("{:016x}", stable_test_hash64(json.as_bytes())),
                contract_value(&contract, "json_fnv1a64"),
                "{contract_path}: json_fnv1a64"
            );
            assert!(
                json.contains(contract_value(&contract, "fixture")),
                "{contract_path}: JSON summary should name fixture"
            );
        }
        checked += 1;
    }

    assert!(
        checked > 0,
        "expected at least one Mol* reference fixture contract"
    );
}

#[test]
fn molstar_reference_json_summaries_diff_all_reference_fixtures() {
    let manifest = include_str!("../../tests/expected/molstar-reference/reference-fixtures.txt");
    let mut checked = 0;

    for contract_path in reference_manifest_entries_for_format(manifest, "json") {
        let contract = read_manifest_file(contract_path);
        let reference = read_manifest_file(contract_value(&contract, "json_reference"));
        let generated = molstar_reference_summary_from_contract(contract_path, &contract);
        let report = crate::diff_text(
            &reference,
            &generated,
            &format!("{contract_path}: json summary"),
        );
        assert!(report.passed, "{}", report.message);
        checked += 1;
    }

    assert!(
        checked > 0,
        "expected at least one Mol* reference JSON summary"
    );
}

#[test]
fn molstar_reference_manifest_accepts_extended_contract_entries() {
    let manifest = "\
# comments are ignored
	tests/expected/molstar-reference/9r1o-molstar-assembly-1.reference.contract
	formats=json,obj contract=tests/expected/molstar-reference/9r1o-molstar-assembly-1.reference.contract # inline comment
	path=tests/expected/molstar-reference/9r1o-molstar-assembly-1.reference.contract formats=stl tag=binary
	";

    assert_eq!(
        reference_manifest_entries(manifest),
        vec![
            "tests/expected/molstar-reference/9r1o-molstar-assembly-1.reference.contract",
            "tests/expected/molstar-reference/9r1o-molstar-assembly-1.reference.contract",
            "tests/expected/molstar-reference/9r1o-molstar-assembly-1.reference.contract"
        ]
    );
    assert_eq!(
        reference_manifest_entries_for_format(manifest, "json").len(),
        2
    );
    assert_eq!(
        reference_manifest_entries_for_format(manifest, "obj").len(),
        2
    );
    assert_eq!(
        reference_manifest_entries_for_format(manifest, "stl").len(),
        2
    );
    assert_eq!(
        reference_manifest_entries_for_formats(manifest, &["obj", "stl"]).len(),
        1
    );
    assert_eq!(
        reference_manifest_entries_with_formats(manifest)[1],
        ReferenceManifestEntry {
            contract_path:
                "tests/expected/molstar-reference/9r1o-molstar-assembly-1.reference.contract",
            formats: Some(vec!["json", "obj"])
        }
    );
}

#[test]
fn alt_loc_selection_remaps_bonds() {
    let pdb = b"ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 10.00           N\nATOM      2  CA BGLY A   1       0.500   0.000   0.000  1.00 10.00           C\nATOM      3  C   GLY A   1       1.000   0.000   0.000  1.00 10.00           C\nATOM      4  O   GLY A   1       2.000   0.000   0.000  1.00 10.00           O\nCONECT    1    3\nEND\n";
    let mol = parse_molecule_with_options(
        pdb,
        &MeshOptions {
            format: InputFormat::Pdb,
            assembly: None,
            alt_loc: "A".to_string(),
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(
        mol.atoms.iter().map(|a| a.id).collect::<Vec<_>>(),
        vec![1, 3, 4]
    );
    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.atoms[mol.bonds[0].a].id, 1);
    assert_eq!(mol.atoms[mol.bonds[0].b].id, 3);
    assert_eq!(mol.bond_metadata.len(), 1);
    assert_eq!(mol.bond_metadata[0].source, BondSource::PdbConect);
    assert_eq!(mol.bond_metadata[0].order, 1);
    assert!(mol.bond_metadata[0].flags.contains(BondFlags::COVALENT));
}

#[test]
fn highest_occupancy_alt_loc_selects_best_site() {
    let pdb = b"ATOM      1  CA AGLY A   1       0.000   0.000   0.000  0.25 10.00           C\nATOM      2  CA BGLY A   1       1.000   0.000   0.000  0.75 10.00           C\nATOM      3  N   GLY A   1       2.000   0.000   0.000  1.00 10.00           N\nEND\n";
    let mol = parse_molecule_with_options(
        pdb,
        &MeshOptions {
            format: InputFormat::Pdb,
            assembly: None,
            alt_loc: "highest-occupancy".to_string(),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.atoms.len(), 2);
    assert_eq!(mol.atoms[0].id, 2);
    assert_eq!(mol.atoms[0].alt_id, "B");
    assert_eq!(mol.atoms[0].occupancy, 0.75);
}

#[test]
fn atomic_protein_no_altloc_fixture_has_no_alternate_sites() {
    let fixtures: &[(&[u8], InputFormat)] = &[
        (
            include_bytes!("../../tests/fixtures/pdb/atomic-protein-no-altloc.pdb"),
            InputFormat::Pdb,
        ),
        (
            include_bytes!("../../tests/fixtures/cif/atomic-protein-no-altloc.cif"),
            InputFormat::Cif,
        ),
    ];

    for (bytes, format) in fixtures {
        let mol = parse_molecule_with_options(
            bytes,
            &MeshOptions {
                format: *format,
                assembly: None,
                alt_loc: "all".to_string(),
                infer_bonds: false,
                ..MeshOptions::default()
            },
        )
        .unwrap();

        assert_eq!(mol.atoms.len(), 7);
        assert_eq!(mol.bonds.len(), 6);
        assert!(mol.atoms.iter().all(|atom| atom.alt_id.is_empty()));
        assert_eq!(mol.atomic_structure().alt_loc_count(), 0);
    }
}

#[test]
fn atomic_protein_altloc_tie_fixture_prefers_a_on_highest_occupancy_tie() {
    let fixtures: &[(&[u8], InputFormat)] = &[
        (
            include_bytes!("../../tests/fixtures/pdb/atomic-protein-altloc-tie.pdb"),
            InputFormat::Pdb,
        ),
        (
            include_bytes!("../../tests/fixtures/cif/atomic-protein-altloc-tie.cif"),
            InputFormat::Cif,
        ),
    ];

    for (bytes, format) in fixtures {
        let all = parse_molecule_with_options(
            bytes,
            &MeshOptions {
                format: *format,
                assembly: None,
                alt_loc: "all".to_string(),
                infer_bonds: false,
                ..MeshOptions::default()
            },
        )
        .unwrap();
        assert_eq!(all.atoms.len(), 7);
        assert_eq!(
            all.atoms
                .iter()
                .filter(|atom| atom.name == "CB")
                .map(|atom| (atom.alt_id.as_str(), atom.occupancy))
                .collect::<Vec<_>>(),
            vec![("A", 0.5), ("B", 0.5)]
        );
        assert_eq!(all.atomic_structure().alt_loc_count(), 2);

        let selected = parse_molecule_with_options(
            bytes,
            &MeshOptions {
                format: *format,
                assembly: None,
                alt_loc: "highest-occupancy".to_string(),
                infer_bonds: false,
                ..MeshOptions::default()
            },
        )
        .unwrap();

        assert_eq!(selected.atoms.len(), 6);
        let cb = selected
            .atoms
            .iter()
            .find(|atom| atom.name == "CB")
            .unwrap();
        assert_eq!(cb.alt_id, "A");
        assert_eq!(cb.occupancy, 0.5);
        assert_eq!(cb.position.z, 1.0);
        assert!(selected.atoms.iter().all(|atom| atom.alt_id != "B"));
    }
}

#[test]
fn computed_bond_flags_match_molstar_metal_and_covalent_cases() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 FE FE HEM A 1 0.000 0.000 0.000\nHETATM 2 O O1 HEM A 1 2.000 0.000 0.000\nATOM 3 C C1 LIG B 1 10.000 0.000 0.000\nATOM 4 C C2 LIG B 1 11.200 0.000 0.000\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: true,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bonds.len(), 2);
    let metadata_for_pair = |a: &str, b: &str| {
        mol.bonds
            .iter()
            .enumerate()
            .find_map(|(index, bond)| {
                let names = [&mol.atoms[bond.a].name, &mol.atoms[bond.b].name];
                ((names[0] == a && names[1] == b) || (names[0] == b && names[1] == a))
                    .then_some(&mol.bond_metadata[index])
            })
            .unwrap()
    };

    let metal = metadata_for_pair("FE", "O1");
    assert_eq!(metal.source, BondSource::Computed);
    assert!(metal.flags.contains(BondFlags::COMPUTED));
    assert!(metal.flags.contains(BondFlags::METALLIC_COORDINATION));
    assert!(!metal.flags.contains(BondFlags::COVALENT));

    let covalent = metadata_for_pair("C1", "C2");
    assert_eq!(covalent.source, BondSource::Computed);
    assert!(covalent.flags.contains(BondFlags::COMPUTED));
    assert!(covalent.flags.contains(BondFlags::COVALENT));
    assert!(!covalent.flags.contains(BondFlags::METALLIC_COORDINATION));
}

#[test]
fn cif_struct_conn_creates_explicit_bonds_with_metadata() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\nATOM 1 C C1 LIG A 1 0.000 0.000 0.000 1.00\nATOM 2 O O1 LIG A 1 9.000 0.000 0.000 1.00\n#\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_symmetry\n_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_symmetry\n_struct_conn.pdbx_value_order\n_struct_conn.pdbx_dist_value\nmetal1 metalc LIG A 1 C1 1_555 LIG A 1 O1 2_666 doub 2.10\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.bonds[0].a, 0);
    assert_eq!(mol.bonds[0].b, 1);
    assert_eq!(mol.bond_metadata.len(), 1);
    let metadata = &mol.bond_metadata[0];
    assert_eq!(metadata.source, BondSource::StructConn);
    assert_eq!(metadata.order, 2);
    assert!(metadata.flags.contains(BondFlags::METALLIC_COORDINATION));
    assert_eq!(metadata.key, 0);
    assert_eq!(metadata.distance, Some(2.10));
    let struct_conn = metadata.struct_conn.as_ref().unwrap();
    assert_eq!(struct_conn.id, "metal1");
    assert_eq!(struct_conn.conn_type_id, "metalc");
    assert_eq!(struct_conn.value_order, "doub");
    assert_eq!(struct_conn.partner_a_symmetry, "1_555");
    assert_eq!(struct_conn.partner_b_symmetry, "2_666");
}

#[test]
fn cif_struct_conn_survives_alt_loc_b_selection_without_alt_ids() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\nATOM 1 C C1 A LIG A 1 0.000 0.000 0.000 0.40\nATOM 2 C C1 B LIG A 1 1.000 0.000 0.000 0.60\nATOM 3 O O1 . LIG A 1 9.000 0.000 0.000 1.00\n#\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_atom_id\n1 covale A 1 C1 A 1 O1\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            alt_loc: "B".to_string(),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(
        mol.atoms.iter().map(|a| a.id).collect::<Vec<_>>(),
        vec![2, 3]
    );
    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.atoms[mol.bonds[0].a].id, 2);
    assert_eq!(mol.atoms[mol.bonds[0].b].id, 3);
    assert_eq!(mol.bond_metadata[0].source, BondSource::StructConn);
    assert_eq!(mol.bond_metadata[0].key, 0);
}

#[test]
fn cif_struct_conn_auth_only_matches_auth_atom_site_ids() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.auth_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C C1 X1 LIG L A 10 100 0.000 0.000 0.000\nATOM 2 O O1 Y1 LIG L A 10 100 8.000 0.000 0.000\n#\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_auth_asym_id\n_struct_conn.ptnr1_auth_seq_id\n_struct_conn.ptnr1_auth_atom_id\n_struct_conn.ptnr2_auth_asym_id\n_struct_conn.ptnr2_auth_seq_id\n_struct_conn.ptnr2_auth_atom_id\n1 covale A 100 X1 A 100 Y1\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.atoms[mol.bonds[0].a].auth_name, "X1");
    assert_eq!(mol.atoms[mol.bonds[0].b].auth_name, "Y1");
    assert_eq!(mol.bond_metadata[0].source, BondSource::StructConn);
}

#[test]
fn cif_struct_conn_label_ids_do_not_match_auth_ids() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.auth_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C C1 AX LIG L A 10 100 0.000 0.000 0.000\nATOM 2 O O1 BY LIG L A 10 100 8.000 0.000 0.000\nATOM 3 N NZ C1 LIG X L 99 10 4.000 0.000 0.000\n#\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_atom_id\nlabel-link covale L 10 C1 L 10 O1\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.atoms[mol.bonds[0].a].id, 1);
    assert_eq!(mol.atoms[mol.bonds[0].b].id, 2);
    assert_eq!(mol.bond_metadata[0].source, BondSource::StructConn);
    assert_eq!(
        mol.bond_metadata[0].struct_conn.as_ref().unwrap().id,
        "label-link"
    );
}

#[test]
fn cif_struct_conn_label_path_prefers_auth_seq_id_like_molstar() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.auth_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C C1 AX LIG L A 10 100 0.000 0.000 0.000\nATOM 2 O O1 BY LIG L A 10 100 8.000 0.000 0.000\n#\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_auth_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_auth_seq_id\n_struct_conn.ptnr2_label_atom_id\nseq-link covale L 999 100 C1 L 999 100 O1\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.atoms[mol.bonds[0].a].id, 1);
    assert_eq!(mol.atoms[mol.bonds[0].b].id, 2);
}

#[test]
fn cif_struct_conn_empty_insertion_code_is_not_wildcard() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C C1 LIG A 1 A 0.000 0.000 0.000\nATOM 2 O O1 LIG A 1 A 8.000 0.000 0.000\n#\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_atom_id\nempty-ins covale A 1 C1 A 1 O1\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert!(mol.bonds.is_empty());
}

#[test]
fn cif_struct_conn_fallback_order_insertion_code_and_symmetry_affect_unit_bonds() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 N NZ LYS A 1 A 0.000 0.000 0.000\nATOM 2 N NZ LYS A 1 B 1.000 0.000 0.000\nHETATM 3 C C15 RET A 2 . 2.000 0.000 0.000\n#\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.pdbx_ptnr1_PDB_ins_code\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_symmetry\n_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n_struct_conn.pdbx_ptnr2_PDB_ins_code\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_symmetry\nret-link covale LYS A 1 B NZ 1_555 RET A 2 . C15 2_666\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.atoms[mol.bonds[0].a].insertion_code, "B");
    assert_eq!(mol.bond_metadata[0].order, 2);
    assert!(mol.bond_metadata[0].flags.contains(BondFlags::COVALENT));
    let structure = mol.atomic_structure();
    assert_eq!(structure.intra_unit_bond_count, 0);
    assert_eq!(structure.inter_unit_bonds.len(), 0);
}

#[test]
fn cif_struct_conn_inter_unit_matching_uses_units_and_distance_not_symmetry_names() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C C1 LIG A 1 0.000 0.000 0.000\nATOM 2 O O1 LIG B 1 -9.000 0.000 0.000\n#\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_symmetry\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_symmetry\nsym-link covale A 1 C1 1_555 B 1 O1 2_666\n#\nloop_\n_pdbx_struct_assembly.id\n_pdbx_struct_assembly.details\n_pdbx_struct_assembly.oligomeric_details\n_pdbx_struct_assembly.oligomeric_count\n1 'author defined assembly' dimer 2\n#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n1 1 A\n1 2 B\n#\nloop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\n1 1 0 0 0 0 1 0 0 0 0 1 0\n2 1 0 0 10 0 1 0 0 0 0 1 0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: Some("1".to_string()),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let structure = mol.atomic_structure();
    assert_eq!(structure.units.len(), 2);
    assert_eq!(structure.intra_unit_bond_count, 0);
    assert_eq!(structure.inter_unit_bonds.len(), 1);
    let inter = &structure.inter_unit_bonds[0];
    assert_eq!(inter.unit_a, 0);
    assert_eq!(inter.unit_b, 1);
    assert_eq!(inter.order, 1);
    assert!(inter.flags.contains(BondFlags::COVALENT));
    assert_eq!(inter.key, 0);
}

#[test]
fn cif_struct_conn_inter_unit_matching_respects_molstar_default_max_radius() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C C1 LIG A 1 0.000 0.000 0.000\nATOM 2 O O1 LIG B 1 -5.000 0.000 0.000\n#\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_symmetry\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_symmetry\nsym-link covale A 1 C1 1_555 B 1 O1 2_666\n#\nloop_\n_pdbx_struct_assembly.id\n_pdbx_struct_assembly.details\n_pdbx_struct_assembly.oligomeric_details\n_pdbx_struct_assembly.oligomeric_count\n1 'author defined assembly' dimer 2\n#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n1 1 A\n1 2 B\n#\nloop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\n1 1 0 0 0 0 1 0 0 0 0 1 0\n2 1 0 0 10 0 1 0 0 0 0 1 0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: Some("1".to_string()),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let structure = mol.atomic_structure();
    assert_eq!(structure.units.len(), 2);
    assert!(structure.inter_unit_bonds.is_empty());
}

#[test]
fn cif_struct_conn_keeps_row_orientation_and_duplicate_rows_like_molstar() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C C1 LIG A 1 0.000 0.000 0.000\nATOM 2 O O1 LIG A 1 1.000 0.000 0.000\n#\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_symmetry\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_symmetry\nrow-a covale A 1 O1 1_555 A 1 C1 1_555\nrow-b covale A 1 O1 2_666 A 1 C1 2_666\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bonds.len(), 2);
    assert_eq!(mol.bonds[0].a, 0);
    assert_eq!(mol.bonds[0].b, 1);
    assert_eq!(mol.bond_metadata[0].key, 0);
    assert_eq!(mol.bond_metadata[1].key, 1);
    let first = mol.bond_metadata[0].struct_conn.as_ref().unwrap();
    assert_eq!(first.id, "row-a");
    assert_eq!(first.partner_a_atom_index, 1);
    assert_eq!(first.partner_b_atom_index, 0);
    assert_eq!(first.partner_a_symmetry, "1_555");
    let second = mol.bond_metadata[1].struct_conn.as_ref().unwrap();
    assert_eq!(second.id, "row-b");
    assert_eq!(second.partner_a_atom_index, 1);
    assert_eq!(second.partner_b_atom_index, 0);
    assert_eq!(second.partner_a_symmetry, "2_666");
}

#[test]
fn cif_molstar_bond_site_creates_index_pair_metadata() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 10 C C1 LIG A 1 0.000 0.000 0.000\nATOM 20 C C2 LIG A 1 1.000 0.000 0.000\n#\nloop_\n_molstar_bond_site.atom_id_1\n_molstar_bond_site.atom_id_2\n_molstar_bond_site.value_order\n_molstar_bond_site.type_id\n_molstar_bond_site.key\n_molstar_bond_site.distance\n_molstar_bond_site.operator_a\n_molstar_bond_site.operator_b\n10 20 arom covale 7 2.5 1 1\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.bond_metadata[0].source, BondSource::IndexPair);
    assert_eq!(mol.bond_metadata[0].order, 1);
    assert_eq!(mol.bond_metadata[0].key, -1);
    assert_eq!(mol.bond_metadata[0].distance, None);
    assert_eq!(mol.bond_metadata[0].operator_a, -1);
    assert_eq!(mol.bond_metadata[0].operator_b, -1);
    assert!(mol.bond_metadata[0].flags.contains(BondFlags::COVALENT));
    assert!(mol.bond_metadata[0].flags.contains(BondFlags::AROMATIC));
    assert!(!mol.bond_metadata[0].flags.contains(BondFlags::RESONANCE));
    let index_pairs = mol.index_pair_bonds.as_ref().unwrap();
    assert_eq!(index_pairs.bonds.vertex_count, 2);
    assert_eq!(index_pairs.bonds.edge_count, 1);
    assert_eq!(index_pairs.bonds.offset, vec![0, 1, 2]);
    assert_eq!(index_pairs.bonds.a, vec![0, 1]);
    assert_eq!(index_pairs.bonds.b, vec![1, 0]);
    assert!(index_pairs.contains_bond(0));
    assert_eq!(index_pairs.bonds.props.key, vec![-1, -1]);
    assert_eq!(index_pairs.bonds.props.order, vec![1, 1]);
    assert_eq!(index_pairs.bonds.props.distance, vec![-1.0, -1.0]);
    assert_eq!(
        index_pairs.bonds.props.flag,
        vec![mol.bond_metadata[0].flags; 2]
    );
    assert!(index_pairs.max_distance.is_infinite());
    assert!(index_pairs.cacheable);
    assert!(!index_pairs.has_operators);
    assert!(index_pairs.by_same_operator.is_empty());
}

#[test]
fn cif_molstar_bond_site_imports_molstar_order_and_type_semantics() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C C1 LIG A 1 0.000 0.000 0.000\nATOM 2 C C2 LIG A 1 1.000 0.000 0.000\nATOM 3 S S1 LIG A 1 2.000 0.000 0.000\nATOM 4 O O1 LIG A 1 3.000 0.000 0.000\nATOM 5 FE FE1 HEM A 2 4.000 0.000 0.000\nATOM 6 N N1 HEM A 2 5.000 0.000 0.000\n#\nloop_\n_molstar_bond_site.atom_id_1\n_molstar_bond_site.atom_id_2\n_molstar_bond_site.value_order\n_molstar_bond_site.type_id\n1 2 sing covale\n2 3 doub disulf\n3 4 trip hydrog\n4 5 quad metalc\n5 6 arom .\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bonds.len(), 5);
    assert!(mol
        .bond_metadata
        .iter()
        .all(|metadata| metadata.source == BondSource::IndexPair
            && metadata.key == -1
            && metadata.distance.is_none()
            && metadata.operator_a == -1
            && metadata.operator_b == -1));

    let orders = mol
        .bond_metadata
        .iter()
        .map(|metadata| metadata.order)
        .collect::<Vec<_>>();
    assert_eq!(orders, vec![1, 2, 3, 4, 1]);

    assert!(mol.bond_metadata[0].flags.contains(BondFlags::COVALENT));
    assert!(mol.bond_metadata[1].flags.contains(BondFlags::COVALENT));
    assert!(mol.bond_metadata[1].flags.contains(BondFlags::DISULFIDE));
    assert!(mol.bond_metadata[2]
        .flags
        .contains(BondFlags::HYDROGEN_BOND));
    assert!(!mol.bond_metadata[2].flags.contains(BondFlags::COVALENT));
    assert!(mol.bond_metadata[3]
        .flags
        .contains(BondFlags::METALLIC_COORDINATION));
    assert!(!mol.bond_metadata[3].flags.contains(BondFlags::COVALENT));
    assert!(mol.bond_metadata[4].flags.contains(BondFlags::AROMATIC));
    assert!(!mol.bond_metadata[4].flags.contains(BondFlags::RESONANCE));

    let index_pairs = mol.index_pair_bonds.as_ref().unwrap();
    assert_eq!(index_pairs.bonds.edge_count, 5);
    assert!(index_pairs.max_distance.is_infinite());
    assert!(index_pairs.cacheable);
    assert!(!index_pairs.has_operators);
    assert_eq!(index_pairs.bonds.props.order.len(), 10);
    assert_eq!(index_pairs.bonds.props.key, vec![-1; 10]);
    assert_eq!(index_pairs.bonds.props.distance, vec![-1.0; 10]);
}

#[test]
fn cif_molstar_bond_site_uses_molstar_schema_fields_only() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 10 C C1 LIG A 1 0.000 0.000 0.000\nATOM 20 C C2 LIG A 1 1.000 0.000 0.000\n#\nloop_\n_molstar_bond_site.atom_id_1\n_molstar_bond_site.atom_id_2\n_molstar_bond_site.operator_a\n_molstar_bond_site.operator_b\n20 10 2 1\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.bonds[0].a, 0);
    assert_eq!(mol.bonds[0].b, 1);
    assert_eq!(mol.bond_metadata[0].key, -1);
    assert_eq!(mol.bond_metadata[0].distance, None);
    assert_eq!(mol.bond_metadata[0].operator_a, -1);
    assert_eq!(mol.bond_metadata[0].operator_b, -1);
    let index_pairs = mol.index_pair_bonds.as_ref().unwrap();
    assert_eq!(index_pairs.bonds.props.operator_a, vec![-1, -1]);
    assert_eq!(index_pairs.bonds.props.operator_b, vec![-1, -1]);
}

#[test]
fn cif_molstar_bond_site_delo_is_not_aromatic_for_molstar_parity() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C C1 LIG A 1 0.000 0.000 0.000\nATOM 2 C C2 LIG A 1 1.000 0.000 0.000\nATOM 3 C C3 LIG A 1 2.000 0.000 0.000\n#\nloop_\n_molstar_bond_site.atom_id_1\n_molstar_bond_site.atom_id_2\n_molstar_bond_site.value_order\n_molstar_bond_site.type_id\n1 2 delo covale\n2 3 arom covale\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bonds.len(), 2);
    assert_eq!(mol.bond_metadata[0].order, 1);
    assert!(mol.bond_metadata[0].flags.contains(BondFlags::COVALENT));
    assert!(!mol.bond_metadata[0].flags.contains(BondFlags::AROMATIC));
    assert!(!mol.bond_metadata[0].flags.contains(BondFlags::RESONANCE));
    assert!(mol.bond_metadata[1].flags.contains(BondFlags::AROMATIC));
    assert!(!mol.bond_metadata[1].flags.contains(BondFlags::RESONANCE));
}

#[test]
fn molstar_bond_site_export_entries_sort_deduplicate_and_assign_values_like_molstar() {
    let molecule = Molecule {
        atoms: vec![
            test_atom(2, "C2", "A", 1, vec3(0.0, 0.0, 0.0)),
            test_atom(1, "C1", "A", 1, vec3(1.0, 0.0, 0.0)),
            test_atom(4, "S1", "A", 1, vec3(2.0, 0.0, 0.0)),
            test_atom(3, "O1", "A", 1, vec3(3.0, 0.0, 0.0)),
            test_atom(5, "FE1", "A", 1, vec3(4.0, 0.0, 0.0)),
        ],
        bonds: vec![
            Bond { a: 0, b: 1 },
            Bond { a: 1, b: 0 },
            Bond { a: 1, b: 2 },
            Bond { a: 2, b: 3 },
            Bond { a: 3, b: 4 },
        ],
        bond_metadata: vec![
            BondMetadata {
                source: BondSource::IndexPair,
                order: 2,
                flags: BondFlags::COVALENT,
                key: -1,
                distance: None,
                operator_a: -1,
                operator_b: -1,
                struct_conn: None,
            },
            BondMetadata {
                source: BondSource::IndexPair,
                order: 3,
                flags: BondFlags::COVALENT.union(BondFlags::AROMATIC),
                key: -1,
                distance: None,
                operator_a: -1,
                operator_b: -1,
                struct_conn: None,
            },
            BondMetadata {
                source: BondSource::IndexPair,
                order: 1,
                flags: BondFlags::COVALENT.union(BondFlags::DISULFIDE),
                key: -1,
                distance: None,
                operator_a: -1,
                operator_b: -1,
                struct_conn: None,
            },
            BondMetadata {
                source: BondSource::IndexPair,
                order: 4,
                flags: BondFlags::METALLIC_COORDINATION,
                key: -1,
                distance: None,
                operator_a: -1,
                operator_b: -1,
                struct_conn: None,
            },
            BondMetadata {
                source: BondSource::IndexPair,
                order: 0,
                flags: BondFlags::HYDROGEN_BOND,
                key: -1,
                distance: None,
                operator_a: -1,
                operator_b: -1,
                struct_conn: None,
            },
        ],
        ..Molecule::default()
    };

    let entries = molecule.molstar_bond_site_entries();
    assert_eq!(
        entries,
        vec![
            MolstarBondSiteEntry {
                atom_id_1: 1,
                atom_id_2: 2,
                value_order: Some("doub"),
                type_id: Some("covale"),
            },
            MolstarBondSiteEntry {
                atom_id_1: 1,
                atom_id_2: 4,
                value_order: Some("sing"),
                type_id: Some("disulf"),
            },
            MolstarBondSiteEntry {
                atom_id_1: 3,
                atom_id_2: 4,
                value_order: Some("quad"),
                type_id: Some("metalc"),
            },
            MolstarBondSiteEntry {
                atom_id_1: 3,
                atom_id_2: 5,
                value_order: None,
                type_id: Some("hydrog"),
            },
        ]
    );
}

#[test]
fn molstar_bond_site_export_entries_include_inter_unit_bonds_like_molstar() {
    let assembly = Assembly {
        id: "1".to_string(),
        details: String::new(),
        oligomeric_details: String::new(),
        oligomeric_count: None,
        asym_ids: vec!["A".to_string()],
        transforms: Vec::new(),
        generators: vec![AssemblyGenerator::from_transforms(
            "1",
            vec!["A".to_string()],
            0,
            vec![
                Transform::identity(),
                Transform {
                    m: [
                        [1.0, 0.0, 0.0, 10.0],
                        [0.0, 1.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0, 0.0],
                    ],
                },
            ],
            vec![vec!["1".to_string()], vec!["2".to_string()]],
        )],
    };
    let molecule = Molecule {
        atoms: vec![
            test_atom(1, "C1", "A", 1, vec3(0.0, 0.0, 0.0)),
            test_atom(2, "C2", "A", 1, vec3(1.0, 0.0, 0.0)),
        ],
        bonds: vec![Bond { a: 0, b: 1 }],
        bond_metadata: vec![BondMetadata {
            source: BondSource::IndexPair,
            order: 1,
            flags: BondFlags::COVALENT.union(BondFlags::AROMATIC),
            key: -1,
            distance: None,
            operator_a: 1,
            operator_b: 2,
            struct_conn: None,
        }],
        selected_assembly: Some(assembly),
        ..Molecule::default()
    };

    assert_eq!(molecule.atomic_structure().intra_unit_bond_count, 0);
    assert_eq!(molecule.atomic_structure().inter_unit_bonds.len(), 1);
    assert_eq!(
        molecule.molstar_bond_site_entries(),
        vec![MolstarBondSiteEntry {
            atom_id_1: 1,
            atom_id_2: 2,
            value_order: Some("arom"),
            type_id: Some("covale"),
        }]
    );
}

#[test]
fn struct_conn_sing_and_uppercase_aliases_match_molstar() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 N NZ LYS A 1 0.000 0.000 0.000\nHETATM 2 C C15 RET A 2 1.000 0.000 0.000\n#\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.pdbx_value_order\nret-link COVALE LYS A 1 NZ RET A 2 C15 sing\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.bond_metadata[0].order, 1);
    assert!(mol.bond_metadata[0].flags.contains(BondFlags::COVALENT));
}

#[test]
fn struct_conn_type_ids_assign_molstar_bond_flags() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C C1 LIG A 1 0.000 0.000 0.000\nATOM 2 C C2 LIG A 1 1.000 0.000 0.000\nATOM 3 S S1 LIG A 1 2.000 0.000 0.000\nATOM 4 O O1 LIG A 1 3.000 0.000 0.000\nATOM 5 FE FE1 HEM A 2 4.000 0.000 0.000\n#\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_atom_id\nc1 covale A 1 C1 A 1 C2\nd1 disulf A 1 C2 A 1 S1\nh1 hydrog A 1 S1 A 1 O1\nm1 metalc A 1 O1 A 2 FE1\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.bond_metadata.len(), 4);
    assert!(mol.bond_metadata[0].flags.contains(BondFlags::COVALENT));
    assert!(mol.bond_metadata[1].flags.contains(BondFlags::COVALENT));
    assert!(mol.bond_metadata[1].flags.contains(BondFlags::DISULFIDE));
    assert!(mol.bond_metadata[2]
        .flags
        .contains(BondFlags::HYDROGEN_BOND));
    assert!(mol.bond_metadata[3]
        .flags
        .contains(BondFlags::METALLIC_COORDINATION));
}

#[test]
fn insertion_codes_are_distinct_alt_loc_sites() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 10 A 0.000 0.000 0.000\nATOM 2 C CA GLY A 10 B 1.000 0.000 0.000\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.atoms.len(), 2);
    assert_eq!(mol.atoms[0].insertion_code, "A");
    assert_eq!(mol.atoms[1].insertion_code, "B");
}

#[test]
fn empty_insertion_codes_are_preserved_and_split_from_non_empty_codes() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 10 . 0.000 0.000 0.000\nATOM 2 C CB GLY A 10 ? 1.000 0.000 0.000\nATOM 3 C CA GLY A 10 A 2.000 0.000 0.000\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            alt_loc: "all".to_string(),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(
        mol.atoms
            .iter()
            .map(|atom| atom.insertion_code.as_str())
            .collect::<Vec<_>>(),
        vec!["", "", "A"]
    );
    let structure = mol.atomic_structure();
    assert_eq!(structure.model.hierarchy.residues.len(), 2);
    assert_eq!(
        structure.properties.pdbx_pdb_ins_code,
        vec!["".to_string(), "".to_string(), "A".to_string()]
    );
}

#[test]
fn cif_assembly_keeps_generator_asym_ids_separate() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 0.000 0.000 0.000\nATOM 2 C CA GLY B 1 1.000 0.000 0.000\n#\nloop_\n_pdbx_struct_assembly.id\n_pdbx_struct_assembly.details\n_pdbx_struct_assembly.oligomeric_details\n_pdbx_struct_assembly.oligomeric_count\n1 'author defined assembly' dimer 2\n#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n1 1 A\n1 2 B\n#\nloop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\n1 1 0 0 0 0 1 0 0 0 0 1 0\n2 1 0 0 10 0 1 0 0 0 0 1 0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: Some("1".to_string()),
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.atoms.len(), 2);
    assert!(mol.selected_assembly.is_some());
    let assembly = mol.selected_assembly.as_ref().unwrap();
    assert_eq!(assembly.details, "author defined assembly");
    assert_eq!(assembly.oligomeric_details, "dimer");
    assert_eq!(assembly.oligomeric_count, Some(2));
    assert_eq!(mol.atoms[0].chain, "A");
    assert_eq!(mol.atoms[0].position.x, 0.0);
    assert_eq!(mol.atoms[1].chain, "B");
    assert_eq!(mol.atoms[1].position.x, 1.0);

    let structure = mol.atomic_structure();
    assert_eq!(structure.coordinate_system.name, "1");
    assert_eq!(structure.coordinate_system.instance_id, "1");
    assert_eq!(structure.coordinate_system.assembly_id, "1");
    assert_eq!(structure.coordinate_system.oper_id, 0);
    assert!(structure.coordinate_system.oper_list_ids.is_empty());
    assert!(structure.coordinate_system.is_identity);
    assert_eq!(structure.coordinate_system.suffix, "");
    assert_eq!(structure.model.hierarchy.atoms.len(), 2);
    assert_eq!(structure.units.len(), 2);
    assert_eq!(structure.element_count, 2);
    assert_eq!(structure.units[0].operator.name, "ASM_1");
    assert_eq!(structure.units[0].operator.instance_id, "ASM-1");
    assert_eq!(structure.units[0].operator.assembly_id, "1");
    assert_eq!(structure.units[0].operator.oper_id, 1);
    assert_eq!(structure.units[0].operator.oper_list_ids, vec!["1"]);
    assert!(structure.units[0].operator.is_identity);
    assert_eq!(structure.units[0].operator.suffix, "");
    assert_eq!(structure.units[1].operator.name, "ASM_2");
    assert_eq!(structure.units[1].operator.instance_id, "ASM-2");
    assert_eq!(structure.units[1].operator.assembly_id, "1");
    assert_eq!(structure.units[1].operator.oper_id, 2);
    assert_eq!(structure.units[1].operator.oper_list_ids, vec!["2"]);
    assert!(!structure.units[1].operator.is_identity);
    assert_eq!(structure.units[1].operator.suffix, "_2");
    assert_eq!(structure.position(1, 0).unwrap().x, 11.0);

    let geometry = mol.expanded_for_geometry();
    assert_eq!(geometry.atoms.len(), 2);
    assert_eq!(geometry.atoms[1].chain, "B");
    assert_eq!(geometry.atoms[1].position.x, 11.0);
}

#[test]
fn cif_assembly_generators_store_molstar_operator_metadata() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 0.000 0.000 0.000\nATOM 2 C CA GLY B 1 1.000 0.000 0.000\n#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n1 (X0)(1-2) A\n1 3 B\n#\nloop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\nX0 1 0 0 10 0 1 0 0 0 0 1 0\n1 1 0 0 0 0 1 0 0 0 0 1 0\n2 1 0 0 20 0 1 0 0 0 0 1 0\n3 1 0 0 30 0 1 0 0 0 0 1 0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: Some("1".to_string()),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let assembly = mol.selected_assembly.as_ref().unwrap();

    assert_eq!(assembly.generators.len(), 2);
    assert_eq!(assembly.generators[0].operators[0].name, "ASM_1");
    assert_eq!(assembly.generators[0].operators[0].instance_id, "ASM-X0-1");
    assert_eq!(
        assembly.generators[0].operators[0].oper_list_ids,
        vec!["X0", "1"]
    );
    assert_eq!(assembly.generators[0].operators[1].name, "ASM_2");
    assert_eq!(assembly.generators[0].operators[1].instance_id, "ASM-X0-2");
    assert_eq!(assembly.generators[1].operators[0].name, "ASM_3");
    assert_eq!(assembly.generators[1].operators[0].instance_id, "ASM-3");

    let structure = mol.atomic_structure();
    assert!(structure.units.iter().any(|unit| {
        unit.chain_index == 1
            && unit.operator.name == "ASM_3"
            && unit.operator.instance_id == "ASM-3"
            && unit.operator.oper_list_ids == vec!["3"]
    }));
}

#[test]
fn cif_assembly_operator_expression_ignores_molstar_trailing_bracket_suffix() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 1.000 0.000 0.000\n#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n1 (X0)(1-2)] A\n#\nloop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\nX0 1 0 0 10 0 1 0 0 0 0 1 0\n1 1 0 0 0 0 1 0 0 0 0 1 0\n2 1 0 0 20 0 1 0 0 0 0 1 0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: Some("1".to_string()),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let assembly = mol.selected_assembly.as_ref().unwrap();

    assert_eq!(assembly.generators.len(), 1);
    assert_eq!(assembly.generators[0].operators.len(), 2);
    assert_eq!(
        assembly.generators[0]
            .operators
            .iter()
            .map(|operator| operator.oper_list_ids.clone())
            .collect::<Vec<_>>(),
        vec![
            vec!["X0".to_string(), "1".to_string()],
            vec!["X0".to_string(), "2".to_string()]
        ]
    );
    assert_eq!(assembly.generators[0].operators[0].instance_id, "ASM-X0-1");
    assert_eq!(assembly.generators[0].operators[1].instance_id, "ASM-X0-2");

    let structure = mol.atomic_structure();
    assert_eq!(structure.units.len(), 2);
    assert_eq!(structure.position(0, 0).unwrap().x, 11.0);
    assert_eq!(structure.position(1, 0).unwrap().x, 31.0);
}

#[test]
fn cif_assembly_unit_order_follows_generator_operator_and_base_unit_order() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 0.000 0.000 0.000\nATOM 2 C CA GLY B 1 1.000 0.000 0.000\n#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n1 2 B\n1 1,3 A\n#\nloop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\n1 1 0 0 0 0 1 0 0 0 0 1 0\n2 1 0 0 20 0 1 0 0 0 0 1 0\n3 1 0 0 30 0 1 0 0 0 0 1 0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: Some("1".to_string()),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = mol.atomic_structure();

    assert_eq!(
        structure
            .units
            .iter()
            .map(|unit| {
                (
                    unit.chain_indices.clone(),
                    unit.operator.name.clone(),
                    unit.operator.oper_list_ids.clone(),
                )
            })
            .collect::<Vec<_>>(),
        vec![
            (vec![1], "ASM_1".to_string(), vec!["2".to_string()]),
            (vec![0], "ASM_2".to_string(), vec!["1".to_string()]),
            (vec![0], "ASM_3".to_string(), vec!["3".to_string()]),
        ]
    );
    assert_eq!(structure.position(0, 0).unwrap().x, 21.0);
    assert_eq!(structure.position(1, 0).unwrap().x, 0.0);
    assert_eq!(structure.position(2, 0).unwrap().x, 30.0);
    assert_eq!(structure.symmetry_groups[0].unit_ids, vec![0]);
    assert_eq!(structure.symmetry_groups[0].transform_hash, 84_696_351);
    assert_eq!(structure.symmetry_groups[1].unit_ids, vec![1, 2]);
    assert_eq!(structure.symmetry_groups[1].transform_hash, 3_983_810_698);
}

#[test]
fn assembly_index_pair_operator_metadata_filters_intra_and_inter_unit_bonds() {
    let assembly = Assembly {
        id: "1".to_string(),
        details: String::new(),
        oligomeric_details: String::new(),
        oligomeric_count: None,
        asym_ids: vec!["A".to_string()],
        transforms: Vec::new(),
        generators: vec![AssemblyGenerator::from_transforms(
            "1",
            vec!["A".to_string()],
            0,
            vec![
                Transform::identity(),
                Transform {
                    m: [
                        [1.0, 0.0, 0.0, 10.0],
                        [0.0, 1.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0, 0.0],
                    ],
                },
            ],
            vec![vec!["1".to_string()], vec!["2".to_string()]],
        )],
    };
    let structure = Molecule {
        atoms: vec![
            test_atom(1, "C1", "A", 1, vec3(0.0, 0.0, 0.0)),
            test_atom(2, "C2", "A", 1, vec3(1.0, 0.0, 0.0)),
            test_atom(3, "C3", "A", 1, vec3(2.0, 0.0, 0.0)),
        ],
        bonds: vec![Bond { a: 0, b: 1 }, Bond { a: 1, b: 2 }],
        bond_metadata: vec![
            BondMetadata {
                source: BondSource::IndexPair,
                order: 1,
                flags: BondFlags::COVALENT,
                key: 101,
                distance: None,
                operator_a: 1,
                operator_b: 1,
                struct_conn: None,
            },
            BondMetadata {
                source: BondSource::IndexPair,
                order: 1,
                flags: BondFlags::COVALENT,
                key: 102,
                distance: None,
                operator_a: 1,
                operator_b: 2,
                struct_conn: None,
            },
        ],
        selected_assembly: Some(assembly),
        ..Molecule::default()
    }
    .atomic_structure();

    assert_eq!(structure.units.len(), 2);
    assert_eq!(structure.intra_unit_bond_count, 1);
    assert_eq!(structure.units[0].props.intra_unit_bond_count, 1);
    assert_eq!(structure.units[1].props.intra_unit_bond_count, 0);
    assert_eq!(structure.inter_unit_bonds.len(), 1);
    let bond = &structure.inter_unit_bonds[0];
    assert_eq!(
        (
            bond.unit_a,
            bond.index_a,
            bond.unit_b,
            bond.index_b,
            bond.source_bond,
            bond.key
        ),
        (0, 1, 1, 2, 1, 102)
    );
}

#[test]
fn inter_unit_bond_storage_uses_molstar_inter_unit_graph_edges_and_metadata() {
    let assembly = Assembly {
        id: "1".to_string(),
        details: String::new(),
        oligomeric_details: String::new(),
        oligomeric_count: None,
        asym_ids: vec!["A".to_string()],
        transforms: Vec::new(),
        generators: vec![AssemblyGenerator::from_transforms(
            "1",
            vec!["A".to_string()],
            0,
            vec![
                Transform::identity(),
                Transform {
                    m: [
                        [1.0, 0.0, 0.0, 10.0],
                        [0.0, 1.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0, 0.0],
                    ],
                },
                Transform {
                    m: [
                        [1.0, 0.0, 0.0, 20.0],
                        [0.0, 1.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0, 0.0],
                    ],
                },
            ],
            vec![
                vec!["1".to_string()],
                vec!["2".to_string()],
                vec!["3".to_string()],
            ],
        )],
    };
    let structure = Molecule {
        atoms: vec![
            test_atom(1, "C1", "A", 1, vec3(0.0, 0.0, 0.0)),
            test_atom(2, "C2", "A", 1, vec3(1.0, 0.0, 0.0)),
            test_atom(3, "C3", "A", 1, vec3(2.0, 0.0, 0.0)),
        ],
        bonds: vec![Bond { a: 0, b: 1 }, Bond { a: 0, b: 2 }],
        bond_metadata: vec![
            BondMetadata {
                source: BondSource::IndexPair,
                order: 2,
                flags: BondFlags::COVALENT.union(BondFlags::AROMATIC),
                key: 7,
                distance: None,
                operator_a: 1,
                operator_b: 2,
                struct_conn: None,
            },
            BondMetadata {
                source: BondSource::IndexPair,
                order: 1,
                flags: BondFlags::METALLIC_COORDINATION,
                key: 8,
                distance: None,
                operator_a: 1,
                operator_b: 3,
                struct_conn: None,
            },
        ],
        selected_assembly: Some(assembly),
        ..Molecule::default()
    }
    .atomic_structure();

    assert_eq!(structure.inter_unit_bonds.len(), 2);
    let graph = &structure.inter_unit_bond_graph;
    assert_eq!(graph.edge_count, 4);
    assert_eq!(
        graph.edges,
        vec![
            InterUnitBondEdge {
                unit_a: 0,
                unit_b: 1,
                index_a: 0,
                index_b: 1,
                props: InterUnitBondProps {
                    order: 2,
                    flag: BondFlags::COVALENT.union(BondFlags::AROMATIC),
                    key: 7,
                },
            },
            InterUnitBondEdge {
                unit_a: 0,
                unit_b: 2,
                index_a: 0,
                index_b: 2,
                props: InterUnitBondProps {
                    order: 1,
                    flag: BondFlags::METALLIC_COORDINATION,
                    key: 8,
                },
            },
            InterUnitBondEdge {
                unit_a: 1,
                unit_b: 0,
                index_a: 1,
                index_b: 0,
                props: InterUnitBondProps {
                    order: 2,
                    flag: BondFlags::COVALENT.union(BondFlags::AROMATIC),
                    key: 7,
                },
            },
            InterUnitBondEdge {
                unit_a: 2,
                unit_b: 0,
                index_a: 2,
                index_b: 0,
                props: InterUnitBondProps {
                    order: 1,
                    flag: BondFlags::METALLIC_COORDINATION,
                    key: 8,
                },
            },
        ]
    );
    assert_eq!(graph.get_edge_index(0, 0, 1, 1), Some(0));
    assert_eq!(graph.get_edge_index(0, 0, 2, 2), Some(1));
    assert_eq!(graph.get_edge_index(1, 1, 0, 0), Some(2));
    assert_eq!(graph.get_edge_index(2, 2, 0, 0), Some(3));
    assert!(graph.has_edge(0, 0, 1, 1));
    assert_eq!(graph.get_edge(0, 0, 1, 1), graph.edges.first());
    assert_eq!(graph.get_edge_indices(0, 0), &[0, 1]);
    assert_eq!(graph.get_edge_indices(1, 1), &[2]);
    assert_eq!(graph.get_edge_indices(2, 2), &[3]);

    let connected = graph.get_connected_units(0);
    assert_eq!(connected.len(), 2);
    assert!(connected[0].are_units_ordered());
    assert_eq!(connected[0].unit_a, 0);
    assert_eq!(connected[0].unit_b, 1);
    assert_eq!(connected[0].edge_count, 1);
    assert_eq!(connected[0].connected_indices, vec![0]);
    assert!(connected[0].has_edges(0));
    assert_eq!(
        connected[0].get_edges(0),
        &[InterUnitBondInfo {
            index_b: 1,
            props: InterUnitBondProps {
                order: 2,
                flag: BondFlags::COVALENT.union(BondFlags::AROMATIC),
                key: 7,
            },
        }]
    );
    assert_eq!(connected[1].unit_a, 0);
    assert_eq!(connected[1].unit_b, 2);
    assert_eq!(connected[1].get_edges(0)[0].index_b, 2);
    assert!(!graph.get_connected_units(1)[0].are_units_ordered());
}

#[test]
fn index_pair_distance_metadata_filters_unit_bonds_like_molstar() {
    let base = Molecule {
        atoms: vec![
            test_atom(1, "C1", "A", 1, vec3(0.0, 0.0, 0.0)),
            test_atom(2, "C2", "A", 1, vec3(10.0, 0.0, 0.0)),
        ],
        bonds: vec![Bond { a: 0, b: 1 }],
        selected_assembly: None,
        ..Molecule::default()
    };

    let mut rejected = base.clone();
    rejected.bond_metadata = vec![BondMetadata {
        source: BondSource::IndexPair,
        order: 1,
        flags: BondFlags::COVALENT,
        key: 1,
        distance: Some(1.0),
        operator_a: -1,
        operator_b: -1,
        struct_conn: None,
    }];
    assert_eq!(rejected.atomic_structure().intra_unit_bond_count, 0);

    let mut accepted = base;
    accepted.bond_metadata = vec![BondMetadata {
        source: BondSource::IndexPair,
        order: 1,
        flags: BondFlags::COVALENT,
        key: 2,
        distance: Some(10.2),
        operator_a: -1,
        operator_b: -1,
        struct_conn: None,
    }];
    assert_eq!(accepted.atomic_structure().intra_unit_bond_count, 1);
}

#[test]
fn index_pair_sidecar_max_distance_filters_bonds_when_pair_distance_is_absent() {
    let mut molecule = Molecule {
        atoms: vec![
            test_atom(1, "C1", "A", 1, vec3(0.0, 0.0, 0.0)),
            test_atom(2, "C2", "A", 1, vec3(4.0, 0.0, 0.0)),
            test_atom(3, "C3", "A", 1, vec3(10.0, 0.0, 0.0)),
        ],
        bonds: vec![Bond { a: 0, b: 1 }, Bond { a: 0, b: 2 }],
        bond_metadata: vec![
            BondMetadata {
                source: BondSource::IndexPair,
                order: 1,
                flags: BondFlags::COVALENT,
                key: 1,
                distance: None,
                operator_a: -1,
                operator_b: -1,
                struct_conn: None,
            },
            BondMetadata {
                source: BondSource::IndexPair,
                order: 1,
                flags: BondFlags::COVALENT,
                key: 2,
                distance: None,
                operator_a: -1,
                operator_b: -1,
                struct_conn: None,
            },
        ],
        ..Molecule::default()
    };
    molecule.index_pair_bonds = IndexPairBonds::from_bonds(
        &molecule.bonds,
        &molecule.bond_metadata,
        molecule.atoms.len(),
        5.0,
        true,
    );
    assert_eq!(molecule.atomic_structure().intra_unit_bond_count, 1);

    molecule.index_pair_bonds = IndexPairBonds::from_bonds(
        &molecule.bonds,
        &molecule.bond_metadata,
        molecule.atoms.len(),
        f32::INFINITY,
        true,
    );
    assert_eq!(molecule.atomic_structure().intra_unit_bond_count, 2);
}

#[test]
fn index_pair_sidecar_derives_operator_groups_like_molstar() {
    let bonds = vec![Bond { a: 0, b: 1 }, Bond { a: 0, b: 2 }];
    let metadata = vec![
        BondMetadata {
            source: BondSource::IndexPair,
            order: 1,
            flags: BondFlags::COVALENT,
            key: 1,
            distance: None,
            operator_a: 1,
            operator_b: 1,
            struct_conn: None,
        },
        BondMetadata {
            source: BondSource::IndexPair,
            order: 1,
            flags: BondFlags::COVALENT,
            key: 2,
            distance: None,
            operator_a: 1,
            operator_b: 2,
            struct_conn: None,
        },
        BondMetadata::computed(),
    ];
    let index_pairs = IndexPairBonds::from_bonds(&bonds, &metadata, 3, -1.0, true).unwrap();

    assert_eq!(index_pairs.bonds.edge_count, 2);
    assert!(index_pairs.contains_bond(0));
    assert!(index_pairs.contains_bond(1));
    assert!(index_pairs.has_operators);
    assert_eq!(index_pairs.by_same_operator.get(&1), Some(&vec![0, 2]));
}

#[test]
fn index_pair_graph_schema_matches_molstar_directed_edge_builder() {
    let metadata = vec![
        BondMetadata {
            source: BondSource::IndexPair,
            order: 2,
            flags: BondFlags::COVALENT.union(BondFlags::AROMATIC),
            key: 7,
            distance: Some(1.5),
            operator_a: 2,
            operator_b: 1,
            struct_conn: None,
        },
        BondMetadata {
            source: BondSource::IndexPair,
            order: 1,
            flags: BondFlags::COVALENT,
            key: 8,
            distance: None,
            operator_a: 1,
            operator_b: 1,
            struct_conn: None,
        },
    ];
    let index_pairs =
        IndexPairBonds::from_pairs(&[1, 0], &[0, 2], &[4, 9], &metadata, 3, -1.0, true).unwrap();

    assert_eq!(index_pairs.bonds.vertex_count, 3);
    assert_eq!(index_pairs.bonds.offset, vec![0, 2, 3, 4]);
    assert_eq!(index_pairs.bonds.a, vec![0, 0, 1, 2]);
    assert_eq!(index_pairs.bonds.b, vec![1, 2, 0, 0]);
    assert_eq!(index_pairs.bonds.edge_count, 2);
    assert!(index_pairs.contains_bond(4));
    assert!(index_pairs.contains_bond(9));
    assert_eq!(index_pairs.bonds.props.key, vec![7, 8, 7, 8]);
    assert_eq!(index_pairs.bonds.props.operator_a, vec![1, 1, 2, 1]);
    assert_eq!(index_pairs.bonds.props.operator_b, vec![2, 1, 1, 1]);
    assert_eq!(index_pairs.bonds.props.order, vec![2, 1, 2, 1]);
    assert_eq!(index_pairs.bonds.props.distance, vec![1.5, -1.0, 1.5, -1.0]);
    assert!(index_pairs.has_operators);
    assert_eq!(index_pairs.by_same_operator.get(&1), Some(&vec![1, 3]));
    assert_eq!(
        index_pairs.get_edge_index_for_operators(1, 0, 2, 1),
        Some(0)
    );
    assert_eq!(
        index_pairs.get_edge_index_for_operators(0, 1, 1, 2),
        Some(0)
    );
    assert_eq!(index_pairs.get_edge_index_for_operators(0, 1, 2, 1), None);
}

#[test]
fn intra_unit_bond_storage_uses_molstar_int_adjacency_graph_slots_and_metadata() {
    let mut molecule = Molecule {
        atoms: vec![
            test_atom(1, "C1", "A", 1, vec3(0.0, 0.0, 0.0)),
            test_atom(2, "C2", "A", 1, vec3(1.0, 0.0, 0.0)),
            test_atom(3, "C3", "A", 1, vec3(2.0, 0.0, 0.0)),
        ],
        bonds: vec![Bond { a: 0, b: 1 }, Bond { a: 0, b: 2 }],
        bond_metadata: vec![
            BondMetadata {
                source: BondSource::IndexPair,
                order: 2,
                flags: BondFlags::COVALENT.union(BondFlags::AROMATIC),
                key: 7,
                distance: None,
                operator_a: -1,
                operator_b: -1,
                struct_conn: None,
            },
            BondMetadata {
                source: BondSource::IndexPair,
                order: 1,
                flags: BondFlags::COVALENT,
                key: 8,
                distance: None,
                operator_a: -1,
                operator_b: -1,
                struct_conn: None,
            },
        ],
        ..Molecule::default()
    };
    molecule.index_pair_bonds = IndexPairBonds::from_bonds(
        &molecule.bonds,
        &molecule.bond_metadata,
        molecule.atoms.len(),
        -1.0,
        true,
    );
    let structure = molecule.atomic_structure();
    let bonds = &structure.units[0].props.intra_unit_bonds;

    assert_eq!(structure.intra_unit_bond_count, 2);
    assert_eq!(structure.units[0].props.intra_unit_bond_count, 2);
    assert_eq!(bonds.vertex_count, 3);
    assert_eq!(bonds.offset, vec![0, 2, 3, 4]);
    assert_eq!(bonds.a, vec![0, 0, 1, 2]);
    assert_eq!(bonds.b, vec![1, 2, 0, 0]);
    assert_eq!(bonds.edge_count, 2);
    assert_eq!(bonds.props.key, vec![7, 8, 7, 8]);
    assert_eq!(bonds.props.order, vec![2, 1, 2, 1]);
    assert_eq!(
        bonds
            .props
            .flags
            .iter()
            .map(|flags| flags.bits)
            .collect::<Vec<_>>(),
        vec![
            BondFlags::COVALENT.union(BondFlags::AROMATIC).bits,
            BondFlags::COVALENT.bits,
            BondFlags::COVALENT.union(BondFlags::AROMATIC).bits,
            BondFlags::COVALENT.bits,
        ]
    );
    assert!(!bonds.can_remap);
    assert!(bonds.cacheable);
}

#[test]
fn index_pair_same_operator_metadata_does_not_create_inter_unit_bonds() {
    let assembly = Assembly {
        id: "1".to_string(),
        details: String::new(),
        oligomeric_details: String::new(),
        oligomeric_count: None,
        asym_ids: vec!["A".to_string(), "B".to_string()],
        transforms: Vec::new(),
        generators: vec![AssemblyGenerator::from_transforms(
            "1",
            vec!["A".to_string(), "B".to_string()],
            0,
            vec![Transform::identity()],
            vec![vec!["1".to_string()]],
        )],
    };
    let structure = Molecule {
        atoms: vec![
            test_atom(1, "C1", "A", 1, vec3(0.0, 0.0, 0.0)),
            test_atom(2, "C2", "B", 1, vec3(1.0, 0.0, 0.0)),
        ],
        bonds: vec![Bond { a: 0, b: 1 }],
        bond_metadata: vec![BondMetadata {
            source: BondSource::IndexPair,
            order: 1,
            flags: BondFlags::COVALENT,
            key: 9,
            distance: None,
            operator_a: 1,
            operator_b: 1,
            struct_conn: None,
        }],
        selected_assembly: Some(assembly),
        ..Molecule::default()
    }
    .atomic_structure();

    assert_eq!(structure.units.len(), 2);
    assert_eq!(structure.intra_unit_bond_count, 0);
    assert!(structure.inter_unit_bonds.is_empty());
}

#[test]
fn cif_assembly_composes_operator_expression_groups() {
    let groups = expand_oper_expression("(X0)(1-2)");
    assert_eq!(
        groups,
        vec![
            vec!["X0".to_string(), "1".to_string()],
            vec!["X0".to_string(), "2".to_string()],
        ]
    );
    assert_eq!(
        expand_oper_expression("(1-2)(3-4)"),
        vec![
            vec!["1".to_string(), "3".to_string()],
            vec!["2".to_string(), "3".to_string()],
            vec!["1".to_string(), "4".to_string()],
            vec!["2".to_string(), "4".to_string()],
        ]
    );
    assert_eq!(
        expand_oper_expression("1,2"),
        vec![vec!["1".to_string()], vec!["2".to_string()]]
    );
    assert_eq!(
        expand_oper_expression("(1-3)"),
        vec![
            vec!["1".to_string()],
            vec!["2".to_string()],
            vec!["3".to_string()],
        ]
    );
    assert_eq!(
        expand_oper_expression("(1-2)(3-4)(5-6)"),
        vec![
            vec!["1".to_string(), "3".to_string(), "5".to_string()],
            vec!["2".to_string(), "3".to_string(), "5".to_string()],
            vec!["1".to_string(), "4".to_string(), "5".to_string()],
            vec!["2".to_string(), "4".to_string(), "5".to_string()],
            vec!["1".to_string(), "3".to_string(), "6".to_string()],
            vec!["2".to_string(), "3".to_string(), "6".to_string()],
            vec!["1".to_string(), "4".to_string(), "6".to_string()],
            vec!["2".to_string(), "4".to_string(), "6".to_string()],
        ]
    );

    let op_map = vec![
        (
            "X0".to_string(),
            Transform {
                m: [
                    [1.0, 0.0, 0.0, 10.0],
                    [0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0],
                ],
            },
        ),
        (
            "2".to_string(),
            Transform {
                m: [
                    [2.0, 0.0, 0.0, 0.0],
                    [0.0, 2.0, 0.0, 0.0],
                    [0.0, 0.0, 2.0, 0.0],
                ],
            },
        ),
    ];
    let transform = compose_operator_transforms(&groups[1], &op_map).unwrap();
    let out = transform.apply(Vec3 {
        x: 1.0,
        y: 1.0,
        z: 1.0,
    });

    assert_eq!(out.x, 12.0);
    assert_eq!(out.y, 2.0);
    assert_eq!(out.z, 2.0);

    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 1.000 1.000 1.000\n#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n1 (X0)(2) A\n#\nloop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\nX0 1 0 0 10 0 1 0 0 0 0 1 0\n2 2 0 0 0 0 2 0 0 0 0 2 0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: Some("1".to_string()),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = mol.atomic_structure();
    assert_eq!(structure.model.hierarchy.atoms.len(), 1);
    assert_eq!(structure.element_count, 1);
    assert_eq!(structure.units[0].operator.oper_list_ids, vec!["X0", "2"]);
    assert_eq!(structure.position(0, 0).unwrap(), vec3(12.0, 2.0, 2.0));
    let geometry = mol.expanded_for_geometry();
    assert_eq!(geometry.atoms.len(), 1);
    assert_eq!(geometry.atoms[0].position, vec3(12.0, 2.0, 2.0));
}

fn parse_operator_matrix_fixture(assembly: &str) -> AtomicStructure {
    parse_molecule_with_options(
        include_bytes!("../../tests/fixtures/cif/assembly-operator-matrix.cif"),
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: Some(assembly.to_string()),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap()
    .atomic_structure()
}

fn operator_fixture_signature(assembly: &str) -> Vec<(i32, Vec<String>, Vec3)> {
    parse_operator_matrix_fixture(assembly)
        .units
        .into_iter()
        .map(|unit| {
            (
                unit.operator.oper_id,
                unit.operator.oper_list_ids,
                unit.operator.transform.apply(vec3(1.0, 2.0, 3.0)),
            )
        })
        .collect()
}

#[test]
fn biological_assembly_single_operator_fixture_expands_one_operator() {
    assert_eq!(
        operator_fixture_signature("single"),
        vec![(1, vec!["2".to_string()], vec3(11.0, 2.0, 3.0))]
    );
}

#[test]
fn biological_assembly_operator_range_fixture_expands_in_molstar_order() {
    assert_eq!(
        operator_fixture_signature("range"),
        vec![
            (1, vec!["1".to_string()], vec3(1.0, 2.0, 3.0)),
            (2, vec!["2".to_string()], vec3(11.0, 2.0, 3.0)),
            (3, vec!["3".to_string()], vec3(21.0, 2.0, 3.0)),
        ]
    );
}

#[test]
fn biological_assembly_cartesian_product_fixture_preserves_operand_ids() {
    assert_eq!(
        operator_fixture_signature("cart"),
        vec![
            (
                1,
                vec!["X0".to_string(), "1".to_string()],
                vec3(101.0, 2.0, 3.0)
            ),
            (
                2,
                vec!["X0".to_string(), "2".to_string()],
                vec3(111.0, 2.0, 3.0)
            ),
        ]
    );
}

#[test]
fn assembly_fixture_reference_summary_matches_molstar_snapshot() {
    let actual = assembly_fixture_reference_summary_json();
    let expected =
        include_str!("../../tests/expected/assembly-fixture-reference-summary.json").trim_end();
    assert_eq!(actual, expected);
}

#[test]
fn cif_secondary_structure_uses_auth_fields_and_filters_non_helices() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA X 10 0.000 0.000 0.000\nATOM 2 C CA GLY X 20 1.000 0.000 0.000\n#\nloop_\n_struct_conf.conf_type_id\n_struct_conf.beg_auth_asym_id\n_struct_conf.beg_auth_seq_id\n_struct_conf.pdbx_beg_PDB_ins_code\n_struct_conf.end_auth_seq_id\n_struct_conf.pdbx_end_PDB_ins_code\nHELX_P X 10 A 20 B\nTURN_P X 30 . 40 ?\n#\nloop_\n_struct_sheet_range.sheet_id\n_struct_sheet_range.beg_auth_asym_id\n_struct_sheet_range.beg_auth_seq_id\n_struct_sheet_range.pdbx_beg_PDB_ins_code\n_struct_sheet_range.end_auth_seq_id\n_struct_sheet_range.pdbx_end_PDB_ins_code\nS1 X 50 C 60 D\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.helices.len(), 1);
    assert_eq!(mol.helices[0].chain, "X");
    assert_eq!(mol.helices[0].start, 10);
    assert_eq!(mol.helices[0].start_insertion_code, "A");
    assert_eq!(mol.helices[0].end, 20);
    assert_eq!(mol.helices[0].end_insertion_code, "B");
    assert_eq!(mol.sheets.len(), 1);
    assert_eq!(mol.sheets[0].start, 50);
    assert_eq!(mol.sheets[0].start_insertion_code, "C");
    assert_eq!(mol.sheets[0].end, 60);
    assert_eq!(mol.sheets[0].end_insertion_code, "D");
}

#[test]
fn cif_struct_conf_strn_maps_to_sheet_range() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A 1 0.000 0.000 0.000\nATOM 2 C CA VAL A 2 1.000 0.000 0.000\nATOM 3 C CA GLY A 3 2.000 0.000 0.000\nATOM 4 C CA SER A 4 3.000 0.000 0.000\nATOM 5 C CA THR A 5 4.000 0.000 0.000\nATOM 6 C CA LEU A 6 5.000 0.000 0.000\n#\nloop_\n_struct_conf.conf_type_id\n_struct_conf.beg_label_asym_id\n_struct_conf.beg_label_seq_id\n_struct_conf.pdbx_beg_PDB_ins_code\n_struct_conf.end_label_seq_id\n_struct_conf.pdbx_end_PDB_ins_code\nHELX_RH_AL_P A 1 . 2 ?\nSTRN A 3 A 4 B\nturn_p A 5 . 5 ?\nBEND A 6 . 6 ?\n#\nloop_\n_struct_sheet_range.sheet_id\n_struct_sheet_range.beg_label_asym_id\n_struct_sheet_range.beg_label_seq_id\n_struct_sheet_range.pdbx_beg_PDB_ins_code\n_struct_sheet_range.end_label_seq_id\n_struct_sheet_range.pdbx_end_PDB_ins_code\nS1 A 5 C 6 D\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.helices.len(), 1);
    assert_eq!(mol.helices[0].chain, "A");
    assert_eq!(mol.helices[0].start, 1);
    assert_eq!(mol.helices[0].end, 2);
    assert_eq!(mol.sheets.len(), 2);
    assert_eq!(mol.sheets[0].chain, "A");
    assert_eq!(mol.sheets[0].start, 3);
    assert_eq!(mol.sheets[0].start_insertion_code, "A");
    assert_eq!(mol.sheets[0].end, 4);
    assert_eq!(mol.sheets[0].end_insertion_code, "B");
    assert_eq!(mol.sheets[1].start, 5);
    assert_eq!(mol.sheets[1].start_insertion_code, "C");
    assert_eq!(mol.sheets[1].end, 6);
    assert_eq!(mol.sheets[1].end_insertion_code, "D");
}

#[test]
fn cif_secondary_structure_falls_back_to_auth_when_label_seq_is_null() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A X . 10 0.000 0.000 0.000\nATOM 2 C CA GLY A X . 20 1.000 0.000 0.000\n#\nloop_\n_struct_conf.conf_type_id\n_struct_conf.beg_label_asym_id\n_struct_conf.beg_label_seq_id\n_struct_conf.beg_auth_asym_id\n_struct_conf.beg_auth_seq_id\n_struct_conf.end_label_seq_id\n_struct_conf.end_auth_seq_id\n_struct_conf.pdbx_beg_PDB_ins_code\n_struct_conf.pdbx_end_PDB_ins_code\nHELX_P A . X 10 ? 20 A B\n#\nloop_\n_struct_sheet_range.sheet_id\n_struct_sheet_range.beg_label_asym_id\n_struct_sheet_range.beg_label_seq_id\n_struct_sheet_range.beg_auth_asym_id\n_struct_sheet_range.beg_auth_seq_id\n_struct_sheet_range.end_label_seq_id\n_struct_sheet_range.end_auth_seq_id\n_struct_sheet_range.pdbx_beg_PDB_ins_code\n_struct_sheet_range.pdbx_end_PDB_ins_code\nS1 A ? X 30 . 40 C D\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.helices.len(), 1);
    assert_eq!(mol.helices[0].chain, "A");
    assert_eq!(mol.helices[0].start, 10);
    assert_eq!(mol.helices[0].start_insertion_code, "A");
    assert_eq!(mol.helices[0].end, 20);
    assert_eq!(mol.helices[0].end_insertion_code, "B");
    assert_eq!(mol.sheets.len(), 1);
    assert_eq!(mol.sheets[0].chain, "A");
    assert_eq!(mol.sheets[0].start, 30);
    assert_eq!(mol.sheets[0].start_insertion_code, "C");
    assert_eq!(mol.sheets[0].end, 40);
    assert_eq!(mol.sheets[0].end_insertion_code, "D");
}

#[test]
fn pdb_secondary_structure_uses_fixed_columns_for_adjacent_fields() {
    fn set_field(line: &mut [u8], start: usize, value: &str) {
        for (offset, byte) in value.bytes().enumerate() {
            line[start + offset] = byte;
        }
    }

    fn pdb_helix_line(chain: &str, start: &str, end: &str) -> String {
        let mut line = vec![b' '; 80];
        set_field(&mut line, 0, "HELIX ");
        set_field(&mut line, 7, "  1");
        set_field(&mut line, 11, "H1 ");
        set_field(&mut line, 15, "MSEB");
        set_field(&mut line, 19, chain);
        set_field(&mut line, 21, start);
        set_field(&mut line, 25, "A");
        set_field(&mut line, 27, "GLY ");
        set_field(&mut line, 31, chain);
        set_field(&mut line, 33, end);
        set_field(&mut line, 37, "B");
        set_field(&mut line, 38, " 1");
        String::from_utf8(line).unwrap()
    }

    fn pdb_sheet_line(chain: &str, start: &str, end: &str) -> String {
        let mut line = vec![b' '; 80];
        set_field(&mut line, 0, "SHEET ");
        set_field(&mut line, 7, "  1");
        set_field(&mut line, 11, "S1 ");
        set_field(&mut line, 14, " 1");
        set_field(&mut line, 17, "SERX");
        set_field(&mut line, 21, chain);
        set_field(&mut line, 22, start);
        set_field(&mut line, 26, "A");
        set_field(&mut line, 28, "THR ");
        set_field(&mut line, 32, chain);
        set_field(&mut line, 33, end);
        set_field(&mut line, 37, "B");
        String::from_utf8(line).unwrap()
    }

    let pdb = format!(
        "{}\n{}\n{}\n{}\nATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 10.00           C\nEND\n",
        pdb_helix_line("B", "  -3", "  12"),
        pdb_sheet_line("C", "   7", "   9"),
        pdb_helix_line(" ", "   1", "   2"),
        pdb_sheet_line(" ", "   3", "   4")
    );
    let mol = parse_molecule_with_options(
        pdb.as_bytes(),
        &MeshOptions {
            format: InputFormat::Pdb,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.helices.len(), 2);
    assert_eq!(mol.helices[0].chain, "B");
    assert_eq!(mol.helices[0].start, -3);
    assert_eq!(mol.helices[0].start_insertion_code, "A");
    assert_eq!(mol.helices[0].end, 12);
    assert_eq!(mol.helices[0].end_insertion_code, "B");
    assert_eq!(mol.helices[1].chain, "");
    assert_eq!(mol.helices[1].start, 1);
    assert_eq!(mol.helices[1].start_insertion_code, "A");
    assert_eq!(mol.helices[1].end, 2);
    assert_eq!(mol.helices[1].end_insertion_code, "B");
    assert_eq!(mol.sheets.len(), 2);
    assert_eq!(mol.sheets[0].chain, "C");
    assert_eq!(mol.sheets[0].start, 7);
    assert_eq!(mol.sheets[0].start_insertion_code, "A");
    assert_eq!(mol.sheets[0].end, 9);
    assert_eq!(mol.sheets[0].end_insertion_code, "B");
    assert_eq!(mol.sheets[1].chain, "");
    assert_eq!(mol.sheets[1].start, 3);
    assert_eq!(mol.sheets[1].start_insertion_code, "A");
    assert_eq!(mol.sheets[1].end, 4);
    assert_eq!(mol.sheets[1].end_insertion_code, "B");
}

#[test]
fn binary_cif_column_masks_roundtrip_to_cif_null_tokens() {
    let values = ColumnData::Str(vec!["A".to_string(), "B".to_string(), "C".to_string()])
        .with_mask(&[0, 1, 2]);

    assert_eq!(values.string_at(0), "A");
    assert_eq!(values.string_at(1), ".");
    assert_eq!(values.string_at(2), "?");
}

#[test]
fn binary_cif_typed_columns_preserve_numeric_access_under_masks() {
    let values = ColumnData::Float(vec![1.25, 2.5, 3.75]).with_mask(&[0, 1, 0]);

    assert_eq!(values.f32_at(0), Some(1.25));
    assert_eq!(values.f32_at(1), None);
    assert_eq!(values.f32_at(2), Some(3.75));
    assert_eq!(values.string_at(1), ".");

    let ids = ColumnData::Int(vec![1, 2, 3]);
    let ranges = ColumnData::Int(vec![10, 20, 30]).with_mask(&[0, 2, 0]);
    assert_eq!(ids.usize_at(2), Some(3));
    assert_eq!(ranges.i32_at(0), Some(10));
    assert_eq!(ranges.i32_at(1), None);
    assert_eq!(ranges.i32_at(2), Some(30));
}

#[test]
fn binary_cif_tables_parse_directly_with_assembly_alt_loc_and_bonds() {
    let bcif = include_bytes!("../../tests/fixtures/bcif/assembly-altloc-helix.bcif");
    let mol = parse_molecule_with_options(
        bcif,
        &MeshOptions {
            format: InputFormat::BinaryCif,
            assembly: Some("1".to_string()),
            alt_loc: "A".to_string(),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.atoms.len(), 13);
    assert_eq!(mol.expanded_for_geometry().atoms.len(), 26);
    assert!(mol.selected_assembly.is_some());
    assert_eq!(mol.helices.len(), 1);
    assert_eq!(mol.helices[0].chain, "A");
    assert_eq!(mol.helices[0].start, 1);
    assert_eq!(mol.helices[0].start_insertion_code, "A");
    assert_eq!(mol.helices[0].end, 3);
    assert_eq!(mol.helices[0].end_insertion_code, "B");
    assert_eq!(mol.sheets.len(), 1);
    assert_eq!(mol.sheets[0].chain, "A");
    assert_eq!(mol.sheets[0].start, 2);
    assert_eq!(mol.sheets[0].start_insertion_code, "");
    assert_eq!(mol.sheets[0].end, 3);
    assert_eq!(mol.sheets[0].end_insertion_code, "B");
    assert_eq!(mol.atoms[0].auth_name, "NX");
    assert_eq!(mol.atoms[0].auth_residue, "ALAX");
    assert_eq!(mol.atoms[0].auth_chain, "X");
    assert_eq!(mol.atoms[0].entity_id, "1");
    assert_eq!(mol.atoms[0].auth_residue_seq, "101");
    assert_eq!(mol.atoms[0].insertion_code, "A");
    assert!(mol
        .atoms
        .iter()
        .any(|atom| atom.residue_seq == "3" && atom.insertion_code == "B"));
    assert!(!mol.bonds.is_empty());
    let structure = mol.atomic_structure();
    assert_eq!(structure.model.hierarchy.chains[0].auth_id, "X");
    assert_eq!(structure.model.hierarchy.residues[0].auth_seq_id, "101");
    assert_eq!(structure.model.hierarchy.residues[0].insertion_code, "A");
    assert_eq!(structure.properties.auth_atom_id[0], "NX");
    assert_eq!(structure.properties.auth_comp_id[0], "ALAX");
    assert_eq!(structure.properties.auth_asym_id[0], "X");
    assert_eq!(structure.properties.auth_seq_id[0], "101");
    assert_eq!(structure.properties.label_entity_id[0], "1");
    assert_eq!(structure.properties.pdbx_pdb_ins_code[0], "A");
    assert_eq!(mol.entities.len(), 1);
    assert_eq!(mol.entries.len(), 1);
    assert_eq!(mol.experiments.len(), 1);
    assert_eq!(mol.entity_polymers.len(), 1);
    assert_eq!(mol.entity_poly_seq.len(), 3);
    assert_eq!(mol.pdbx_entity_branch.len(), 1);
    assert_eq!(mol.pdbx_entity_branch[0].entity_id, "2");
    assert_eq!(mol.pdbx_entity_branch_links.len(), 1);
    assert_eq!(mol.pdbx_entity_branch_links[0].link_id, 1);
    assert_eq!(mol.pdbx_entity_branch_links[0].value_order, "sing");
    assert_eq!(mol.pdbx_branch_scheme.len(), 1);
    assert_eq!(mol.pdbx_branch_scheme[0].num, 1);
    assert_eq!(mol.pdbx_nonpoly_scheme.len(), 1);
    assert_eq!(mol.pdbx_nonpoly_scheme[0].pdb_ins_code, "");
    assert_eq!(mol.pdbx_poly_seq_scheme.len(), 1);
    assert_eq!(mol.pdbx_poly_seq_scheme[0].seq_id, 1);
    assert_eq!(mol.struct_asym.len(), 1);
    assert_eq!(mol.chemical_components.len(), 3);
    assert_eq!(mol.chemical_components[0].name, "ALANINE");
    assert_eq!(mol.chemical_components[0].formula_weight, Some(89.09));
    assert_eq!(mol.chemical_components[0].mon_nstd_flag, "n");
    assert_eq!(mol.chemical_components[0].pdbx_release_status, "REL");
    assert_eq!(mol.chemical_component_atoms.len(), 5);
    assert_eq!(mol.chemical_component_bonds.len(), 4);
    assert_eq!(mol.chemical_component_angles.len(), 2);
    assert_eq!(mol.chemical_component_bonds[2].order, 2);
    assert_eq!(mol.chemical_component_bonds[3].ordinal, Some(104));
    assert!(mol.chemical_component_bonds[3]
        .flags
        .contains(BondFlags::AROMATIC));
    assert_eq!(mol.chemical_component_atoms[1].charge, 1);
    assert_eq!(
        mol.chemical_component_atoms[1].model_cartn,
        Some(vec3(1.45, 0.05, 0.1))
    );
    assert_eq!(mol.chemical_component_angles[0].value_angle, Some(111.0));
    assert_eq!(mol.atoms[0].entity_id, "1");
    assert_eq!(mol.atoms[1].formal_charge, 1);
    assert_eq!(mol.atoms[1].b_iso, 11.0);
    assert_eq!(mol.atom_site_anisotrop.len(), 1);
    assert_eq!(mol.atom_site_anisotrop[0].atom_id, 2);
    let structure = mol.atomic_structure();
    assert_eq!(
        structure.model.conformation.atom_ids,
        vec![1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14]
    );
    assert_eq!(
        structure.model.conformation.element_to_anisotrop,
        vec![-1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    );
    assert_eq!(
        structure.model.conformation.anisotropic_displacement[1],
        Some([[0.10, 0.01, 0.02], [0.01, 0.11, 0.03], [0.02, 0.03, 0.12]])
    );
}

#[test]
fn binary_cif_ihm_metadata_parses_directly() {
    let bcif = include_bytes!("../../tests/fixtures/bcif/ihm-only.bcif");
    let mol = parse_molecule_with_options(
        bcif,
        &MeshOptions {
            format: InputFormat::BinaryCif,
            center: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(mol.atoms.len(), 0);
    assert_eq!(mol.ihm_model_list.len(), 1);
    assert_eq!(mol.ihm_model_list[0].model_id, 1);
    assert_eq!(mol.ihm_model_list[0].model_name, "model one");
    assert_eq!(mol.ihm_model_groups.len(), 1);
    assert_eq!(mol.ihm_model_groups[0].id, 1);
    assert_eq!(mol.ihm_model_group_links.len(), 1);
    assert_eq!(mol.ihm_model_group_links[0].model_id, 1);
    assert_eq!(mol.ihm_cross_link_restraints.len(), 1);
    assert_eq!(mol.ihm_cross_link_restraints[0].seq_id_2, 10);
    assert_eq!(
        mol.ihm_cross_link_restraints[0].distance_threshold,
        Some(25.0)
    );

    let info = String::from_utf8(
        molecule_info(
            bcif,
            br#"{"format":"bcif","center":false,"assembly":"asymmetric-unit"}"#,
        )
        .unwrap(),
    )
    .unwrap();
    assert!(info.contains(r#""ihm_model_count":1"#));
    assert!(info.contains(r#""ihm_model_group_count":1"#));
    assert!(info.contains(r#""ihm_model_group_link_count":1"#));
    assert!(info.contains(r#""ihm_cross_link_restraint_count":1"#));
}

#[test]
fn info_exposes_semantic_render_object_summary() {
    let cif = include_bytes!("../../tests/fixtures/cif/assembly-altloc-helix.cif");
    let info = String::from_utf8(
            molecule_info(
                cif,
                br#"{"format":"cif","representation":"cartoon","assembly":"asymmetric-unit","alt-loc":"A"}"#,
            )
            .unwrap(),
        )
        .unwrap();
    assert!(info.contains("\"render_objects\""));
    assert!(info.contains("\"structure\":{\"model_count\":1"));
    assert!(info.contains("\"conformation\":{\"atom_id_count\":"));
    assert!(info.contains("\"secondary_type\":\"helix\""));
    assert!(info.contains("\"geometry_type\":\"tube\""));
    assert!(info.contains("\"representation\":\"cartoon\""));
    assert!(info.contains(r#""representation":{"name":"cartoon""#));
    assert!(info.contains(r#""selected_visuals":["polymer-trace"]"#));
    assert!(info.contains(r#""realized_visuals":["polymer-trace"]"#));
}

#[test]
fn atomic_structure_groups_atoms_into_hierarchy_units_and_ranges() {
    let mut atoms = vec![
        test_atom(1, "CA", "A", 1, vec3(0.0, 0.0, 0.0)),
        test_atom(2, "C", "A", 1, vec3(0.3, 0.0, 0.0)),
        test_atom(3, "CA", "A", 3, vec3(2.0, 0.0, 0.0)),
        test_atom(4, "O", "B", 1, vec3(0.0, 2.0, 0.0)),
    ];
    atoms[1].alt_id = "A".to_string();
    atoms[3].het = true;
    atoms[3].residue = "HOH".to_string();
    let molecule = Molecule {
        atoms,
        ..Molecule::default()
    };

    let structure = molecule.atomic_structure();
    assert_eq!(structure.model.hierarchy.atoms.len(), 4);
    assert_eq!(structure.model.hierarchy.residues.len(), 3);
    assert_eq!(structure.model.hierarchy.chains.len(), 2);
    assert_eq!(structure.model.conformation.positions.len(), 4);
    assert_eq!(structure.model.conformation.occupancies.len(), 4);
    assert_eq!(structure.units.len(), 2);
    assert_eq!(structure.units[0].kind, AtomicUnitKind::Atomic);
    assert_eq!(structure.units[0].elements, vec![0, 1, 2]);
    assert_eq!(structure.units[0].atom_indices, vec![0, 1, 2]);
    assert_eq!(
        structure.units[0].residue_index_by_element,
        structure.model.hierarchy.residue_atom_segments.index
    );
    assert_eq!(
        structure.units[0].chain_index_by_element,
        structure.model.hierarchy.chain_atom_segments.index
    );
    assert_eq!(structure.units[1].residue_index_by_element[3], 2);
    assert_eq!(structure.units[1].chain_index_by_element[3], 1);
    assert_eq!(structure.units[0].operator.name, "1_555");
    assert_eq!(structure.units[0].operator.instance_id, "1_555");
    assert_eq!(structure.units[0].operator.oper_id, -1);
    assert!(structure.units[0].operator.oper_list_ids.is_empty());
    assert!(structure.units[0].operator.is_identity);
    assert_eq!(
        structure.model.hierarchy.atom_source_index,
        vec![0, 1, 2, 3]
    );
    assert_eq!(
        structure.model.hierarchy.residue_source_index,
        vec![0, 1, 2]
    );
    assert_eq!(structure.model.hierarchy.residue_atom_segments.count, 3);
    assert_eq!(
        structure.model.hierarchy.residue_atom_segments.offsets,
        vec![0, 2, 3, 4]
    );
    assert_eq!(
        structure.model.hierarchy.residue_atom_segments.index,
        vec![0, 0, 1, 2]
    );
    assert_eq!(structure.model.hierarchy.chain_atom_segments.count, 2);
    assert_eq!(
        structure.model.hierarchy.chain_atom_segments.offsets,
        vec![0, 3, 4]
    );
    assert_eq!(
        structure.model.hierarchy.chain_atom_segments.index,
        vec![0, 0, 0, 1]
    );
    assert_eq!(structure.properties.residue_index, vec![0, 0, 1, 2]);
    assert_eq!(structure.properties.chain_index, vec![0, 0, 0, 1]);
    assert_eq!(structure.properties.label_atom_id[0], "CA");
    assert_eq!(structure.properties.label_comp_id[3], "HOH");
    assert_eq!(structure.properties.label_asym_id[3], "B");
    assert_eq!(structure.properties.len(), 4);
    assert!(!structure.properties.is_empty());
    assert_eq!(structure.properties.atom_key(1), Some(1));
    assert_eq!(structure.properties.atom_id(1), Some(2));
    assert_eq!(structure.properties.atom_source_index(1), Some(1));
    assert_eq!(structure.properties.atom_type_symbol(1), Some("C"));
    assert_eq!(structure.properties.atom_label_atom_id(1), Some("C"));
    assert_eq!(structure.properties.atom_auth_atom_id(1), Some("C"));
    assert_eq!(structure.properties.atom_label_alt_id(1), Some("A"));
    assert_eq!(structure.properties.atom_label_comp_id(3), Some("HOH"));
    assert_eq!(structure.properties.atom_auth_comp_id(3), Some("ALA"));
    assert_eq!(structure.properties.atom_formal_charge(1), Some(0));
    assert_eq!(structure.properties.residue_key(3), Some(2));
    assert_eq!(structure.properties.residue_group_pdb(3), Some("ATOM"));
    assert_eq!(structure.properties.residue_label_comp_id(3), Some("HOH"));
    assert_eq!(structure.properties.residue_auth_comp_id(3), Some("ALA"));
    assert_eq!(structure.properties.residue_label_seq_id(3), Some("1"));
    assert_eq!(structure.properties.residue_auth_seq_id(3), Some("1"));
    assert_eq!(structure.properties.residue_pdb_ins_code(3), Some(""));
    assert_eq!(structure.properties.chain_key(3), Some(1));
    assert_eq!(structure.properties.chain_label_asym_id(3), Some("B"));
    assert_eq!(structure.properties.chain_auth_asym_id(3), Some("B"));
    assert_eq!(structure.properties.chain_label_entity_id(3), Some(""));
    assert_eq!(structure.properties.atom_id(4), None);
    assert_eq!(structure.properties.chain_label_asym_id(4), None);
    assert_eq!(
        structure.model.hierarchy.derived.residue.molecule_type,
        vec![
            MoleculeType::Protein,
            MoleculeType::Protein,
            MoleculeType::Water
        ]
    );
    assert_eq!(
        structure
            .model
            .hierarchy
            .derived
            .residue
            .trace_element_index,
        vec![Some(0), Some(2), None]
    );
    assert_eq!(
        structure
            .model
            .hierarchy
            .derived
            .residue
            .direction_from_element_index[0],
        Some(1)
    );
    assert_eq!(structure.units[0].props.residue_count, 2);
    assert_eq!(structure.units[0].props.protein_elements, vec![0, 2]);
    assert_eq!(structure.units[0].props.polymer_elements, vec![0]);
    assert_eq!(structure.units[0].props.gap_elements, vec![0, 2]);
    assert!(structure.units[0].props.nucleotide_elements.is_empty());
    assert_eq!(structure.units[1].props.water_elements, vec![3]);
    assert_eq!(structure.ranges.polymer_ranges, vec![0, 1]);
    assert_eq!(structure.ranges.gap_ranges, vec![0, 2]);
    assert_eq!(structure.element_count, 4);
    assert_eq!(structure.polymer_gap_count, 1);
    assert_vec3_close(
        structure.units[0].props.boundary.box_min,
        vec3(-1.7, -1.7, -1.7),
        0.000_001,
    );
    assert_vec3_close(
        structure.units[0].props.boundary.box_max,
        vec3(3.7, 1.7, 1.7),
        0.000_001,
    );
    assert_eq!(structure.units[0].props.lookup3d.len(), 3);
    assert!(structure.units[0]
        .props
        .lookup3d
        .check(vec3(0.1, 0.0, 0.0), 0.2));
    assert_eq!(
        structure.units[0]
            .props
            .lookup3d
            .nearest(vec3(2.1, 0.0, 0.0), 1)[0]
            .index,
        2
    );
    assert_vec3_close(
        structure.boundary.box_min,
        vec3(-1.7, -1.7, -1.7),
        0.000_001,
    );
    assert_vec3_close(structure.boundary.box_max, vec3(3.7, 3.52, 1.7), 0.000_001);
    assert!(structure.lookup3d.check(vec3(0.0, 2.0, 0.0), 0.01));
    assert_eq!(
        structure.lookup3d.nearest(vec3(0.0, 2.1, 0.0), 1)[0].unit_id,
        1
    );
    assert_eq!(structure.symmetry_groups.len(), 2);
    assert_eq!(structure.symmetry_groups[0].model_id, 0);
    assert_eq!(
        structure.symmetry_groups[0].operator_instance_ids,
        vec!["1_555"]
    );
    assert_eq!(structure.symmetry_groups[0].unit_index_map, vec![(0, 0)]);
    assert_ne!(structure.symmetry_groups[0].hash_code, 0);
    assert_ne!(structure.symmetry_groups[0].transform_hash, 0);
    assert_eq!(structure.position(0, 2).unwrap(), vec3(2.0, 0.0, 0.0));
    assert_eq!(structure.alt_loc_count(), 1);
}

#[test]
fn boundary_union_contains_off_axis_input_spheres() {
    let a = Boundary {
        box_min: vec3(0.0, -2.0, 0.0),
        box_max: vec3(0.0, 2.0, 0.0),
        sphere: BoundingSphere {
            center: vec3(0.0, 0.0, 0.0),
            radius: 2.0,
            extrema: Vec::new(),
            center64: None,
            radius64: None,
            extrema64: Vec::new(),
        },
    };
    let b = Boundary {
        box_min: vec3(4.0, 0.0, -2.0),
        box_max: vec3(4.0, 0.0, 2.0),
        sphere: BoundingSphere {
            center: vec3(4.0, 0.0, 0.0),
            radius: 2.0,
            extrema: Vec::new(),
            center64: None,
            radius64: None,
            extrema64: Vec::new(),
        },
    };

    let union = a.union(b);
    for point in [
        vec3(0.0, -2.0, 0.0),
        vec3(0.0, 2.0, 0.0),
        vec3(4.0, 0.0, -2.0),
        vec3(4.0, 0.0, 2.0),
    ] {
        assert!(
            union.sphere.center.distance(point) <= union.sphere.radius + 0.000_1,
            "point {point:?} outside {union:?}"
        );
    }
}

#[test]
fn boundary_union_preserves_molstar_double_precision_extrema() {
    let sphere = BoundingSphere {
        center: vec3(0.5, 0.0, 0.0),
        radius: 0.5,
        extrema: vec![vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0)],
        center64: Some([0.5, 0.0, 0.0]),
        radius64: Some(0.5),
        extrema64: vec![[0.0, 0.0004, 0.0], [1.0, 0.0004, 0.0]],
    };

    let boundary = Boundary::from_bounding_spheres(&[sphere]);

    assert!((boundary.sphere.center.y - 0.0004).abs() < 1e-7);
    assert!(boundary
        .sphere
        .extrema64
        .iter()
        .any(|point| (point[1] - 0.0004).abs() < 1e-12));
}

#[test]
fn boundary_from_positions_matches_molstar_epos_centroid() {
    let boundary = Boundary::from_positions(&[
        vec3(0.0, 0.0, 0.0),
        vec3(3.0, 1.0, 0.0),
        vec3(-1.0, 2.0, 4.0),
        vec3(2.0, -2.0, 1.0),
    ]);

    assert_vec3_close(
        boundary.sphere.center,
        vec3(0.826_530_64, 0.459_183_66, 1.622_449),
        0.000_001,
    );
    assert!((boundary.sphere.radius - 3.370_916_4).abs() < 0.000_001);
    assert_eq!(boundary.box_min, vec3(-1.0, -2.0, 0.0));
    assert_eq!(boundary.box_max, vec3(3.0, 2.0, 4.0));
    assert_eq!(boundary.sphere.extrema.len(), 4);
}

#[test]
fn unit_lookup3d_uses_molstar_grid_bucket_layout_and_query_order() {
    let positions = vec![
        vec3(0.0, 0.0, 0.0),
        vec3(4.0, 0.0, 0.0),
        vec3(10.0, 0.0, 0.0),
        vec3(4.2, 0.0, 0.0),
        vec3(4.4, 0.0, 0.0),
    ];
    let lookup = UnitLookup3D::new(positions.clone(), Boundary::from_positions(&positions));

    assert_eq!(lookup.grid_size(), [5, 1, 1]);
    assert!((lookup.grid_delta().x - 2.2).abs() < 0.000_001);
    assert_eq!(lookup.bucket_offsets(), &[0, 1, 4]);
    assert_eq!(lookup.bucket_counts(), &[1, 3, 1]);
    assert_eq!(lookup.bucket_array(), &[0, 1, 3, 4, 2]);

    let hits = lookup.find(vec3(4.25, 0.0, 0.0), 0.3);
    assert_eq!(
        hits.iter().map(|hit| hit.index).collect::<Vec<_>>(),
        vec![1, 3, 4]
    );
    assert!(lookup.check(vec3(4.25, 0.0, 0.0), 0.051));
    assert!(!lookup.check(vec3(4.25, 0.0, 0.0), 0.04));
    assert_eq!(lookup.nearest(vec3(4.25, 0.0, 0.0), 1)[0].index, 3);
}

#[test]
fn structure_lookup3d_uses_molstar_unit_grid_layout_and_radius_filter() {
    let cif = b"data_demo\nloop_\n_ihm_sphere_obj_site.id\n_ihm_sphere_obj_site.entity_id\n_ihm_sphere_obj_site.asym_id\n_ihm_sphere_obj_site.seq_id_begin\n_ihm_sphere_obj_site.seq_id_end\n_ihm_sphere_obj_site.Cartn_x\n_ihm_sphere_obj_site.Cartn_y\n_ihm_sphere_obj_site.Cartn_z\n_ihm_sphere_obj_site.object_radius\n1 1 A 1 1 0.0 0.0 0.0 0.0\n2 1 B 1 1 4.0 0.0 0.0 0.0\n3 1 C 1 1 10.0 0.0 0.0 0.0\n4 1 D 1 1 4.2 0.0 0.0 1.0\n5 1 E 1 1 4.4 0.0 0.0 0.0\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = mol.atomic_structure();

    assert_eq!(structure.units.len(), 5);
    assert_eq!(structure.lookup3d.unit_grid_size(), [3, 1, 1]);
    assert_eq!(structure.lookup3d.unit_bucket_offsets(), &[0, 1, 4]);
    assert_eq!(structure.lookup3d.unit_bucket_counts(), &[1, 3, 1]);
    assert_eq!(structure.lookup3d.unit_bucket_array(), &[0, 1, 3, 4, 2]);
    assert_eq!(
        structure
            .lookup3d
            .close_unit_indices(vec3(5.15, 0.0, 0.0), 0.0),
        vec![3]
    );
    assert!(!structure.lookup3d.check(vec3(5.15, 0.0, 0.0), 0.0));

    let hits = structure.lookup3d.find(vec3(4.25, 0.0, 0.0), 0.3);
    assert_eq!(
        hits.iter().map(|hit| hit.unit_id).collect::<Vec<_>>(),
        vec![1, 3, 4]
    );
    assert_eq!(
        hits.iter().map(|hit| hit.index).collect::<Vec<_>>(),
        vec![0, 0, 0]
    );
}

#[test]
fn atomic_hierarchy_tables_match_molstar_column_shapes() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.auth_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.auth_comp_id\n_atom_site.pdbx_formal_charge\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA CAA A GLY GLY 1 A X 1 10 110 . 0.000 0.000 0.000\nATOM 2 C CB CBA B ALA ALB -1 A X 1 10 110 . 1.000 0.000 0.000\nHETATM 3 O O OA . HOH HOH 0 B Y 2 5 205 Z 2.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let structure = molecule.atomic_structure();
    let tables = structure.model.hierarchy.tables();

    assert_eq!(
        crate::model::AtomicAtomsTable::COLUMN_NAMES,
        &[
            "type_symbol",
            "label_atom_id",
            "auth_atom_id",
            "label_alt_id",
            "label_comp_id",
            "auth_comp_id",
            "pdbx_formal_charge"
        ]
    );
    assert_eq!(
        crate::model::AtomicResiduesTable::COLUMN_NAMES,
        &[
            "group_PDB",
            "label_seq_id",
            "auth_seq_id",
            "pdbx_PDB_ins_code"
        ]
    );
    assert_eq!(
        crate::model::AtomicChainsTable::COLUMN_NAMES,
        &["label_asym_id", "auth_asym_id", "label_entity_id"]
    );

    assert_eq!(tables.atoms.row_count(), 3);
    assert_eq!(tables.atoms.type_symbol, vec!["C", "C", "O"]);
    assert_eq!(tables.atoms.label_atom_id, vec!["CA", "CB", "O"]);
    assert_eq!(tables.atoms.auth_atom_id, vec!["CAA", "CBA", "OA"]);
    assert_eq!(tables.atoms.label_alt_id, vec!["A", "B", ""]);
    assert_eq!(tables.atoms.label_comp_id, vec!["GLY", "ALA", "HOH"]);
    assert_eq!(tables.atoms.auth_comp_id, vec!["GLY", "ALB", "HOH"]);
    assert_eq!(tables.atoms.pdbx_formal_charge, vec![1, -1, 0]);

    assert_eq!(tables.residues.row_count(), 2);
    assert_eq!(tables.residues.group_pdb, vec!["ATOM", "HETATM"]);
    assert_eq!(tables.residues.label_seq_id, vec!["10", "5"]);
    assert_eq!(tables.residues.auth_seq_id, vec!["110", "205"]);
    assert_eq!(tables.residues.pdbx_pdb_ins_code, vec!["", "Z"]);
    assert_eq!(structure.model.hierarchy.residues[0].comp_id, "GLY");
    assert_eq!(structure.model.hierarchy.residues[0].auth_comp_id, "GLY");
    assert_eq!(structure.model.hierarchy.atoms[1].label_comp_id, "ALA");
    assert_eq!(structure.model.hierarchy.atoms[1].auth_comp_id, "ALB");

    assert_eq!(tables.chains.row_count(), 2);
    assert_eq!(tables.chains.label_asym_id, vec!["A", "B"]);
    assert_eq!(tables.chains.auth_asym_id, vec!["X", "Y"]);
    assert_eq!(tables.chains.label_entity_id, vec!["1", "2"]);
}

#[test]
fn atomic_ranges_use_backbone_connectivity_and_cyclic_polymer_map() {
    let mut atoms = Vec::new();
    for (seq, n_x, ca_x, c_x, o_x) in [
        (1, 0.0, 0.5, 1.0, 1.1),
        (2, 8.0, 8.5, 9.0, 9.1),
        (3, 9.2, 9.5, 0.4, 0.5),
    ] {
        atoms.push(test_atom(
            atoms.len() + 1,
            "N",
            "A",
            seq,
            vec3(n_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "CA",
            "A",
            seq,
            vec3(ca_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "C",
            "A",
            seq,
            vec3(c_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "O",
            "A",
            seq,
            vec3(o_x, 0.0, 0.0),
        ));
    }
    for atom in &mut atoms {
        atom.entity_id = "1".to_string();
    }
    let entities = vec![Entity {
        id: "1".to_string(),
        type_name: "polymer".to_string(),
        description: String::new(),
    }];
    let entity_poly_seq = (1..=3)
        .map(|num| EntityPolySeq {
            entity_id: "1".to_string(),
            num,
            mon_id: "ALA".to_string(),
            hetero: "n".to_string(),
        })
        .collect::<Vec<_>>();
    let structure = Molecule {
        atoms,
        entity_index: EntityIndexMap::from_entities(&entities, &[], &[]),
        entities,
        entity_poly_seq,
        ..Molecule::default()
    }
    .atomic_structure();

    assert_eq!(structure.ranges.polymer_ranges, vec![0, 3, 4, 11]);
    assert_eq!(structure.ranges.gap_ranges, vec![0, 7]);
    assert_eq!(structure.ranges.cyclic_polymer_map.get(&0), Some(&2));
    assert_eq!(structure.ranges.cyclic_polymer_map.get(&2), Some(&0));
}

#[test]
fn atomic_ranges_match_molstar_polymer_range_generation() {
    let mut atoms = Vec::new();
    for (seq, n_x, ca_x, c_x) in [
        (1, 0.0, 0.5, 1.0),
        (2, 1.2, 1.7, 2.2),
        (3, 25.0, 25.5, 26.0),
        (4, 26.2, 26.7, 27.2),
    ] {
        atoms.push(test_atom(
            atoms.len() + 1,
            "N",
            "A",
            seq,
            vec3(n_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "CA",
            "A",
            seq,
            vec3(ca_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "C",
            "A",
            seq,
            vec3(c_x, 0.0, 0.0),
        ));
    }
    let structure = Molecule {
        atoms,
        ..Molecule::default()
    }
    .atomic_structure();

    assert_eq!(structure.ranges.polymer_ranges, vec![0, 5, 6, 11]);
    assert_eq!(structure.ranges.gap_ranges, vec![3, 8]);
    assert_eq!(structure.units[0].props.polymer_elements, vec![1, 4, 7, 10]);
    assert_eq!(structure.units[0].props.gap_elements, vec![4, 7]);
}

#[test]
fn atomic_gap_ranges_match_molstar_prev_to_current_residue_span() {
    let mut atoms = Vec::new();
    for (seq, n_x, ca_x, c_x) in [
        (1, 0.0, 0.5, 1.0),
        (3, 3.0, 3.5, 4.0),
        (4, 24.0, 24.5, 25.0),
        (5, 25.2, 25.7, 26.2),
    ] {
        atoms.push(test_atom(
            atoms.len() + 1,
            "N",
            "A",
            seq,
            vec3(n_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "CA",
            "A",
            seq,
            vec3(ca_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "C",
            "A",
            seq,
            vec3(c_x, 0.0, 0.0),
        ));
    }
    let structure = Molecule {
        atoms,
        ..Molecule::default()
    }
    .atomic_structure();

    assert_eq!(structure.ranges.polymer_ranges, vec![0, 2, 3, 5, 6, 11]);
    assert_eq!(structure.ranges.gap_ranges, vec![0, 5, 3, 8]);
    assert_eq!(structure.units[0].props.gap_elements, vec![1, 4, 4, 7]);
}

#[test]
fn atomic_ranges_cyclic_polymer_uses_entity_sequence_max_seq_id() {
    let mut atoms = Vec::new();
    for (seq, n_x, ca_x, c_x, o_x) in [
        (1, 0.0, 0.5, 1.0, 1.1),
        (2, 2.0, 2.5, 3.0, 3.1),
        (3, 3.2, 3.5, 0.4, 0.5),
    ] {
        atoms.push(test_atom(
            atoms.len() + 1,
            "N",
            "A",
            seq,
            vec3(n_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "CA",
            "A",
            seq,
            vec3(ca_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "C",
            "A",
            seq,
            vec3(c_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "O",
            "A",
            seq,
            vec3(o_x, 0.0, 0.0),
        ));
    }
    for atom in &mut atoms {
        atom.entity_id = "1".to_string();
    }
    let entities = vec![Entity {
        id: "1".to_string(),
        type_name: "polymer".to_string(),
        description: String::new(),
    }];
    let entity_poly_seq = (1..=4)
        .map(|num| EntityPolySeq {
            entity_id: "1".to_string(),
            num,
            mon_id: "ALA".to_string(),
            hetero: "n".to_string(),
        })
        .collect::<Vec<_>>();

    let structure = Molecule {
        atoms,
        entity_index: EntityIndexMap::from_entities(&entities, &[], &[]),
        entities,
        entity_poly_seq,
        ..Molecule::default()
    }
    .atomic_structure();

    assert!(structure.ranges.cyclic_polymer_map.is_empty());
}

#[test]
fn atomic_ranges_cyclic_polymer_uses_range_backed_coarse_sequence_like_molstar() {
    let mut atoms = Vec::new();
    for (seq, n_x, ca_x, c_x, o_x) in [
        (1, 0.0, 0.5, 1.0, 1.1),
        (2, 2.0, 2.5, 3.0, 3.1),
        (3, 3.2, 3.5, 0.4, 0.5),
    ] {
        atoms.push(test_atom(
            atoms.len() + 1,
            "N",
            "A",
            seq,
            vec3(n_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "CA",
            "A",
            seq,
            vec3(ca_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "C",
            "A",
            seq,
            vec3(c_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "O",
            "A",
            seq,
            vec3(o_x, 0.0, 0.0),
        ));
    }
    for atom in &mut atoms {
        atom.entity_id = "1".to_string();
    }
    let entities = vec![Entity {
        id: "1".to_string(),
        type_name: "polymer".to_string(),
        description: String::new(),
    }];

    let structure = Molecule {
        atoms,
        coarse_spheres: vec![CoarseSphere {
            id: 1,
            model_num: 1,
            entity_id: "1".to_string(),
            asym_id: "B".to_string(),
            seq_id_begin: 1,
            seq_id_end: 3,
            position: vec3(10.0, 0.0, 0.0),
            radius: 1.0,
            rmsf: 0.0,
        }],
        entity_index: EntityIndexMap::from_entities(&entities, &[], &[]),
        entities,
        ..Molecule::default()
    }
    .atomic_structure();

    assert!(structure.ranges.cyclic_polymer_map.is_empty());
    assert_eq!(structure.model.sequence.sequences.len(), 2);
    let entity_key = structure
        .model
        .hierarchy
        .index
        .entity_from_chain(0)
        .unwrap();
    let sequence_index = structure.model.sequence.by_entity_key[&entity_key];
    assert!(structure.model.sequence.sequences[sequence_index]
        .residues
        .is_empty());
    assert_eq!(
        structure.model.sequence.sequences[sequence_index].ranges,
        vec![SequenceRange {
            seq_id_begin: 1,
            seq_id_end: 3,
        }]
    );
}

#[test]
fn atomic_ranges_use_coarse_backbone_threshold_when_direction_atoms_are_missing() {
    let structure = Molecule {
        atoms: [0.0, 8.0, 16.0]
            .into_iter()
            .enumerate()
            .map(|(i, x)| test_atom(i + 1, "CA", "A", i as i32 + 1, vec3(x, 0.0, 0.0)))
            .collect(),
        ..Molecule::default()
    }
    .atomic_structure();

    assert_eq!(structure.ranges.polymer_ranges.len() / 2, 1);
    assert_eq!(structure.ranges.gap_ranges.len() / 2, 0);
}

#[test]
fn atomic_ranges_match_molstar_terminal_single_residue_behavior() {
    let structure = Molecule {
        atoms: vec![
            test_atom(1, "N", "A", 1, vec3(0.0, 0.0, 0.0)),
            test_atom(2, "CA", "A", 1, vec3(1.0, 0.0, 0.0)),
            test_atom(3, "C", "A", 1, vec3(2.0, 0.0, 0.0)),
        ],
        ..Molecule::default()
    }
    .atomic_structure();

    assert!(structure.ranges.polymer_ranges.is_empty());
    assert!(structure.ranges.gap_ranges.is_empty());
    assert!(structure.units[0].props.polymer_elements.is_empty());
}

#[test]
fn atomic_ranges_skip_coordinate_dependent_checks_when_xyz_is_missing_like_molstar() {
    let mut atoms = Vec::new();
    for (seq, x) in [(1, 0.0), (2, 25.0), (3, 50.0)] {
        atoms.push(test_atom(
            atoms.len() + 1,
            "CA",
            "A",
            seq,
            vec3(x, 0.0, 0.0),
        ));
    }
    let structure = Molecule {
        atom_site_columns: AtomSiteColumnPresence {
            xyz_defined: false,
            ..AtomSiteColumnPresence::default()
        },
        atoms,
        ..Molecule::default()
    }
    .atomic_structure();

    assert_eq!(structure.ranges.polymer_ranges, vec![0, 2]);
    assert!(structure.ranges.gap_ranges.is_empty());
}

#[test]
fn atomic_ranges_skip_cyclic_detection_when_xyz_is_missing_like_molstar() {
    let atoms = [0.0, 2.0, 1.0]
        .into_iter()
        .enumerate()
        .map(|(i, x)| test_atom(i + 1, "CA", "A", i as i32 + 1, vec3(x, 0.0, 0.0)))
        .collect::<Vec<_>>();
    let entities = vec![Entity {
        id: String::new(),
        type_name: "polymer".to_string(),
        description: String::new(),
    }];
    let entity_poly_seq = (1..=3)
        .map(|num| EntityPolySeq {
            entity_id: String::new(),
            num,
            mon_id: "ALA".to_string(),
            hetero: "n".to_string(),
        })
        .collect::<Vec<_>>();
    let structure = Molecule {
        atom_site_columns: AtomSiteColumnPresence {
            xyz_defined: false,
            ..AtomSiteColumnPresence::default()
        },
        atoms,
        entity_index: EntityIndexMap::from_entities(&entities, &[], &[]),
        entities,
        entity_poly_seq,
        ..Molecule::default()
    }
    .atomic_structure();

    assert!(structure.ranges.cyclic_polymer_map.is_empty());
}

#[test]
fn atomic_units_get_molstar_coarse_grained_trait_for_trace_only_models() {
    let structure = Molecule {
        atoms: [0.0, 2.0, 4.0]
            .into_iter()
            .enumerate()
            .map(|(i, x)| test_atom(i + 1, "CA", "A", i as i32 + 1, vec3(x, 0.0, 0.0)))
            .collect(),
        ..Molecule::default()
    }
    .atomic_structure();

    assert!(structure.units[0]
        .traits
        .contains(UnitTraits::COARSE_GRAINED));
}

#[test]
fn atomic_ranges_use_molstar_strict_backbone_distance_thresholds() {
    let coarse = Molecule {
        atoms: [0.0, 10.0, 20.0]
            .into_iter()
            .enumerate()
            .map(|(i, x)| test_atom(i + 1, "CA", "A", i as i32 + 1, vec3(x, 0.0, 0.0)))
            .collect(),
        ..Molecule::default()
    }
    .atomic_structure();
    assert_eq!(coarse.ranges.polymer_ranges.len() / 2, 2);
    assert_eq!(coarse.ranges.gap_ranges.len() / 2, 1);

    let mut atoms = Vec::new();
    for (seq, n_x, ca_x, c_x, o_x) in [
        (1, -2.0, -1.0, 0.0, 0.1),
        (2, 3.0, 3.5, 4.0, 4.1),
        (3, 5.0, 5.5, 6.0, 6.1),
    ] {
        atoms.push(test_atom(
            atoms.len() + 1,
            "N",
            "A",
            seq,
            vec3(n_x, 1.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "CA",
            "A",
            seq,
            vec3(ca_x, 1.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "C",
            "A",
            seq,
            vec3(c_x, 1.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "O",
            "A",
            seq,
            vec3(o_x, 1.0, 0.0),
        ));
    }
    let atomic = Molecule {
        atoms,
        ..Molecule::default()
    }
    .atomic_structure();
    assert_eq!(atomic.ranges.polymer_ranges.len() / 2, 2);
    assert_eq!(atomic.ranges.gap_ranges.len() / 2, 1);
}

#[test]
fn molstar_gamma_beta_backbone_fallback_uses_only_ca_for_connectivity() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nGAM 'l-gamma-peptide, c-delta linking'\nBET 'l-beta-peptide, c-gamma linking'\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C BB GAM A 1 0.000 0.000 0.000\nATOM 2 C BB BET A 2 1.000 0.000 0.000\nATOM 3 C BB GAM A 3 2.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = molecule.atomic_structure();

    assert_eq!(
        structure.model.hierarchy.derived.residue.polymer_type,
        vec![
            PolymerType::GammaPeptide,
            PolymerType::BetaPeptide,
            PolymerType::GammaPeptide,
        ]
    );
    assert_eq!(structure.ranges.polymer_ranges.len() / 2, 2);
    assert_eq!(structure.ranges.gap_ranges.len() / 2, 1);
    assert_eq!(structure.units[0].props.polymer_elements, vec![0, 1, 2]);
    assert_eq!(structure.units[0].props.gap_elements, vec![0, 1]);
}

#[test]
fn molstar_gamma_beta_backbone_end_roles_use_cd_and_cg() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nGAM 'l-gamma-peptide, c-delta linking'\nBET 'l-beta-peptide, c-gamma linking'\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 N N GAM G 1 -2.000 0.000 0.000\nATOM 2 C CA GAM G 1 -1.000 0.000 0.000\nATOM 3 C CD GAM G 1 0.000 0.000 0.000\nATOM 4 C C GAM G 1 20.000 0.000 0.000\nATOM 5 O O GAM G 1 20.100 0.000 0.000\nATOM 6 N N GAM G 2 1.000 0.000 0.000\nATOM 7 C CA GAM G 2 1.500 0.000 0.000\nATOM 8 C CD GAM G 2 2.000 0.000 0.000\nATOM 9 C C GAM G 2 21.000 0.000 0.000\nATOM 10 O O GAM G 2 21.100 0.000 0.000\nATOM 11 N N GAM G 3 3.000 0.000 0.000\nATOM 12 C CA GAM G 3 3.500 0.000 0.000\nATOM 13 C CD GAM G 3 4.000 0.000 0.000\nATOM 14 C C GAM G 3 22.000 0.000 0.000\nATOM 15 O O GAM G 3 22.100 0.000 0.000\nATOM 16 N N BET B 1 -2.000 1.000 0.000\nATOM 17 C CA BET B 1 -1.000 1.000 0.000\nATOM 18 C CG BET B 1 0.000 1.000 0.000\nATOM 19 C C BET B 1 20.000 1.000 0.000\nATOM 20 O O BET B 1 20.100 1.000 0.000\nATOM 21 N N BET B 2 1.000 1.000 0.000\nATOM 22 C CA BET B 2 1.500 1.000 0.000\nATOM 23 C CG BET B 2 2.000 1.000 0.000\nATOM 24 C C BET B 2 21.000 1.000 0.000\nATOM 25 O O BET B 2 21.100 1.000 0.000\nATOM 26 N N BET B 3 3.000 1.000 0.000\nATOM 27 C CA BET B 3 3.500 1.000 0.000\nATOM 28 C CG BET B 3 4.000 1.000 0.000\nATOM 29 C C BET B 3 22.000 1.000 0.000\nATOM 30 O O BET B 3 22.100 1.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = molecule.atomic_structure();

    assert_eq!(structure.ranges.polymer_ranges.len() / 2, 2);
    assert!(structure.ranges.gap_ranges.is_empty());
    assert_eq!(
        structure
            .model
            .hierarchy
            .derived
            .residue
            .polymer_type
            .iter()
            .filter(|ty| **ty == PolymerType::GammaPeptide)
            .count(),
        3
    );
    assert_eq!(
        structure
            .model
            .hierarchy
            .derived
            .residue
            .polymer_type
            .iter()
            .filter(|ty| **ty == PolymerType::BetaPeptide)
            .count(),
        3
    );
}

#[test]
fn molstar_structure_builder_merges_water_and_single_atom_chains_and_tracks_inter_unit_bonds() {
    let mut atoms = vec![
        test_atom(1, "O", "W1", 1, vec3(0.0, 0.0, 0.0)),
        test_atom(2, "O", "W2", 1, vec3(1.0, 0.0, 0.0)),
        test_atom(3, "ZN", "I1", 1, vec3(3.0, 0.0, 0.0)),
        test_atom(4, "CL", "I2", 1, vec3(4.0, 0.0, 0.0)),
        test_atom(5, "CA", "P", 1, vec3(6.0, 0.0, 0.0)),
    ];
    atoms[0].residue = "HOH".to_string();
    atoms[1].residue = "HOH".to_string();
    atoms[0].entity_id = "W".to_string();
    atoms[1].entity_id = "W".to_string();
    atoms[2].entity_id = "ION".to_string();
    atoms[3].entity_id = "ION".to_string();
    atoms[4].entity_id = "P".to_string();
    atoms[2].auth_chain = "I".to_string();
    atoms[3].auth_chain = "I".to_string();
    let entities = vec![
        Entity {
            id: "W".to_string(),
            type_name: "water".to_string(),
            description: String::new(),
        },
        Entity {
            id: "ION".to_string(),
            type_name: "non-polymer".to_string(),
            description: String::new(),
        },
        Entity {
            id: "P".to_string(),
            type_name: "polymer".to_string(),
            description: String::new(),
        },
    ];
    let molecule = Molecule {
        atoms,
        entity_index: EntityIndexMap::from_entities(&entities, &[], &[]),
        entities,
        bonds: vec![Bond { a: 0, b: 4 }, Bond { a: 2, b: 3 }],
        bond_metadata: vec![
            BondMetadata {
                source: BondSource::StructConn,
                order: 3,
                flags: BondFlags::COVALENT.union(BondFlags::DISULFIDE),
                key: 11,
                distance: Some(2.0),
                operator_a: -1,
                operator_b: -1,
                struct_conn: None,
            },
            BondMetadata::pdb_conect(12),
        ],
        ..Molecule::default()
    };

    let structure = molecule.atomic_structure();
    assert_eq!(structure.units.len(), 3);
    let water_unit = structure
        .units
        .iter()
        .find(|unit| unit.traits.contains(UnitTraits::WATER))
        .unwrap();
    assert_eq!(water_unit.chain_indices.len(), 2);
    assert_eq!(water_unit.elements, vec![0, 1]);
    assert_eq!(
        water_unit.residue_index_by_element,
        structure.model.hierarchy.residue_atom_segments.index
    );
    assert_eq!(
        water_unit.chain_index_by_element,
        structure.model.hierarchy.chain_atom_segments.index
    );
    assert!(water_unit.traits.contains(UnitTraits::MULTI_CHAIN));
    let ion_unit = structure
        .units
        .iter()
        .find(|unit| {
            unit.traits.contains(UnitTraits::MULTI_CHAIN)
                && !unit.traits.contains(UnitTraits::WATER)
        })
        .unwrap();
    assert_eq!(ion_unit.chain_indices.len(), 2);
    assert_eq!(ion_unit.elements, vec![2, 3]);
    assert_eq!(ion_unit.residue_index_by_element[2], 2);
    assert_eq!(ion_unit.residue_index_by_element[3], 3);
    assert_eq!(ion_unit.chain_index_by_element[2], 2);
    assert_eq!(ion_unit.chain_index_by_element[3], 3);
    let protein_unit = structure
        .units
        .iter()
        .find(|unit| {
            !unit.traits.contains(UnitTraits::MULTI_CHAIN)
                && !unit.traits.contains(UnitTraits::WATER)
        })
        .unwrap();
    assert_eq!(structure.intra_unit_bond_count, 1);
    assert_eq!(structure.inter_unit_bonds.len(), 1);
    assert_eq!(
        sorted_pair(
            structure.inter_unit_bonds[0].unit_a,
            structure.inter_unit_bonds[0].unit_b
        ),
        sorted_pair(water_unit.id, protein_unit.id)
    );
    assert_eq!(structure.inter_unit_bonds[0].order, 3);
    assert!(structure.inter_unit_bonds[0]
        .flags
        .contains(BondFlags::DISULFIDE));
    assert_eq!(structure.inter_unit_bonds[0].key, 11);
    assert_eq!(water_unit.props.inter_unit_bond_count, 1);
    assert_eq!(ion_unit.props.intra_unit_bond_count, 1);
    assert_eq!(protein_unit.props.inter_unit_bond_count, 1);

    let expanded = molecule.expanded_for_geometry();
    assert_eq!(expanded.bonds.len(), 2);
}

#[test]
fn molstar_single_atom_grouping_rejects_different_entity_auth_asym_operator_and_multi_atom_neighbor(
) {
    let mut different_entity = vec![
        test_atom(1, "ZN", "A", 1, vec3(0.0, 0.0, 0.0)),
        test_atom(2, "CL", "B", 1, vec3(1.0, 0.0, 0.0)),
    ];
    different_entity[0].entity_id = "E1".to_string();
    different_entity[1].entity_id = "E2".to_string();
    different_entity[0].auth_chain = "I".to_string();
    different_entity[1].auth_chain = "I".to_string();
    assert_eq!(
        Molecule {
            atoms: different_entity,
            ..Molecule::default()
        }
        .atomic_structure()
        .units
        .len(),
        2
    );

    let mut different_auth = vec![
        test_atom(1, "ZN", "A", 1, vec3(0.0, 0.0, 0.0)),
        test_atom(2, "CL", "B", 1, vec3(1.0, 0.0, 0.0)),
    ];
    for atom in &mut different_auth {
        atom.entity_id = "ION".to_string();
    }
    different_auth[0].auth_chain = "I".to_string();
    different_auth[1].auth_chain = "J".to_string();
    assert_eq!(
        Molecule {
            atoms: different_auth,
            ..Molecule::default()
        }
        .atomic_structure()
        .units
        .len(),
        2
    );

    let mut different_operator = vec![
        test_atom(1, "ZN", "A", 1, vec3(0.0, 0.0, 0.0)),
        test_atom(2, "CL", "B", 1, vec3(1.0, 0.0, 0.0)),
    ];
    for atom in &mut different_operator {
        atom.entity_id = "ION".to_string();
        atom.auth_chain = "I".to_string();
    }
    different_operator[0].operator_name = "1_555".to_string();
    different_operator[1].operator_name = "2_666".to_string();
    assert_eq!(
        Molecule {
            atoms: different_operator,
            ..Molecule::default()
        }
        .atomic_structure()
        .units
        .len(),
        2
    );

    let mut multi_atom_neighbor = vec![
        test_atom(1, "ZN", "A", 1, vec3(0.0, 0.0, 0.0)),
        test_atom(2, "CL", "B", 1, vec3(1.0, 0.0, 0.0)),
        test_atom(3, "CL", "B", 1, vec3(2.0, 0.0, 0.0)),
    ];
    for atom in &mut multi_atom_neighbor {
        atom.entity_id = "ION".to_string();
        atom.auth_chain = "I".to_string();
    }
    assert_eq!(
        Molecule {
            atoms: multi_atom_neighbor,
            ..Molecule::default()
        }
        .atomic_structure()
        .units
        .len(),
        2
    );
}

#[test]
fn assembly_filters_merged_multi_chain_units_to_selected_asym_ids() {
    let cif = b"data_demo\nloop_\n_entity.id\n_entity.type\n1 water\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 O O HOH W1 1 1 0.0 0.0 0.0\nHETATM 2 O O HOH W2 1 1 1.0 0.0 0.0\n#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n1 1 W1\n#\nloop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\n1 1 0 0 0 0 1 0 0 0 0 1 0\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: Some("1".to_string()),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = molecule.atomic_structure();

    assert_eq!(structure.units.len(), 1);
    assert_eq!(structure.units[0].chain_indices, vec![0]);
    assert_eq!(structure.units[0].elements, vec![0]);
    assert_eq!(structure.element_count, 1);
}

#[test]
fn assembly_suppresses_duplicate_single_atom_units_with_signed_zero_coordinates() {
    let mut atom = test_atom(1, "ZN", "A", 1, vec3(0.0, 0.0, 0.0));
    atom.entity_id = "ION".to_string();
    let assembly = Assembly {
        id: "1".to_string(),
        details: String::new(),
        oligomeric_details: String::new(),
        oligomeric_count: None,
        asym_ids: vec!["A".to_string()],
        transforms: Vec::new(),
        generators: vec![AssemblyGenerator::from_transforms(
            "1",
            vec!["A".to_string()],
            0,
            vec![
                Transform::identity(),
                Transform {
                    m: [
                        [1.0, 0.0, 0.0, -0.0],
                        [0.0, 1.0, 0.0, -0.0],
                        [0.0, 0.0, 1.0, -0.0],
                    ],
                },
            ],
            vec![vec!["1".to_string()], vec!["2".to_string()]],
        )],
    };

    let structure = Molecule {
        atoms: vec![atom],
        selected_assembly: Some(assembly),
        ..Molecule::default()
    }
    .atomic_structure();

    assert_eq!(structure.units.len(), 1);
    assert_eq!(structure.element_count, 1);
}

#[test]
fn molstar_structure_builder_uses_entity_type_for_water_unit_grouping() {
    let cif = b"data_demo\nloop_\n_entity.id\n_entity.type\n1 water\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 O O UNK W1 1 1 0.0 0.0 0.0\nHETATM 2 O O UNK W2 1 1 1.0 0.0 0.0\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let structure = molecule.atomic_structure();
    assert_eq!(structure.units.len(), 1);
    assert_eq!(structure.units[0].chain_indices, vec![0, 1]);
    assert!(structure.units[0].traits.contains(UnitTraits::WATER));
    assert!(structure.units[0].traits.contains(UnitTraits::MULTI_CHAIN));
}

#[test]
fn atomic_derived_data_uses_atomic_numbers_and_chem_comp_types() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nXAA 'l-peptide linking'\nNH4 ion\nRCX 'rna linking'\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA XAA A 1 0.000 0.000 0.000\nATOM 2 C C XAA A 1 1.000 0.000 0.000\nATOM 3 O O XAA A 1 1.200 0.000 0.000\nHETATM 4 Na NA NH4 B 1 4.000 0.000 0.000\nATOM 5 O \"O3'\" RCX C 1 5.000 0.000 0.000\nATOM 6 C \"C3'\" RCX C 1 5.500 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(
        molecule.chemical_components,
        vec![
            ChemicalComponent {
                id: "XAA".to_string(),
                type_name: "l-peptide linking".to_string(),
                ..ChemicalComponent::default()
            },
            ChemicalComponent {
                id: "NH4".to_string(),
                type_name: "ion".to_string(),
                ..ChemicalComponent::default()
            },
            ChemicalComponent {
                id: "RCX".to_string(),
                type_name: "rna linking".to_string(),
                ..ChemicalComponent::default()
            },
        ]
    );
    let structure = molecule.atomic_structure();
    let derived = &structure.model.hierarchy.derived;
    assert_eq!(derived.atom.atomic_number, vec![6, 6, 8, 11, 8, 6]);
    assert_eq!(
        derived.residue.molecule_type,
        vec![MoleculeType::Protein, MoleculeType::Ion, MoleculeType::Rna]
    );
    assert_eq!(
        derived.residue.polymer_type,
        vec![PolymerType::PeptideL, PolymerType::None, PolymerType::Rna]
    );
    assert_eq!(
        derived.residue.trace_element_index,
        vec![Some(0), None, Some(4)]
    );
    assert_eq!(derived.residue.direction_from_element_index[0], Some(1));
    assert_eq!(derived.residue.direction_to_element_index[0], Some(2));
    assert_eq!(derived.residue.direction_from_element_index[2], None);
    assert_eq!(derived.residue.direction_to_element_index[2], Some(5));
    assert!(structure.ranges.polymer_ranges.is_empty());
}

#[test]
fn atomic_derived_data_uses_molstar_generated_component_name_tables() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nNAG non-polymer\nXSC saccharide\nXLI lipid\nXIO ion\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 O O OHX A 1 0.000 0.000 0.000\nHETATM 2 C C1 CHL B 1 1.000 0.000 0.000\nHETATM 3 C C1 SQD C 1 2.000 0.000 0.000\nHETATM 4 O O HOH D 1 3.000 0.000 0.000\nHETATM 5 C C1 NAG E 1 4.000 0.000 0.000\nHETATM 6 C C1 XSC F 1 5.000 0.000 0.000\nHETATM 7 C C1 XLI G 1 6.000 0.000 0.000\nHETATM 8 Na NA XIO H 1 7.000 0.000 0.000\nHETATM 9 X X ZZZ I 1 8.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let structure = molecule.atomic_structure();
    let derived = &structure.model.hierarchy.derived;
    assert_eq!(
        derived.residue.molecule_type,
        vec![
            MoleculeType::Ion,
            MoleculeType::Lipid,
            MoleculeType::Saccharide,
            MoleculeType::Water,
            MoleculeType::Saccharide,
            MoleculeType::Saccharide,
            MoleculeType::Unknown,
            MoleculeType::Unknown,
            MoleculeType::Other,
        ]
    );
    assert_eq!(
        derived.residue.trace_element_index,
        vec![None, None, None, None, None, None, None, None, None]
    );
}

#[test]
fn saccharide_constants_and_component_maps_match_molstar_defaults() {
    assert_eq!(
        get_saccharide_name(SaccharideType::DiDeoxyhexose),
        "Di-deoxyhexose"
    );
    assert_eq!(
        get_saccharide_shape(SaccharideType::HexNAc, 6),
        SaccharideShape::FilledCube
    );
    assert_eq!(
        get_saccharide_shape(SaccharideType::Unknown, 4),
        SaccharideShape::DiamondPrism
    );
    assert_eq!(
        get_saccharide_shape(SaccharideType::Unknown, 8),
        SaccharideShape::FlatHexagon
    );

    let nag = saccharide_component("NAG").unwrap();
    assert_eq!(nag.abbr, "GlcNAc");
    assert_eq!(nag.name, "N-Acetyl Glucosamine");
    assert_eq!(nag.color, 0x0090bc);
    assert_eq!(nag.component_type, SaccharideType::HexNAc);

    let glc = saccharide_component("GLC").unwrap();
    assert_eq!(glc.abbr, "Glc");
    assert_eq!(glc.component_type, SaccharideType::Hexose);

    let charmm = saccharide_component("AGLC").unwrap();
    assert_eq!(charmm.abbr, "Glc");
    assert_eq!(charmm.component_type, SaccharideType::Hexose);

    let generated_only = saccharide_component("SQD").unwrap();
    assert_eq!(generated_only.abbr, "Unk");
    assert_eq!(generated_only.name, "Unknown");
    assert_eq!(generated_only.color, 0xf1ece1);
    assert_eq!(generated_only.component_type, SaccharideType::Unknown);

    assert!(saccharide_component("0GA").is_none());
    assert_eq!(
        saccharide_component_with_map("0GA", SaccharideCompIdMapType::Glycam)
            .unwrap()
            .abbr,
        "Glc"
    );
    assert_eq!(
        saccharide_component_with_map("4GL", SaccharideCompIdMapType::Glycam)
            .unwrap()
            .abbr,
        "Gul"
    );
    assert_eq!(
        saccharide_component_with_map("SIA", SaccharideCompIdMapType::Glycam)
            .unwrap()
            .abbr,
        "Neu5Ac"
    );
}

#[test]
fn carbohydrate_residue_detection_uses_default_saccharide_component_map() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nBET 'd-saccharide, beta linking'\nOLD 'l-saccharide 1,4 and 1,6 linking'\nGYC other\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 AGLC A 1 0.000 0.000 0.000\nHETATM 2 C C1 SQD B 1 1.000 0.000 0.000\nHETATM 3 C C1 BET C 1 2.000 0.000 0.000\nHETATM 4 C C1 OLD D 1 3.000 0.000 0.000\nHETATM 5 C C1 0GA E 1 4.000 0.000 0.000\nHETATM 6 C C1 GYC F 1 5.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(
        molecule
            .atomic_structure()
            .model
            .hierarchy
            .derived
            .residue
            .molecule_type,
        vec![
            MoleculeType::Saccharide,
            MoleculeType::Saccharide,
            MoleculeType::Saccharide,
            MoleculeType::Saccharide,
            MoleculeType::Other,
            MoleculeType::Other,
        ]
    );
}

#[test]
fn atomic_derived_molecule_type_priority_matches_molstar() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nSOL other\nNAG non-polymer\nSAC saccharide\nION ion\nLIP lipid\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 O O W A 1 0.000 0.000 0.000\nHETATM 2 O O SOL B 1 1.000 0.000 0.000\nHETATM 3 O O OHX C 1 2.000 0.000 0.000\nHETATM 4 N N NH4 D 1 3.000 0.000 0.000\nHETATM 5 C C1 CHL E 1 4.000 0.000 0.000\nHETATM 6 C C1 NAG F 1 5.000 0.000 0.000\nHETATM 7 C C1 SQD G 1 6.000 0.000 0.000\nHETATM 8 C C1 SAC H 1 7.000 0.000 0.000\nHETATM 9 C C1 ION I 1 8.000 0.000 0.000\nHETATM 10 C C1 LIP J 1 9.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = molecule.atomic_structure();

    assert_eq!(
        structure.model.hierarchy.derived.residue.molecule_type,
        vec![
            MoleculeType::Water,
            MoleculeType::Water,
            MoleculeType::Ion,
            MoleculeType::Ion,
            MoleculeType::Lipid,
            MoleculeType::Saccharide,
            MoleculeType::Saccharide,
            MoleculeType::Saccharide,
            MoleculeType::Unknown,
            MoleculeType::Unknown,
        ]
    );
    assert!(structure.model.hierarchy.derived.atom.is_water[0]);
    assert!(structure.model.hierarchy.derived.atom.is_water[1]);
}

#[test]
fn atomic_derived_polymer_type_classification_matches_molstar() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nPEL 'l-peptide linking'\nTER 'l-peptide cooh carboxy terminus'\nGAM 'l-gamma-peptide, c-delta linking'\nBET 'l-beta-peptide, c-gamma linking'\nRNA 'rna linking'\nDNA 'dna linking'\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA PEL A 1 0.000 0.000 0.000\nATOM 2 C C TER A 2 1.000 0.000 0.000\nATOM 3 C CA GAM B 1 2.000 0.000 0.000\nATOM 4 C CA BET C 1 3.000 0.000 0.000\nATOM 5 O \"O3'\" RNA D 1 4.000 0.000 0.000\nATOM 6 O \"O3'\" DNA E 1 5.000 0.000 0.000\nATOM 7 N \"N4'\" APN F 1 6.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let derived = molecule.atomic_structure().model.hierarchy.derived.residue;

    assert_eq!(
        derived.molecule_type,
        vec![
            MoleculeType::Protein,
            MoleculeType::Protein,
            MoleculeType::Protein,
            MoleculeType::Protein,
            MoleculeType::Rna,
            MoleculeType::Dna,
            MoleculeType::Pna,
        ]
    );
    assert_eq!(
        derived.polymer_type,
        vec![
            PolymerType::PeptideL,
            PolymerType::None,
            PolymerType::GammaPeptide,
            PolymerType::BetaPeptide,
            PolymerType::Rna,
            PolymerType::Dna,
            PolymerType::Pna,
        ]
    );
    assert_eq!(
        expand_oper_expression("(X0)(1-2)]"),
        vec![
            vec!["X0".to_string(), "1".to_string()],
            vec!["X0".to_string(), "2".to_string()]
        ]
    );
}

#[test]
fn atomic_derived_atom_role_selection_matches_molstar_atom_order() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nPRT 'l-peptide linking'\nRNA 'rna linking'\nDNA 'dna linking'\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 O OXT PRT A 1 0.000 0.000 0.000\nATOM 2 O O PRT A 1 0.100 0.000 0.000\nATOM 3 C C PRT A 1 0.200 0.000 0.000\nATOM 4 C BB PRT A 1 0.300 0.000 0.000\nATOM 5 O \"O3*\" RNA B 1 1.000 0.000 0.000\nATOM 6 O \"O3'\" RNA B 1 1.100 0.000 0.000\nATOM 7 C \"C3'\" RNA B 1 1.200 0.000 0.000\nATOM 8 C \"C4'\" RNA B 1 1.300 0.000 0.000\nATOM 9 O \"O3'\" DNA C 1 2.000 0.000 0.000\nATOM 10 C \"C1*\" DNA C 1 2.100 0.000 0.000\nATOM 11 C \"C3*\" DNA C 1 2.200 0.000 0.000\nATOM 12 N \"N4*\" APN D 1 3.000 0.000 0.000\nATOM 13 C \"C7*\" APN D 1 3.100 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = molecule.atomic_structure();
    let derived = &structure.model.hierarchy.derived.residue;

    assert_eq!(
        derived.polymer_type,
        vec![
            PolymerType::PeptideL,
            PolymerType::Rna,
            PolymerType::Dna,
            PolymerType::Pna,
        ]
    );
    assert_eq!(
        derived.trace_element_index,
        vec![Some(3), Some(4), Some(8), Some(11)]
    );
    assert_eq!(
        derived.direction_from_element_index,
        vec![Some(2), Some(7), Some(10), Some(11)]
    );
    assert_eq!(
        derived.direction_to_element_index,
        vec![Some(0), Some(6), Some(9), Some(12)]
    );
}

#[test]
fn atomic_derived_atom_role_names_are_case_sensitive_like_molstar() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nPRT 'l-peptide linking'\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 N ca PRT A 1 0.000 0.000 0.000\nATOM 2 N c PRT A 1 0.100 0.000 0.000\nATOM 3 N o PRT A 1 0.200 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = molecule.atomic_structure();
    let derived = &structure.model.hierarchy.derived.residue;

    assert_eq!(derived.molecule_type, vec![MoleculeType::Protein]);
    assert_eq!(derived.trace_element_index, vec![None]);
    assert_eq!(derived.direction_from_element_index, vec![None]);
    assert_eq!(derived.direction_to_element_index, vec![None]);
}

#[test]
fn atomic_derived_default_component_type_matches_molstar_fallback_names() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA LSN A 1 0.000 0.000 0.000\nATOM 2 C CA ASPP A 2 1.000 0.000 0.000\nATOM 3 C CA GLUP A 3 2.000 0.000 0.000\nATOM 4 C CA HID A 4 3.000 0.000 0.000\nATOM 5 C CA HIE A 5 4.000 0.000 0.000\nATOM 6 C CA HIP A 6 5.000 0.000 0.000\nATOM 7 C CA LYN A 7 6.000 0.000 0.000\nATOM 8 C CA ASH A 8 7.000 0.000 0.000\nATOM 9 C CA GLH A 9 8.000 0.000 0.000\nATOM 10 C C1 T B 1 9.000 0.000 0.000\nATOM 11 C C1 N B 2 10.000 0.000 0.000\nATOM 12 C C1 DN C 1 11.000 0.000 0.000\nHETATM 13 C C1 AGLC D 1 12.000 0.000 0.000\nATOM 14 C C1 RA E 1 13.000 0.000 0.000\nATOM 15 C C1 RC E 2 14.000 0.000 0.000\nATOM 16 C C1 RG E 3 15.000 0.000 0.000\nATOM 17 C C1 RU E 4 16.000 0.000 0.000\nATOM 18 C C1 RI E 5 17.000 0.000 0.000\nATOM 19 C C1 ADE F 1 18.000 0.000 0.000\nATOM 20 C C1 CYT F 2 19.000 0.000 0.000\nATOM 21 C C1 GUA F 3 20.000 0.000 0.000\nATOM 22 C C1 THY F 4 21.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let molecule_type = &molecule
        .atomic_structure()
        .model
        .hierarchy
        .derived
        .residue
        .molecule_type;
    assert_eq!(&molecule_type[0..9], [MoleculeType::Protein; 9]);
    assert_eq!(molecule_type[9], MoleculeType::Rna);
    assert_eq!(molecule_type[10], MoleculeType::Rna);
    assert_eq!(molecule_type[11], MoleculeType::Dna);
    assert_eq!(molecule_type[12], MoleculeType::Saccharide);
    assert_eq!(
        &molecule_type[13..],
        [
            MoleculeType::Other,
            MoleculeType::Other,
            MoleculeType::Other,
            MoleculeType::Ion,
            MoleculeType::Other,
            MoleculeType::Other,
            MoleculeType::Other,
            MoleculeType::Other,
            MoleculeType::Other,
        ]
    );
}

#[test]
fn atomic_derived_unusual_modified_residues_match_molstar_component_tables() {
    let residues = ["MSE", "SEC", "PYL", "SEP", "TPO", "PTR", "PCA", "HYP"];
    let mut cif = String::from(
        "data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\n\
         MSE non-polymer\nSEC other\nPYL non-polymer\nSEP other\n\
         TPO non-polymer\nPTR other\nPCA non-polymer\nHYP other\n\
         #\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
         _atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
         _atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
    );
    for (i, residue) in residues.iter().enumerate() {
        let seq_id = i + 1;
        let atom_id = i * 3 + 1;
        let x = i as f32;
        cif.push_str(&format!(
            "HETATM {atom_id} C CA {residue} A {seq_id} {x:.3} 0.000 0.000\n\
             HETATM {} C C {residue} A {seq_id} {:.3} 0.000 0.000\n\
             HETATM {} O O {residue} A {seq_id} {:.3} 0.000 0.000\n",
            atom_id + 1,
            x + 0.2,
            atom_id + 2,
            x + 0.4
        ));
    }
    cif.push_str("#\n");

    let molecule = parse_molecule_with_options(
        cif.as_bytes(),
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = molecule.atomic_structure();
    let derived = &structure.model.hierarchy.derived;

    assert_eq!(
        structure
            .model
            .hierarchy
            .residues
            .iter()
            .map(|residue| residue.comp_id.as_str())
            .collect::<Vec<_>>(),
        residues
    );
    assert_eq!(
        derived.residue.molecule_type,
        vec![MoleculeType::Protein; residues.len()]
    );
    assert_eq!(
        derived.residue.polymer_type,
        vec![PolymerType::PeptideL; residues.len()]
    );
    assert_eq!(
        derived.residue.trace_element_index,
        (0..residues.len()).map(|i| Some(i * 3)).collect::<Vec<_>>()
    );
    assert_eq!(
        derived.residue.direction_from_element_index,
        (0..residues.len())
            .map(|i| Some(i * 3 + 1))
            .collect::<Vec<_>>()
    );
    assert_eq!(
        derived.residue.direction_to_element_index,
        (0..residues.len())
            .map(|i| Some(i * 3 + 2))
            .collect::<Vec<_>>()
    );
    assert!(derived.atom.is_protein.iter().all(|is_protein| *is_protein));
}

#[test]
fn atomic_derived_chem_comp_lookup_is_exact_and_later_rows_win_like_molstar() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nqqq 'l-peptide linking'\nQQQ 'rna linking'\nDUP 'l-peptide linking'\nDUP 'rna linking'\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA qqq A 1 0.000 0.000 0.000\nATOM 2 O \"O3'\" QQQ B 1 1.000 0.000 0.000\nATOM 3 O \"O3'\" DUP C 1 2.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = molecule.atomic_structure();

    assert_eq!(
        molecule
            .chemical_components
            .iter()
            .map(|component| (component.id.as_str(), component.type_name.as_str()))
            .collect::<Vec<_>>(),
        vec![
            ("qqq", "l-peptide linking"),
            ("QQQ", "rna linking"),
            ("DUP", "rna linking"),
        ]
    );
    assert_eq!(
        structure.model.hierarchy.derived.residue.molecule_type,
        vec![MoleculeType::Protein, MoleculeType::Rna, MoleculeType::Rna]
    );
    assert_eq!(
        structure.model.hierarchy.derived.residue.polymer_type,
        vec![PolymerType::PeptideL, PolymerType::Rna, PolymerType::Rna]
    );
}

#[test]
fn atomic_derived_molecule_type_does_not_treat_ion_lipid_comp_type_as_molecule_type() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nLIG ion\nLPD lipid\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 LIG A 1 0.000 0.000 0.000\nHETATM 2 C C1 LPD B 1 1.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let structure = molecule.atomic_structure();
    assert_eq!(
        structure.model.hierarchy.derived.residue.molecule_type,
        vec![MoleculeType::Unknown, MoleculeType::Unknown]
    );
}

#[test]
fn parsers_preserve_model_entity_and_operator_metadata_for_structure_units() {
    let pdb = b"MODEL        2\nATOM      1  CA  GLY A   1       1.000   2.000   3.000  1.00 10.00           C\nENDMDL\nMODEL        9\nATOM      2  CA  ALA A   2       4.000   5.000   6.000  1.00 11.00           C\nENDMDL\n";
    let pdb_mol = parse_molecule(pdb, InputFormat::Pdb).unwrap();
    assert_eq!(
        pdb_mol
            .atoms
            .iter()
            .map(|atom| atom.model_num)
            .collect::<Vec<_>>(),
        vec![1, 2]
    );

    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.pdbx_PDB_model_num\nATOM 1 C CA GLY A 1 1 0.000 0.000 0.000 1.00 1\nATOM 2 C CA GLY A 1 2 1.000 0.000 0.000 1.00 2\n#\n";
    let cif_mol = parse_molecule(cif, InputFormat::Cif).unwrap();
    assert_eq!(cif_mol.atoms[0].entity_id, "1");
    assert_eq!(cif_mol.atoms[1].model_num, 2);
    let structure = cif_mol.atomic_structure();
    assert_eq!(structure.models.len(), 2);
    assert_eq!(structure.models[0].model_num, 1);
    assert_eq!(structure.models[1].model_num, 2);
    assert_eq!(structure.properties.label_entity_id, vec!["1"]);

    let ihm_cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\n_atom_site.ihm_model_id\nATOM 1 C CA GLY A 1 0.000 0.000 0.000 9 101\nATOM 2 C CA GLY A 1 1.000 0.000 0.000 9 102\n#\n";
    let ihm_mol = parse_molecule(ihm_cif, InputFormat::Cif).unwrap();
    assert_eq!(ihm_mol.atoms[0].model_num, 101);
    assert_eq!(ihm_mol.atoms[1].model_num, 102);

    let non_contiguous = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\nATOM 1 C CA GLY A 1 0.000 0.000 0.000 1\nATOM 2 C CA GLY A 1 1.000 0.000 0.000 2\nATOM 3 C CA GLY A 2 2.000 0.000 0.000 1\n#\n";
    let non_contiguous = parse_molecule(non_contiguous, InputFormat::Cif)
        .unwrap()
        .atomic_structure();
    assert_eq!(non_contiguous.models.len(), 3);
    assert_eq!(
        non_contiguous
            .models
            .iter()
            .map(|model| model.model_num)
            .collect::<Vec<_>>(),
        vec![1, 2, 1]
    );
    assert_eq!(
        non_contiguous
            .models
            .iter()
            .map(|model| model.hierarchy.atoms.len())
            .collect::<Vec<_>>(),
        vec![1, 1, 1]
    );

    let assembled = parse_molecule_with_options(
        include_bytes!("../../tests/fixtures/cif/assembly-altloc-helix.cif"),
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: Some("1".to_string()),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    assert!(
        assembled
            .atomic_structure()
            .units
            .iter()
            .any(|unit| unit.operator.name == "ASM_1"
                && unit.operator.instance_id.starts_with("ASM-"))
    );
}

#[test]
fn atomic_models_allocate_molstar_shaped_model_and_conformation_ids() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\nATOM 1 C CA GLY A 1 1 0.000 0.000 0.000 1\nATOM 2 C CA GLY A 1 1 1.000 0.000 0.000 2\n#\n";
    let structure = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap()
    .atomic_structure();

    assert_eq!(structure.models.len(), 2);
    assert_eq!(structure.model.id, structure.models[0].id);
    assert_eq!(
        structure.model.conformation.id,
        structure.models[0].conformation.id
    );

    let ids = structure
        .models
        .iter()
        .flat_map(|model| [model.id.as_str(), model.conformation.id.as_str()])
        .collect::<Vec<_>>();
    assert!(ids.iter().all(|id| is_molstar_uuid22(id)));
    assert_eq!(
        ids.iter()
            .copied()
            .collect::<std::collections::BTreeSet<_>>()
            .len(),
        ids.len()
    );
}

fn is_molstar_uuid22(id: &str) -> bool {
    id.len() == 22
        && id
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || byte == b'-' || byte == b'_')
}

#[test]
fn atomic_conformation_preserves_molstar_column_definedness() {
    let without_b = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 1.000 2.000 3.000\n#\n";
    let structure = parse_molecule(without_b, InputFormat::Cif)
        .unwrap()
        .atomic_structure();
    assert_eq!(structure.model.conformation.b_iso, vec![0.0]);
    assert!(structure.model.conformation.occupancy_defined);
    assert!(!structure.model.conformation.b_iso_defined);
    assert!(structure.model.conformation.xyz_defined);

    let with_b = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.B_iso_or_equiv\nATOM 1 C CA GLY A 1 1.000 2.000 3.000 12.50\n#\n";
    let structure = parse_molecule(with_b, InputFormat::Cif)
        .unwrap()
        .atomic_structure();
    assert_eq!(structure.model.conformation.b_iso, vec![12.5]);
    assert!(structure.model.conformation.b_iso_defined);
}

#[test]
fn pdb_anisou_rows_create_molstar_scaled_anisotropic_mapping() {
    let pdb = b"ATOM      1  CA  GLY A   1       1.000   2.000   3.000  1.00 10.00           C\nANISOU    1  CA  GLY A   1    1000   2000   3000    100    200    300       C\nEND\n";
    let mol = parse_molecule(pdb, InputFormat::Pdb).unwrap();
    assert_eq!(mol.atom_site_anisotrop.len(), 1);
    assert_eq!(mol.atom_site_anisotrop[0].atom_id, 1);
    assert_eq!(
        mol.atom_site_anisotrop[0].u,
        [[0.1, 0.01, 0.02], [0.01, 0.2, 0.03], [0.02, 0.03, 0.3]]
    );

    let structure = mol.atomic_structure();
    assert_eq!(structure.model.conformation.element_to_anisotrop, vec![0]);
    assert_eq!(
        structure.model.conformation.anisotropic_displacement,
        vec![Some([
            [0.1, 0.01, 0.02],
            [0.01, 0.2, 0.03],
            [0.02, 0.03, 0.3]
        ])]
    );
}

#[test]
fn atomic_hierarchy_keeps_microheterogeneous_residues_in_one_offset() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.auth_atom_id\n_atom_site.label_comp_id\n_atom_site.auth_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA CA GLY GLY A X 1 10 110 A 0.000 0.000 0.000\nATOM 2 C CB CB ALA ALA A X 1 10 110 A 1.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let structure = molecule.atomic_structure();
    assert_eq!(structure.model.hierarchy.residues.len(), 1);
    assert_eq!(structure.model.hierarchy.residue_atom_segments.count, 1);
    assert_eq!(
        structure.model.hierarchy.residue_atom_segments.offsets,
        vec![0, 2]
    );
    assert_eq!(structure.model.hierarchy.atoms[0].residue_index, 0);
    assert_eq!(structure.model.hierarchy.atoms[1].residue_index, 0);
    assert_eq!(structure.model.hierarchy.atoms[0].label_comp_id, "GLY");
    assert_eq!(structure.model.hierarchy.atoms[0].auth_comp_id, "GLY");
    assert_eq!(structure.model.hierarchy.atoms[1].label_comp_id, "ALA");
    assert_eq!(structure.model.hierarchy.atoms[1].auth_comp_id, "ALA");
    assert_eq!(structure.model.hierarchy.residues[0].comp_id, "GLY");
    assert_eq!(structure.model.hierarchy.residues[0].auth_comp_id, "GLY");
    assert_eq!(structure.model.hierarchy.residues[0].group_pdb, "ATOM");
    assert_eq!(structure.properties.residue_index, vec![0, 0]);
    assert_eq!(structure.properties.group_pdb, vec!["ATOM", "ATOM"]);
    assert_eq!(structure.properties.label_comp_id, vec!["GLY", "ALA"]);
    assert_eq!(structure.properties.auth_comp_id, vec!["GLY", "ALA"]);
}

#[test]
fn atomic_residue_group_pdb_uses_first_atom_like_molstar_table_view() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 LIG A 1 0.000 0.000 0.000\nATOM 2 C C2 LIG A 1 1.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let structure = molecule.atomic_structure();
    assert_eq!(structure.model.hierarchy.residues.len(), 1);
    assert_eq!(structure.model.hierarchy.residues[0].group_pdb, "HETATM");
    assert_eq!(structure.properties.group_pdb, vec!["HETATM", "HETATM"]);
}

#[test]
fn atomic_models_have_molstar_shaped_custom_property_placeholders() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\nATOM 1 C CA GLY A 1 0.000 0.000 0.000 1\nATOM 2 C CA GLY A 2 1.000 0.000 0.000 2\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let mut structure = molecule.atomic_structure();
    assert_eq!(structure.models.len(), 2);
    assert!(structure.model.custom_properties.is_empty());
    assert!(structure.model.static_property_data.is_empty());
    assert!(structure.model.dynamic_property_data.is_empty());
    assert!(structure.models[0].custom_properties.is_empty());
    assert!(structure.models[1].custom_properties.is_empty());

    structure.models[0]
        .custom_properties
        .register("molstar.example.static");
    structure.models[0]
        .static_property_data
        .insert("molstar.example.static", "first-model");
    structure.models[0]
        .dynamic_property_data
        .insert("molstar.example.dynamic", "frame-1");

    assert!(structure.models[0]
        .custom_properties
        .contains("molstar.example.static"));
    assert_eq!(
        structure.models[0]
            .static_property_data
            .get("molstar.example.static"),
        Some("first-model")
    );
    assert_eq!(
        structure.models[0]
            .dynamic_property_data
            .get("molstar.example.dynamic"),
        Some("frame-1")
    );
    assert!(!structure.models[1]
        .custom_properties
        .contains("molstar.example.static"));
    assert!(!structure.models[1]
        .static_property_data
        .contains_key("molstar.example.static"));
    assert!(!structure.models[1]
        .dynamic_property_data
        .contains_key("molstar.example.dynamic"));
}

#[test]
fn molstar_global_model_transform_info_is_parsed_and_attached() {
    let cif = b"data_demo\nloop_\n_molstar_global_model_transform_info.matrix[1][1]\n_molstar_global_model_transform_info.matrix[1][2]\n_molstar_global_model_transform_info.matrix[1][3]\n_molstar_global_model_transform_info.matrix[1][4]\n_molstar_global_model_transform_info.matrix[2][1]\n_molstar_global_model_transform_info.matrix[2][2]\n_molstar_global_model_transform_info.matrix[2][3]\n_molstar_global_model_transform_info.matrix[2][4]\n_molstar_global_model_transform_info.matrix[3][1]\n_molstar_global_model_transform_info.matrix[3][2]\n_molstar_global_model_transform_info.matrix[3][3]\n_molstar_global_model_transform_info.matrix[3][4]\n_molstar_global_model_transform_info.matrix[4][1]\n_molstar_global_model_transform_info.matrix[4][2]\n_molstar_global_model_transform_info.matrix[4][3]\n_molstar_global_model_transform_info.matrix[4][4]\n1.0 0.0 0.0 10.0 0.0 2.0 0.0 20.0 0.0 0.0 3.0 30.0 0.0 0.0 0.0 1.0\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\nATOM 1 C CA GLY A 1 0.000 0.000 0.000 1\nATOM 2 C CA GLY A 1 1.000 0.000 0.000 2\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let expected = [
        [1.0, 0.0, 0.0, 10.0],
        [0.0, 2.0, 0.0, 20.0],
        [0.0, 0.0, 3.0, 30.0],
        [0.0, 0.0, 0.0, 1.0],
    ];

    assert_eq!(
        molecule.global_model_transform.as_ref().unwrap().matrix,
        expected
    );
    let structure = molecule.atomic_structure();
    assert_eq!(
        structure
            .model
            .global_model_transform
            .as_ref()
            .unwrap()
            .matrix,
        expected
    );
    assert!(structure
        .model
        .custom_properties
        .contains(GlobalModelTransform::DESCRIPTOR));
    assert_eq!(
        structure
            .model
            .static_property_data
            .get(GlobalModelTransform::DESCRIPTOR),
        Some("1,0,0,10;0,2,0,20;0,0,3,30;0,0,0,1")
    );
    assert_eq!(
        structure.models[1]
            .global_model_transform
            .as_ref()
            .unwrap()
            .matrix,
        expected
    );
    assert!(structure.models[1]
        .custom_properties
        .contains(GlobalModelTransform::DESCRIPTOR));
}

#[test]
fn model_source_data_preserves_molstar_format_and_category_references() {
    let cif = include_bytes!("../../tests/fixtures/cif/water.cif");
    let cif_molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    assert_eq!(cif_molecule.source_data.kind, "mmCIF");
    assert_eq!(cif_molecule.source_data.name, "water");
    assert!(cif_molecule
        .source_data
        .categories
        .iter()
        .any(|category| category.name == "atom_site" && category.row_count == 3));
    assert!(cif_molecule
        .source_data
        .db_categories
        .iter()
        .any(|category| category.name == "struct_conn" && category.row_count == 2));
    assert_eq!(
        cif_molecule.source_data.frame_categories,
        cif_molecule.source_data.categories
    );
    assert_eq!(
        cif_molecule.atomic_structure().model.source_data,
        cif_molecule.source_data
    );

    let bcif = include_bytes!("../../tests/fixtures/bcif/water.bcif");
    let bcif_molecule = parse_molecule_with_options(
        bcif,
        &MeshOptions {
            format: InputFormat::BinaryCif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    assert_eq!(bcif_molecule.source_data.kind, "mmCIF");
    assert_eq!(bcif_molecule.source_data.name, "water");
    assert!(bcif_molecule
        .source_data
        .categories
        .iter()
        .any(|category| category.name == "atom_site" && category.column_count >= 10));

    let info = String::from_utf8(molecule_info(cif, br#"{"format":"cif"}"#).unwrap()).unwrap();
    assert!(info.contains(r#""source_data":{"kind":"mmCIF","name":"water""#));
    assert!(info.contains(r#""categories":[{"name":"entry""#));
    assert!(info.contains(r#""db":{"category_count":"#));
    assert!(info.contains(r#""frame":{"category_count":"#));
    assert!(info.contains(r#""name":"atom_site","row_count":3"#));

    let custom = b"data_custom\nloop_\n_molfig_extra.id\n_molfig_extra.value\n1 test\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 0.000 0.000 0.000\n#\n";
    let custom = parse_molecule_with_options(
        custom,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    assert!(custom
        .source_data
        .frame_categories
        .iter()
        .any(|category| category.name == "molfig_extra"));
    assert!(!custom
        .source_data
        .db_categories
        .iter()
        .any(|category| category.name == "molfig_extra"));
}

#[test]
fn binary_cif_counterparts_match_text_fixture_structure_counts() {
    let fixtures = [
        "assembly-altloc-helix",
        "assembly-altloc-secondary",
        "assembly-operator-matrix",
        "atomic-protein-altloc-tie",
        "atomic-protein-no-altloc",
        "carbohydrate-branched",
        "covalent-cross-link",
        "ihm-gaussian-only",
        "ihm-sphere-only",
        "ligand-metal-aromatic",
        "mixed-atomic-coarse-ihm",
        "mixed-protein-nucleic",
        "multi-model-atomic",
        "nucleic-acid-rna-dna",
        "tiny-peptide",
        "water",
    ];
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));

    for fixture in fixtures {
        let cif = std::fs::read(root.join(format!("tests/fixtures/cif/{fixture}.cif"))).unwrap();
        let bcif = std::fs::read(root.join(format!("tests/fixtures/bcif/{fixture}.bcif"))).unwrap();
        let options = |format| MeshOptions {
            format,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        };
        let text = parse_molecule_with_options(&cif, &options(InputFormat::Cif)).unwrap();
        let binary = parse_molecule_with_options(&bcif, &options(InputFormat::BinaryCif)).unwrap();

        assert_eq!(binary.atoms.len(), text.atoms.len(), "{fixture} atom count");
        assert_eq!(binary.bonds.len(), text.bonds.len(), "{fixture} bond count");
        assert_eq!(
            binary.coarse_spheres.len(),
            text.coarse_spheres.len(),
            "{fixture} sphere count"
        );
        assert_eq!(
            binary.coarse_gaussians.len(),
            text.coarse_gaussians.len(),
            "{fixture} gaussian count"
        );
        for text_category in &text.source_data.categories {
            let binary_category = binary
                .source_data
                .categories
                .iter()
                .find(|category| category.name == text_category.name)
                .unwrap_or_else(|| {
                    panic!(
                        "{fixture} missing BinaryCIF category {}",
                        text_category.name
                    )
                });
            assert_eq!(
                binary_category.row_count, text_category.row_count,
                "{fixture} {} row count",
                text_category.name
            );
            assert!(
                binary_category.column_count >= text_category.column_count,
                "{fixture} {} column count",
                text_category.name
            );
        }
        for (text_atom, binary_atom) in text.atoms.iter().zip(&binary.atoms) {
            assert_eq!(binary_atom.id, text_atom.id, "{fixture} atom id");
            assert_eq!(binary_atom.name, text_atom.name, "{fixture} atom name");
            assert_eq!(binary_atom.residue, text_atom.residue, "{fixture} residue");
            assert_eq!(binary_atom.chain, text_atom.chain, "{fixture} chain");
            assert_eq!(
                binary_atom.position, text_atom.position,
                "{fixture} atom position"
            );
        }
    }
}

#[test]
fn every_text_cif_fixture_has_binary_cif_counterpart() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let cif_dir = root.join("tests/fixtures/cif");
    let bcif_dir = root.join("tests/fixtures/bcif");
    let mut missing = Vec::new();

    for entry in std::fs::read_dir(&cif_dir).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().and_then(|ext| ext.to_str()) != Some("cif") {
            continue;
        }
        let stem = path.file_stem().and_then(|stem| stem.to_str()).unwrap();
        if !bcif_dir.join(format!("{stem}.bcif")).is_file() {
            missing.push(stem.to_string());
        }
    }

    assert!(
        missing.is_empty(),
        "missing BinaryCIF counterparts for text mmCIF fixtures: {missing:?}"
    );
}

#[test]
fn entity_table_mapping_synthesizes_validates_and_indexes_ids() {
    let synthesized = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 1 0.000 0.000 0.000\nATOM 2 C CA ALA A 2 1 1.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        synthesized,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    assert_eq!(
        molecule
            .entities
            .iter()
            .map(|entity| entity.id.as_str())
            .collect::<Vec<_>>(),
        vec!["1", "2"]
    );
    assert_eq!(molecule.entity_index.get_entity_index("1"), Some(0));
    assert_eq!(molecule.entity_index.get_entity_index("2"), Some(1));

    let structure = molecule.atomic_structure();
    assert_eq!(structure.model.hierarchy.chains.len(), 2);
    assert_eq!(
        structure.model.hierarchy.index.chain_entity_index,
        vec![Some(0), Some(1)]
    );
    assert_eq!(
        structure
            .model
            .hierarchy
            .index
            .chain_by_entity_and_label_asym_id(0, "A"),
        Some(0)
    );
    assert_eq!(
        structure
            .model
            .hierarchy
            .index
            .chain_by_entity_and_label_asym_id(1, "A"),
        Some(1)
    );

    let struct_asym_fallback = b"data_demo\nloop_\n_struct_asym.id\n_struct_asym.entity_id\nA 42\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 0.000 0.000 0.000\n#\n";
    let fallback = parse_molecule_with_options(
        struct_asym_fallback,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    assert_eq!(fallback.atoms[0].entity_id, "42");
    assert_eq!(fallback.entities[0].id, "42");
    assert_eq!(fallback.entity_index.get_entity_index("42"), Some(0));

    let missing = b"data_demo\nloop_\n_entity.id\n_entity.type\n1 polymer\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 2 1 0.000 0.000 0.000\n#\n";
    let error = parse_molecule_with_options(
        missing,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap_err();
    assert_eq!(error, "missing _entity row for referenced entity id 2");
}

#[test]
fn unit_kind_discriminants_match_molstar_const_enum() {
    assert_eq!(AtomicUnitKind::Atomic as u8, 0);
    assert_eq!(AtomicUnitKind::Spheres as u8, 1);
    assert_eq!(AtomicUnitKind::Gaussians as u8, 2);
}

#[test]
fn atomic_index_finds_entity_chain_residue_and_atoms_like_molstar() {
    let cif = b"data_demo\nloop_\n_entity.id\n_entity.type\n1 polymer\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.auth_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.auth_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA CAA A GLY GLY A X 1 10 101 A 0.000 0.000 0.000\nATOM 2 C CA CAB B GLY GLY A X 1 10 101 A 1.000 0.000 0.000\nATOM 3 C CB CBB . GLY GLY A X 1 10 101 A 2.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            alt_loc: "all".to_string(),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = molecule.atomic_structure();
    let hierarchy = &structure.model.hierarchy;
    let index = &hierarchy.index;

    assert_eq!(index.entity_by_label_asym_id("A"), Some(0));
    assert_eq!(index.chain_by_label_asym_id("A"), Some(0));
    assert_eq!(index.chain_by_auth_asym_and_seq_id("X", "101"), Some(0));
    assert_eq!(index.residue_by_label_key(0, "10", "A"), Some(0));
    assert_eq!(index.residue_by_label_key(0, "10", ""), None);
    assert_eq!(index.residue_by_auth_key(0, "101", "A"), Some(0));
    assert_eq!(
        index.residue_by_entity_label_asym_and_auth_seq(0, "A", "101", "A"),
        Some(0)
    );
    assert_eq!(
        index.atom_by_label_key(hierarchy, 0, "10", "A", "CA", None),
        Some(0)
    );
    assert_eq!(
        index.atom_by_label_key(hierarchy, 0, "10", "A", "CA", Some("B")),
        Some(1)
    );
    assert_eq!(
        index.atom_by_label_key(hierarchy, 0, "10", "A", "CA", Some("C")),
        None
    );
    assert_eq!(
        index.atom_by_auth_key(hierarchy, "X", "101", "A", "CAB", Some("B")),
        Some(1)
    );

    let properties = &structure.properties;
    assert_eq!(properties.atom_label_atom_id(0), Some("CA"));
    assert_eq!(properties.atom_auth_atom_id(1), Some("CAB"));
    assert_eq!(properties.atom_label_alt_id(1), Some("B"));
    assert_eq!(properties.atom_label_comp_id(0), Some("GLY"));
    assert_eq!(properties.atom_auth_comp_id(0), Some("GLY"));
    assert_eq!(properties.residue_label_comp_id(2), Some("GLY"));
    assert_eq!(properties.residue_auth_comp_id(2), Some("GLY"));
    assert_eq!(properties.residue_label_seq_id(0), Some("10"));
    assert_eq!(properties.residue_auth_seq_id(0), Some("101"));
    assert_eq!(properties.residue_pdb_ins_code(0), Some("A"));
    assert_eq!(properties.chain_label_asym_id(0), Some("A"));
    assert_eq!(properties.chain_auth_asym_id(0), Some("X"));
    assert_eq!(properties.chain_label_entity_id(0), Some("1"));
    assert_eq!(properties.atom_label_atom_id(99), None);
    assert_eq!(properties.chain_auth_asym_id(99), None);
}

#[test]
fn sequence_mapping_falls_back_to_atomic_polymer_hierarchy() {
    let cif = b"data_demo\nloop_\n_entity.id\n_entity.type\n1 polymer\n2 non-polymer\n#\nloop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\nL 2\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A 5 0.000 0.000 0.000\nATOM 2 C CA GLY A 6 1.000 0.000 0.000\nHETATM 3 C C1 LIG L 1 5.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = molecule.atomic_structure();

    assert_eq!(structure.model.sequence.sequences.len(), 1);
    assert_eq!(structure.model.sequence.sequences[0].entity_id, "1");
    assert_eq!(
        structure.model.sequence.sequences[0]
            .residues
            .iter()
            .map(|residue| (residue.comp_id.as_str(), residue.seq_id))
            .collect::<Vec<_>>(),
        vec![("ALA", 5), ("GLY", 6)]
    );
    assert!(structure
        .model
        .sequence
        .by_entity_key
        .contains_key(&molecule.entity_index.get_entity_index("1").unwrap()));
    assert!(!structure
        .model
        .sequence
        .by_entity_key
        .contains_key(&molecule.entity_index.get_entity_index("2").unwrap()));
}

#[test]
fn sequence_mapping_collapses_entity_poly_seq_microheterogeneity() {
    let cif = b"data_demo\nloop_\n_entity.id\n_entity.type\n1 polymer\n#\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n_entity_poly_seq.hetero\n1 1 ALA n\n1 2 MSE y\n1 2 MET y\n1 3 GLY n\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A 1 1 0.000 0.000 0.000\nATOM 2 C CA MSE A 1 2 1.000 0.000 0.000\nATOM 3 C CA GLY A 1 3 2.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = molecule.atomic_structure();
    let sequence = &structure.model.sequence.sequences[0];

    assert_eq!(
        sequence
            .residues
            .iter()
            .map(|residue| (residue.comp_id.as_str(), residue.seq_id))
            .collect::<Vec<_>>(),
        vec![("ALA", 1), ("MSE", 2), ("GLY", 3)]
    );
    assert_eq!(sequence.index_by_seq_id.get(&2), Some(&1));
    assert_eq!(
        sequence.micro_het.get(&2).cloned().unwrap_or_default(),
        vec!["MSE".to_string(), "MET".to_string()]
    );
}

#[test]
fn cif_atom_site_normalizes_label_auth_pairs_and_type_symbols_like_molstar() {
    let fallback = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.auth_atom_id\n_atom_site.label_comp_id\n_atom_site.auth_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 ? . FE1 . HEM ? X . 42 0.000 0.000 0.000\n#\n";
    let fallback_mol = parse_molecule(fallback, InputFormat::Cif).unwrap();
    assert_eq!(fallback_mol.atoms[0].name, "FE1");
    assert_eq!(fallback_mol.atoms[0].auth_name, "FE1");
    assert_eq!(fallback_mol.atoms[0].residue, "HEM");
    assert_eq!(fallback_mol.atoms[0].auth_residue, "HEM");
    assert_eq!(fallback_mol.atoms[0].chain, "X");
    assert_eq!(fallback_mol.atoms[0].auth_chain, "X");
    assert_eq!(fallback_mol.atoms[0].entity_id, "X");
    assert_eq!(fallback_mol.atoms[0].residue_seq, "42");
    assert_eq!(fallback_mol.atoms[0].auth_residue_seq, "42");
    assert_eq!(fallback_mol.atoms[0].type_symbol, "FE");
    assert_eq!(fallback_mol.atoms[0].element, "Fe");
    let fallback_structure = fallback_mol.atomic_structure();
    assert_eq!(fallback_structure.properties.type_symbol, vec!["FE"]);
    assert_eq!(fallback_structure.properties.auth_comp_id, vec!["HEM"]);
    assert_eq!(fallback_structure.properties.label_entity_id, vec!["X"]);

    let symbols = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 fe FE1 HEM A 1 0.000 0.000 0.000\nATOM 2 CL CL1 LIG A 2 1.000 0.000 0.000\nATOM 3 C1 C1 LIG A 3 2.000 0.000 0.000\n#\n";
    let symbols = parse_molecule(symbols, InputFormat::Cif).unwrap();
    let structure = symbols.atomic_structure();
    assert_eq!(structure.properties.type_symbol, vec!["FE", "CL", "C1"]);
    assert_eq!(
        symbols
            .atoms
            .iter()
            .map(|atom| atom.element.as_str())
            .collect::<Vec<_>>(),
        vec!["Fe", "Cl", "C"]
    );
}

#[test]
fn cif_atom_site_label_auth_fallback_is_column_level_like_molstar() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.auth_atom_id\n_atom_site.label_comp_id\n_atom_site.auth_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ? GLY ? A ? 1 ? 0.000 0.000 0.000\nATOM 2 C ? CAX ? GLA ? AA ? 101 1.000 0.000 0.000\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    assert_eq!(mol.atoms[0].name, "CA");
    assert_eq!(mol.atoms[0].auth_name, "");
    assert_eq!(mol.atoms[0].residue, "GLY");
    assert_eq!(mol.atoms[0].auth_residue, "");
    assert_eq!(mol.atoms[0].chain, "A");
    assert_eq!(mol.atoms[0].auth_chain, "");
    assert_eq!(mol.atoms[0].residue_seq, "1");
    assert_eq!(mol.atoms[0].auth_residue_seq, "");
    assert_eq!(mol.atoms[1].name, "");
    assert_eq!(mol.atoms[1].auth_name, "CAX");
    assert_eq!(mol.atoms[1].residue, "");
    assert_eq!(mol.atoms[1].auth_residue, "GLA");
    assert_eq!(mol.atoms[1].chain, "");
    assert_eq!(mol.atoms[1].auth_chain, "AA");
    assert_eq!(mol.atoms[1].residue_seq, "");
    assert_eq!(mol.atoms[1].auth_residue_seq, "101");

    let structure = mol.atomic_structure();
    assert_eq!(
        structure.properties.label_atom_id,
        vec!["CA".to_string(), "".to_string()]
    );
    assert_eq!(
        structure.properties.auth_atom_id,
        vec!["".to_string(), "CAX".to_string()]
    );
}

#[test]
fn cif_atom_site_empty_alt_ids_are_not_counted_as_alt_locs() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA . GLY A 1 0.000 0.000 0.000\nATOM 2 C CB ? GLY A 1 1.000 0.000 0.000\nATOM 3 C CG A GLY A 1 2.000 0.000 0.000\n#\n";
    let mol = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            alt_loc: "all".to_string(),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    assert_eq!(
        mol.atoms
            .iter()
            .map(|atom| atom.alt_id.as_str())
            .collect::<Vec<_>>(),
        vec!["", "", "A"]
    );
    let structure = mol.atomic_structure();
    assert_eq!(structure.alt_loc_count(), 1);
    assert_eq!(
        structure.properties.label_alt_id,
        vec!["".to_string(), "".to_string(), "A".to_string()]
    );
}

#[test]
fn cartoon_summary_marks_terminal_coils_with_boundary_residues() {
    let molecule = Molecule {
        atoms: (1..=5)
            .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
            .collect(),
        helices: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 2,
            start_insertion_code: String::new(),
            end: 4,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""geometry_type":"tube","visual":"polymer-trace","representation":"cartoon","secondary_type":"helix","chain":"A","residue_start":2,"residue_end":2,"group_id":1"#));
    assert!(summary.contains(r#""geometry_type":"tube","visual":"polymer-trace","representation":"cartoon","secondary_type":"helix","chain":"A","residue_start":4,"residue_end":4,"group_id":3"#));
    assert!(summary.contains(r#""geometry_type":"tube","visual":"polymer-trace","representation":"cartoon","secondary_type":"coil","chain":"A","residue_start":1,"residue_end":1,"group_id":0"#));
    assert!(summary.contains(r#""geometry_type":"tube","visual":"polymer-trace","representation":"cartoon","secondary_type":"coil","chain":"A","residue_start":5,"residue_end":5,"group_id":4"#));
    assert!(summary.contains(r#""value_cell":{"group_id":0"#));
    assert!(summary.contains(r#""valueCell":{"drawCount":"#));
    assert!(summary.contains(r#""uVertexCount":"#));
    assert!(summary.contains(r#""uGroupCount":1"#));
    assert!(summary.contains(r#""polymer_trace":{"initial":false,"final":false,"sec_struc_first":true,"sec_struc_last":false}"#));
    assert!(summary.contains(r#""polymer_trace":{"initial":true,"final":false,"sec_struc_first":false,"sec_struc_last":true}"#));
}

#[test]
fn cartoon_secondary_trace_flags_match_molstar_transitions() {
    let molecule = Molecule {
        atoms: (1..=6)
            .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
            .collect(),
        helices: vec![
            SecondaryRange {
                chain: "A".to_string(),
                start: 1,
                start_insertion_code: String::new(),
                end: 2,
                end_insertion_code: String::new(),
            },
            SecondaryRange {
                chain: "A".to_string(),
                start: 3,
                start_insertion_code: String::new(),
                end: 4,
                end_insertion_code: String::new(),
            },
        ],
        sheets: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 5,
            start_insertion_code: String::new(),
            end: 6,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""secondary_type":"helix","chain":"A","residue_start":1,"residue_end":1,"group_id":0,"polymer_trace":{"initial":true,"final":false,"sec_struc_first":false,"sec_struc_last":false}"#));
    assert!(summary.contains(r#""secondary_type":"helix","chain":"A","residue_start":4,"residue_end":4,"group_id":3,"polymer_trace":{"initial":false,"final":false,"sec_struc_first":false,"sec_struc_last":true}"#));
    assert!(summary.contains(r#""secondary_type":"sheet","chain":"A","residue_start":5,"residue_end":5,"group_id":4,"polymer_trace":{"initial":false,"final":false,"sec_struc_first":true,"sec_struc_last":false}"#));
    assert!(summary.contains(r#""secondary_type":"sheet","chain":"A","residue_start":6,"residue_end":6,"group_id":5,"polymer_trace":{"initial":false,"final":true,"sec_struc_first":false,"sec_struc_last":false}"#));
}

#[test]
fn molstar_terminal_caps_follow_trace_and_secondary_boundaries() {
    let molecule = Molecule {
        atoms: (1..=6)
            .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
            .collect(),
        helices: vec![
            SecondaryRange {
                chain: "A".to_string(),
                start: 1,
                start_insertion_code: String::new(),
                end: 2,
                end_insertion_code: String::new(),
            },
            SecondaryRange {
                chain: "A".to_string(),
                start: 3,
                start_insertion_code: String::new(),
                end: 4,
                end_insertion_code: String::new(),
            },
        ],
        sheets: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 5,
            start_insertion_code: String::new(),
            end: 6,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        tubular_helices: true,
        round_cap: true,
        helix_profile: PolymerProfile::Square,
        ..MeshOptions::default()
    };

    let objects = build_render_objects(&molecule, &options);
    let helix_caps = objects
        .iter()
        .filter_map(|object| {
            if let RenderObject::PolymerTraceSegment {
                kind: PolymerTraceSegmentKind::Tube { profile, round_cap },
                start_cap,
                end_cap,
                ..
            } = object
            {
                Some((*profile, *start_cap, *end_cap, *round_cap))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    let sheet_caps = objects
        .iter()
        .filter_map(|object| {
            if let RenderObject::PolymerTraceSegment {
                kind: PolymerTraceSegmentKind::Sheet { .. },
                start_cap,
                end_cap,
                ..
            } = object
            {
                Some((*start_cap, *end_cap))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    assert_eq!(
        helix_caps,
        vec![
            (PolymerProfile::Elliptical, true, false, true),
            (PolymerProfile::Elliptical, false, false, true),
            (PolymerProfile::Elliptical, false, false, true),
            (PolymerProfile::Elliptical, false, true, true),
        ]
    );
    assert_eq!(sheet_caps, vec![(true, false), (false, true)]);
}

#[test]
fn molstar_round_cap_only_applies_to_tubular_helices() {
    let molecule = Molecule {
        atoms: (1..=2)
            .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
            .collect(),
        helices: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 1,
            start_insertion_code: String::new(),
            end: 2,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        round_cap: true,
        tubular_helices: false,
        helix_profile: PolymerProfile::Rounded,
        ..MeshOptions::default()
    };

    let objects = build_render_objects(&molecule, &options);
    let RenderObject::PolymerTraceSegment {
        kind: PolymerTraceSegmentKind::Tube { profile, round_cap },
        ..
    } = &objects[0]
    else {
        panic!("expected helix trace segment");
    };
    assert_eq!(*profile, PolymerProfile::Rounded);
    assert!(!*round_cap);
}

#[test]
fn molstar_trace_uses_elliptical_tube_before_square_profile_when_width_equals_height() {
    let molecule = Molecule {
        atoms: (1..=3)
            .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
            .collect(),
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        helix_profile: PolymerProfile::Square,
        radial_segments: 20,
        ..MeshOptions::default()
    };

    let objects = build_render_objects(&molecule, &options);
    let RenderObject::PolymerTraceSegment {
        widths,
        heights,
        kind: PolymerTraceSegmentKind::Tube { profile, .. },
        ..
    } = &objects[0]
    else {
        panic!("expected coil trace tube");
    };
    assert_eq!(widths[1], heights[1]);
    assert_eq!(*profile, PolymerProfile::Elliptical);
}

#[test]
fn molstar_square_profile_still_routes_non_tubular_helix_to_sheet() {
    let molecule = Molecule {
        atoms: (1..=3)
            .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
            .collect(),
        helices: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 1,
            start_insertion_code: String::new(),
            end: 3,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        tubular_helices: false,
        helix_profile: PolymerProfile::Square,
        radial_segments: 20,
        ..MeshOptions::default()
    };

    let objects = build_render_objects(&molecule, &options);
    assert!(objects.iter().any(|object| matches!(
        object,
        RenderObject::PolymerTraceSegment {
            kind: PolymerTraceSegmentKind::Sheet { .. },
            ..
        }
    )));
    assert!(!objects.iter().any(|object| matches!(
        object,
        RenderObject::PolymerTraceSegment {
            kind: PolymerTraceSegmentKind::Tube { .. },
            ..
        }
    )));
}

#[test]
fn molstar_sec_struc_flags_clamp_at_missing_trace_polymer_ranges() {
    let atoms = vec![
        test_atom(1, "CA", "A", 1, vec3(0.0, 0.0, 0.0)),
        test_atom(2, "N", "A", 2, vec3(1.0, 0.0, 0.0)),
        test_atom(3, "CA", "A", 3, vec3(4.0, 0.0, 0.0)),
        test_atom(4, "CA", "A", 4, vec3(5.0, 0.0, 0.0)),
    ];
    let molecule = Molecule {
        atoms,
        helices: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 1,
            start_insertion_code: String::new(),
            end: 1,
            end_insertion_code: String::new(),
        }],
        sheets: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 3,
            start_insertion_code: String::new(),
            end: 4,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let structure = molecule.atomic_structure();
    assert_eq!(structure.ranges.polymer_ranges, vec![0, 0, 2, 3]);
    assert_eq!(
        structure
            .model
            .hierarchy
            .derived
            .residue
            .trace_element_index,
        vec![Some(0), None, Some(2), Some(3)]
    );

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""secondary_type":"helix","chain":"A","residue_start":1,"residue_end":1,"group_id":0,"polymer_trace":{"initial":true,"final":true,"sec_struc_first":false,"sec_struc_last":false}"#));
    assert!(summary.contains(r#""secondary_type":"sheet","chain":"A","residue_start":3,"residue_end":3,"group_id":1,"polymer_trace":{"initial":true,"final":false,"sec_struc_first":false,"sec_struc_last":false}"#));
    assert!(summary.contains(r#""secondary_type":"sheet","chain":"A","residue_start":4,"residue_end":4,"group_id":2,"polymer_trace":{"initial":false,"final":true,"sec_struc_first":false,"sec_struc_last":false}"#));
}

#[test]
fn molstar_trace_initial_final_flags_follow_polymer_ranges_not_chain_boundaries() {
    fn peptide_atom(id: &mut usize, name: &str, seq: i32, x: f32) -> Atom {
        let atom = test_atom(*id, name, "A", seq, vec3(x, 0.0, 0.0));
        *id += 1;
        atom
    }

    let mut id = 1;
    let atoms = vec![
        peptide_atom(&mut id, "N", 1, 0.0),
        peptide_atom(&mut id, "CA", 1, 0.5),
        peptide_atom(&mut id, "C", 1, 1.0),
        peptide_atom(&mut id, "N", 2, 1.2),
        peptide_atom(&mut id, "CA", 2, 1.5),
        peptide_atom(&mut id, "C", 2, 2.0),
        peptide_atom(&mut id, "N", 3, 20.0),
        peptide_atom(&mut id, "CA", 3, 3.0),
        peptide_atom(&mut id, "C", 3, 3.5),
        peptide_atom(&mut id, "N", 4, 3.7),
        peptide_atom(&mut id, "CA", 4, 4.0),
        peptide_atom(&mut id, "C", 4, 4.5),
    ];
    let molecule = Molecule {
        atoms,
        helices: vec![
            SecondaryRange {
                chain: "A".to_string(),
                start: 1,
                start_insertion_code: String::new(),
                end: 2,
                end_insertion_code: String::new(),
            },
            SecondaryRange {
                chain: "A".to_string(),
                start: 3,
                start_insertion_code: String::new(),
                end: 4,
                end_insertion_code: String::new(),
            },
        ],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let structure = molecule.atomic_structure();
    assert_eq!(structure.ranges.polymer_ranges.len(), 4);

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""secondary_type":"helix","chain":"A","residue_start":1,"residue_end":1,"group_id":0,"polymer_trace":{"initial":true,"final":false,"sec_struc_first":false,"sec_struc_last":false}"#));
    assert!(summary.contains(r#""secondary_type":"helix","chain":"A","residue_start":2,"residue_end":2,"group_id":1,"polymer_trace":{"initial":false,"final":true,"sec_struc_first":false,"sec_struc_last":false}"#));
    assert!(summary.contains(r#""secondary_type":"helix","chain":"A","residue_start":3,"residue_end":3,"group_id":2,"polymer_trace":{"initial":true,"final":false,"sec_struc_first":false,"sec_struc_last":false}"#));
    assert!(summary.contains(r#""secondary_type":"helix","chain":"A","residue_start":4,"residue_end":4,"group_id":3,"polymer_trace":{"initial":false,"final":true,"sec_struc_first":false,"sec_struc_last":false}"#));
}

#[test]
fn secondary_trace_respects_insertion_code_boundaries() {
    fn inserted_atom(id: usize, seq: i32, insertion_code: &str, x: f32) -> Atom {
        let mut atom = test_atom(id, "CA", "A", seq, vec3(x, 0.0, 0.0));
        atom.insertion_code = insertion_code.to_string();
        atom
    }

    let molecule = Molecule {
        atoms: vec![
            inserted_atom(1, 1, "", 0.0),
            inserted_atom(2, 1, "A", 1.0),
            inserted_atom(3, 2, "", 2.0),
            inserted_atom(4, 2, "A", 3.0),
        ],
        helices: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 1,
            start_insertion_code: "A".to_string(),
            end: 2,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(
        r#""secondary_type":"helix","chain":"A","residue_start":1,"residue_end":1,"group_id":1"#
    ));
    assert!(summary.contains(
        r#""secondary_type":"helix","chain":"A","residue_start":2,"residue_end":2,"group_id":2"#
    ));
    assert!(summary.contains(
        r#""secondary_type":"coil","chain":"A","residue_start":1,"residue_end":1,"group_id":0"#
    ));
    assert!(!summary
        .contains(r#""secondary_type":"coil","chain":"A","residue_start":2,"residue_end":2"#));

    let trace_centers = build_render_objects(&molecule, &options)
        .into_iter()
        .filter_map(|object| {
            if let RenderObject::PolymerTraceSegment { controls, .. } = object {
                Some(controls.p2.to_vec3())
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    assert_eq!(trace_centers.len(), 2);
    assert_eq!(trace_centers[0], vec3(1.0, 0.0, 0.0));
    assert_eq!(trace_centers[1], vec3(2.0, 0.0, 0.0));
}

#[test]
fn model_secondary_structure_defaults_to_none_without_categories() {
    let molecule = Molecule {
        atoms: (1..=3)
            .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
            .collect(),
        ..Molecule::default()
    };

    let structure = molecule.atomic_structure();
    let secondary = &structure.model.secondary_structure;
    assert_eq!(
        secondary.residue_type,
        vec![SecondaryStructureType::NONE; 3]
    );
    assert_eq!(secondary.key, vec![0, 0, 0]);
    assert_eq!(secondary.elements, vec![SecondaryStructureElement::None]);
}

#[test]
fn model_secondary_structure_assigns_helix_and_sheet_types_by_residue() {
    let molecule = Molecule {
        atoms: (1..=5)
            .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
            .collect(),
        helices: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 2,
            start_insertion_code: String::new(),
            end: 3,
            end_insertion_code: String::new(),
        }],
        sheets: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 4,
            start_insertion_code: String::new(),
            end: 5,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };

    let structure = molecule.atomic_structure();
    let secondary = &structure.model.secondary_structure;
    assert_eq!(
        secondary.residue_type,
        vec![
            SecondaryStructureType::NONE,
            SecondaryStructureType::HELIX,
            SecondaryStructureType::HELIX,
            SecondaryStructureType::BETA_SHEET,
            SecondaryStructureType::BETA_SHEET,
        ]
    );
    assert_eq!(secondary.key, vec![0, 1, 1, 2, 2]);
    assert_eq!(
        secondary.elements,
        vec![
            SecondaryStructureElement::None,
            SecondaryStructureElement::Helix,
            SecondaryStructureElement::Sheet,
        ]
    );
}

#[test]
fn secondary_structure_assignment_summary_matches_molstar_reference() {
    let Some(molstar_secondary_structure) =
        read_molstar_source("mol-model-formats/structure/property/secondary-structure.ts")
    else {
        eprintln!("skipping pinned Mol* secondary-structure source audit; artifacts is absent");
        return;
    };
    let Some(molstar_secondary_types) = read_molstar_source("mol-model/structure/model/types.ts")
    else {
        eprintln!("skipping pinned Mol* model types source audit; artifacts is absent");
        return;
    };

    assert!(molstar_secondary_structure
        .contains("elements: SecondaryStructure.Element[] = [{ kind: 'none' }]"));
    assert!(molstar_secondary_structure.contains("data.type[rI] = type;"));
    assert!(molstar_secondary_structure.contains("data.key[rI] = key;"));
    assert!(molstar_secondary_types.contains("Helix = 0x2"));
    assert!(molstar_secondary_types.contains("Beta = 0x4"));
    assert!(molstar_secondary_types.contains("BetaSheet = 0x800000"));

    let molecule = Molecule {
        atoms: (1..=5)
            .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
            .collect(),
        helices: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 2,
            start_insertion_code: String::new(),
            end: 3,
            end_insertion_code: String::new(),
        }],
        sheets: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 4,
            start_insertion_code: String::new(),
            end: 4,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };

    assert_eq!(
        secondary_structure_assignment_summary_json(&molecule),
        include_str!("../../tests/expected/secondary-structure-assignment-summary.json").trim_end()
    );
}

#[test]
fn render_object_summary_matches_molstar_representative_snapshot() {
    let molecule = Molecule {
        atoms: (1..=5)
            .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
            .collect(),
        helices: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 2,
            start_insertion_code: String::new(),
            end: 3,
            end_insertion_code: String::new(),
        }],
        sheets: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 4,
            start_insertion_code: String::new(),
            end: 4,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    assert_eq!(
        render_object_summary_json(&molecule, &options),
        include_str!("../../tests/expected/render-object-summary-molstar-representative.json")
            .trim_end()
    );
}

#[test]
fn render_object_span_summary_maps_mesh_faces_to_molstar_sparse_stl_slots() {
    let molecule = Molecule {
        atoms: (1..=5)
            .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
            .collect(),
        helices: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 2,
            start_insertion_code: String::new(),
            end: 3,
            end_insertion_code: String::new(),
        }],
        sheets: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 4,
            start_insertion_code: String::new(),
            end: 4,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = render_object_span_summary_json(&molecule, &options);

    assert!(summary
        .contains(r#""face_start":0,"face_end":160,"stl_facet_start":0,"stl_facet_end":480"#));
    assert!(summary
        .contains(r#""face_start":160,"face_end":432,"stl_facet_start":480,"stl_facet_end":1296"#));
    assert!(summary.contains(r#""first_face":{"indices":["#));
    assert!(summary.contains(r#""vertices":[["#));
    assert!(summary.contains(r#""normals":[["#));
}

#[test]
fn model_secondary_structure_boundaries_match_molstar_numeric_scan() {
    fn boundary_atom(id: usize, label_seq_id: &str, insertion_code: &str) -> Atom {
        let mut atom = test_atom(id, "CA", "A", id as i32, vec3(id as f32, 0.0, 0.0));
        atom.residue_seq = label_seq_id.to_string();
        atom.auth_residue_seq = label_seq_id.to_string();
        atom.insertion_code = insertion_code.to_string();
        atom
    }

    let molecule = Molecule {
        atoms: vec![
            boundary_atom(1, "01", ""),
            boundary_atom(2, "02", ""),
            boundary_atom(3, "02", "A"),
            boundary_atom(4, "03", ""),
        ],
        helices: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 2,
            start_insertion_code: String::new(),
            end: 2,
            end_insertion_code: "Z".to_string(),
        }],
        ..Molecule::default()
    };

    let structure = molecule.atomic_structure();
    let secondary = &structure.model.secondary_structure;
    assert_eq!(
        secondary.residue_type,
        vec![
            SecondaryStructureType::NONE,
            SecondaryStructureType::HELIX,
            SecondaryStructureType::HELIX,
            SecondaryStructureType::HELIX,
        ]
    );
    assert_eq!(secondary.key, vec![0, 1, 1, 1]);
}

#[test]
fn molstar_single_residue_trace_keeps_initial_final_flags() {
    let molecule = Molecule {
        atoms: vec![test_atom(1, "CA", "A", 1, vec3(0.0, 0.0, 0.0))],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""geometry_type":"sphere""#));
    assert!(summary.contains(
        r#""valueCell":{"drawCount":960,"uVertexCount":162,"uGroupCount":1,"instanceCount":1,"uInstanceCount":1}"#
    ));
    assert!(summary.contains(r#""polymer_trace":{"initial":true,"final":true,"sec_struc_first":false,"sec_struc_last":false}"#));

    let mesh = build_mesh(&molecule, &options);
    assert!(!mesh.faces.is_empty());
    assert_eq!(mesh.faces.len(), mesh.face_groups.len());
}

#[test]
fn molstar_cyclic_trace_suppresses_terminal_initial_final_flags() {
    let mut atoms = Vec::new();
    for (seq, n_x, ca_x, c_x, o_x) in [
        (1, 0.0, 0.5, 1.0, 1.1),
        (2, 2.0, 2.5, 3.0, 3.1),
        (3, 3.2, 3.5, 0.4, 0.5),
    ] {
        atoms.push(test_atom(
            atoms.len() + 1,
            "N",
            "A",
            seq,
            vec3(n_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "CA",
            "A",
            seq,
            vec3(ca_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "C",
            "A",
            seq,
            vec3(c_x, 0.0, 0.0),
        ));
        atoms.push(test_atom(
            atoms.len() + 1,
            "O",
            "A",
            seq,
            vec3(o_x, 0.0, 0.0),
        ));
    }
    for atom in &mut atoms {
        atom.entity_id = "1".to_string();
    }
    let entities = vec![Entity {
        id: "1".to_string(),
        type_name: "polymer".to_string(),
        description: String::new(),
    }];
    let entity_poly_seq = (1..=3)
        .map(|num| EntityPolySeq {
            entity_id: "1".to_string(),
            num,
            mon_id: "ALA".to_string(),
            hetero: "n".to_string(),
        })
        .collect::<Vec<_>>();
    let molecule = Molecule {
        atoms,
        entity_index: EntityIndexMap::from_entities(&entities, &[], &[]),
        entities,
        entity_poly_seq,
        ..Molecule::default()
    };
    let structure = molecule.atomic_structure();
    assert_eq!(structure.ranges.cyclic_polymer_map.get(&0), Some(&2));
    assert_eq!(structure.ranges.cyclic_polymer_map.get(&2), Some(&0));

    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };
    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""secondary_type":"coil","chain":"A","residue_start":1,"residue_end":1,"group_id":0,"polymer_trace":{"initial":false,"final":false,"sec_struc_first":false,"sec_struc_last":false}"#));
    assert!(summary.contains(r#""secondary_type":"coil","chain":"A","residue_start":3,"residue_end":3,"group_id":2,"polymer_trace":{"initial":false,"final":false,"sec_struc_first":false,"sec_struc_last":false}"#));
    assert!(!summary.contains(r#""polymer_trace":{"initial":true,"final":true"#));
}

#[test]
fn molstar_coil_bridges_missing_trace_residues_as_polymer_gap() {
    let molecule = Molecule {
        atoms: [1, 2, 3, 6, 7, 8]
            .into_iter()
            .map(|seq| {
                test_atom(
                    seq as usize,
                    "CA",
                    "A",
                    seq,
                    vec3(seq as f32 * 2.0, 0.0, 0.0),
                )
            })
            .collect(),
        helices: vec![
            SecondaryRange {
                chain: "A".to_string(),
                start: 1,
                start_insertion_code: String::new(),
                end: 2,
                end_insertion_code: String::new(),
            },
            SecondaryRange {
                chain: "A".to_string(),
                start: 7,
                start_insertion_code: String::new(),
                end: 8,
                end_insertion_code: String::new(),
            },
        ],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(
            r#""geometry_type":"tube","visual":"polymer-trace","representation":"molstar","secondary_type":"coil","chain":"A","residue_start":3,"residue_end":3"#
        ));
    assert!(summary.contains(
            r#""geometry_type":"tube","visual":"polymer-trace","representation":"molstar","secondary_type":"coil","chain":"A","residue_start":6,"residue_end":6"#
        ));
}

#[test]
fn polymer_trace_iterator_reference_json_matches_molstar_snapshot() {
    let Some(molstar_trace_iterator) =
        read_molstar_source("mol-repr/structure/visual/util/polymer/trace-iterator.ts")
    else {
        eprintln!("skipping pinned Mol* polymer trace iterator source audit; artifacts is absent");
        return;
    };
    let Some(molstar_helix_orientation) =
        read_molstar_source("mol-model-props/computed/helix-orientation/helix-orientation.ts")
    else {
        eprintln!("skipping pinned Mol* helix orientation source audit; artifacts is absent");
        return;
    };
    for field in [
        "centerPrev: StructureElement.Location",
        "centerNext: StructureElement.Location",
        "first: boolean, last: boolean",
        "initial: boolean, final: boolean",
        "secStrucFirst: boolean, secStrucLast: boolean",
        "coarseBackboneFirst: boolean, coarseBackboneLast: boolean",
        "p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, p4: Vec3",
        "d12: Vec3, d23: Vec3",
        "private p6 = Vec3()",
        "private d34 = Vec3()",
    ] {
        assert!(molstar_trace_iterator.contains(field));
    }
    for snippet in [
        "private getResidueIndex(residueIndex: number) {",
        "const cyclicIndex = this.cyclicPolymerMap.get(this.residueSegmentMin);",
        "residueIndex = cyclicIndex - (this.residueSegmentMin - residueIndex - 1);",
        "const cyclicIndex = this.cyclicPolymerMap.get(this.residueSegmentMax);",
        "residueIndex = cyclicIndex + (residueIndex - this.residueSegmentMax - 1);",
        "value.initial = residueIndex === residueIndexPrev1;",
        "value.final = residueIndex === residueIndexNext1;",
        "this.setControlPoint(value.p0, this.p0, this.p1, this.p2, ssPrev2);",
        "this.setControlPoint(value.p4, this.p4, this.p5, this.p6, ssNext2);",
        "if (this.helixOrientationCenters && isHelixSS(ss))",
        "Vec3.set(out, 1, 0, 0);",
        "handle positions for tubular helices",
        "const helixFlag = isHelix && this.helixOrientationCenters;",
        "const enum CoarsePolymerTraceIteratorState { nextPolymer, nextElement }",
        "export class CoarsePolymerTraceIterator implements Iterator<PolymerTraceElement>",
        "private getElementIndex(elementIndex: number) {",
        "Math.min(Math.max(this.polymerSegment.start, elementIndex), this.polymerSegment.end - 1)",
        "const f = 0.5;",
        "Vec3.set(this.value.d12, 1, 0, 0);",
        "this.hasNext = this.elementIndex + 1 < this.polymerSegment.end || this.polymerIt.hasNext;",
    ] {
        assert!(
            molstar_trace_iterator.contains(snippet),
            "missing Mol* trace iterator snippet: {snippet}"
        );
    }
    for snippet in [
        "export function calcHelixOrientation(model: Model): HelixOrientation",
        "const centers = new Float32Array(n * 3);",
        "j = (index - 2);",
        "Vec3.projectPointOnVector(vt, vt, axis, v1);",
    ] {
        assert!(
            molstar_helix_orientation.contains(snippet),
            "missing Mol* helix orientation snippet: {snippet}"
        );
    }

    let missing_bridge = Molecule {
        atoms: [1, 2, 3, 6, 7, 8]
            .into_iter()
            .map(|seq| {
                test_atom(
                    seq as usize,
                    "CA",
                    "A",
                    seq,
                    vec3(seq as f32 * 2.0, 0.0, 0.0),
                )
            })
            .collect(),
        helices: vec![
            SecondaryRange {
                chain: "A".to_string(),
                start: 1,
                start_insertion_code: String::new(),
                end: 2,
                end_insertion_code: String::new(),
            },
            SecondaryRange {
                chain: "A".to_string(),
                start: 7,
                start_insertion_code: String::new(),
                end: 8,
                end_insertion_code: String::new(),
            },
        ],
        ..Molecule::default()
    };

    let mut cyclic_atoms = Vec::new();
    for (seq, n_x, ca_x, c_x, o_x) in [
        (1, 0.0, 0.5, 1.0, 1.1),
        (2, 2.0, 2.5, 3.0, 3.1),
        (3, 3.2, 3.5, 0.4, 0.5),
    ] {
        cyclic_atoms.push(test_atom(
            cyclic_atoms.len() + 1,
            "N",
            "A",
            seq,
            vec3(n_x, 0.0, 0.0),
        ));
        cyclic_atoms.push(test_atom(
            cyclic_atoms.len() + 1,
            "CA",
            "A",
            seq,
            vec3(ca_x, 0.0, 0.0),
        ));
        cyclic_atoms.push(test_atom(
            cyclic_atoms.len() + 1,
            "C",
            "A",
            seq,
            vec3(c_x, 0.0, 0.0),
        ));
        cyclic_atoms.push(test_atom(
            cyclic_atoms.len() + 1,
            "O",
            "A",
            seq,
            vec3(o_x, 0.0, 0.0),
        ));
    }
    for atom in &mut cyclic_atoms {
        atom.entity_id = "1".to_string();
    }
    let entities = vec![Entity {
        id: "1".to_string(),
        type_name: "polymer".to_string(),
        description: String::new(),
    }];
    let entity_poly_seq = (1..=3)
        .map(|num| EntityPolySeq {
            entity_id: "1".to_string(),
            num,
            mon_id: "ALA".to_string(),
            hetero: "n".to_string(),
        })
        .collect::<Vec<_>>();
    let cyclic = Molecule {
        atoms: cyclic_atoms,
        entity_index: EntityIndexMap::from_entities(&entities, &[], &[]),
        entities,
        entity_poly_seq,
        ..Molecule::default()
    };

    let mut oriented_atoms = Vec::new();
    for (seq, ca) in [
        (1, vec3(1.0, 0.0, 0.0)),
        (2, vec3(0.0, 1.0, 1.5)),
        (3, vec3(-1.0, 0.0, 3.0)),
        (4, vec3(0.0, -1.0, 4.5)),
        (5, vec3(1.0, 0.0, 6.0)),
    ] {
        oriented_atoms.push(test_atom(
            oriented_atoms.len() + 1,
            "N",
            "H",
            seq,
            ca - vec3(0.25, 0.0, 0.0),
        ));
        oriented_atoms.push(test_atom(oriented_atoms.len() + 1, "CA", "H", seq, ca));
        oriented_atoms.push(test_atom(
            oriented_atoms.len() + 1,
            "C",
            "H",
            seq,
            ca + vec3(0.25, 0.0, 0.0),
        ));
        oriented_atoms.push(test_atom(
            oriented_atoms.len() + 1,
            "O",
            "H",
            seq,
            ca + vec3(0.25, 0.6, 0.0),
        ));
    }
    let helix_orientation = Molecule {
        atoms: oriented_atoms,
        helices: vec![SecondaryRange {
            chain: "H".to_string(),
            start: 1,
            start_insertion_code: String::new(),
            end: 5,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };

    let coarse_entities = vec![
        Entity {
            id: "1".to_string(),
            type_name: "polymer".to_string(),
            description: String::new(),
        },
        Entity {
            id: "2".to_string(),
            type_name: "polymer".to_string(),
            description: String::new(),
        },
    ];
    let coarse_trace = Molecule {
        entity_index: EntityIndexMap::from_entities(&coarse_entities, &[], &[]),
        entities: coarse_entities,
        coarse_spheres: vec![
            CoarseSphere {
                id: 1,
                model_num: 1,
                entity_id: "1".to_string(),
                asym_id: "S".to_string(),
                seq_id_begin: 1,
                seq_id_end: 10,
                position: vec3(0.0, 0.0, 0.0),
                radius: 1.0,
                rmsf: 0.0,
            },
            CoarseSphere {
                id: 2,
                model_num: 1,
                entity_id: "1".to_string(),
                asym_id: "S".to_string(),
                seq_id_begin: 11,
                seq_id_end: 20,
                position: vec3(2.0, 0.0, 0.0),
                radius: 1.0,
                rmsf: 0.0,
            },
            CoarseSphere {
                id: 3,
                model_num: 1,
                entity_id: "1".to_string(),
                asym_id: "S".to_string(),
                seq_id_begin: 25,
                seq_id_end: 30,
                position: vec3(6.0, 0.0, 0.0),
                radius: 1.0,
                rmsf: 0.0,
            },
        ],
        coarse_gaussians: vec![
            CoarseGaussian {
                id: 1,
                model_num: 1,
                entity_id: "2".to_string(),
                asym_id: "G".to_string(),
                seq_id_begin: 1,
                seq_id_end: 4,
                position: vec3(0.0, 5.0, 0.0),
                weight: 1.0,
                covariance: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            },
            CoarseGaussian {
                id: 2,
                model_num: 1,
                entity_id: "2".to_string(),
                asym_id: "G".to_string(),
                seq_id_begin: 5,
                seq_id_end: 8,
                position: vec3(0.0, 8.0, 0.0),
                weight: 1.0,
                covariance: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            },
        ],
        ..Molecule::default()
    };

    let actual = format!(
        "[{},{},{},{},{}]",
        polymer_trace_iterator_reference_json("missing-residue-bridge", &missing_bridge),
        polymer_trace_iterator_reference_json("cyclic-tripeptide", &cyclic),
        polymer_trace_iterator_reference_json_with_helix_orientation(
            "helix-orientation-centers",
            &helix_orientation
        ),
        coarse_polymer_trace_iterator_reference_json(
            "coarse-sphere-trace",
            &coarse_trace,
            CoarseElementKind::Spheres
        ),
        coarse_polymer_trace_iterator_reference_json(
            "coarse-gaussian-trace",
            &coarse_trace,
            CoarseElementKind::Gaussians
        )
    );
    for snippet in [
        r#""case":"missing-residue-bridge""#,
        r#""polymer_ranges":[[0,2],[3,5]]"#,
        r#""center_prev":1,"center":2,"center_next":2,"first":false,"last":true,"initial":false,"final":true"#,
        r#""p0":[2,0,0],"p1":[4,0,0],"p2":[6,0,0],"p3":[7.5,0,0],"p4":[9,0,0]"#,
        r#""center_prev":3,"center":3,"center_next":4,"first":true,"last":false,"initial":true,"final":false"#,
        r#""case":"cyclic-tripeptide""#,
        r#""cyclic_polymer_map":[[0,2],[2,0]]"#,
        r#""center_prev":9,"center":1,"center_next":5,"first":true,"last":false,"initial":false,"final":false"#,
        r#""p0":[2.5,0,0],"p1":[3.5,0,0],"p2":[0.5,0,0],"p3":[2.5,0,0],"p4":[3.5,0,0]"#,
        r#""center_prev":5,"center":9,"center_next":1,"first":false,"last":true,"initial":false,"final":false"#,
        r#""case":"helix-orientation-centers","use_helix_orientation":true"#,
        r#""d12":[1,0,0],"d23":[1,0,0]"#,
        r#""case":"coarse-sphere-trace","unit_kind":"spheres""#,
        r#""polymer_ranges":[[0,1],[2,2]]"#,
        r#""unit_id":0,"element_index":0,"source_element":0,"chain":"S","seq_begin":1,"seq_end":10,"center_prev":0,"center":0,"center_next":1,"first":true,"last":false,"initial":false,"final":false"#,
        r#""p0":[-2,0,0],"p1":[-1,0,0],"p2":[0,0,0],"p3":[2,0,0],"p4":[3,0,0]"#,
        r#""unit_id":0,"element_index":2,"source_element":2,"chain":"S","seq_begin":25,"seq_end":30,"center_prev":2,"center":2,"center_next":2,"first":true,"last":true,"initial":false,"final":false"#,
        r#""case":"coarse-gaussian-trace","unit_kind":"gaussians""#,
        r#""unit_id":1,"element_index":0,"source_element":0,"chain":"G","seq_begin":1,"seq_end":4,"center_prev":0,"center":0,"center_next":1,"first":true,"last":false,"initial":false,"final":false"#,
        r#""p0":[0,2,0],"p1":[0,3.5,0],"p2":[0,5,0],"p3":[0,8,0],"p4":[0,9.5,0]"#,
    ] {
        assert!(
            actual.contains(snippet),
            "missing trace iterator reference snippet: {snippet}"
        );
    }
    assert_eq!(
        actual,
        include_str!("../../tests/expected/polymer-trace-iterator-reference.json").trim_end()
    );
}

#[test]
fn cartoon_representation_summary_selects_trace_gap_and_nucleotide_ring_visuals() {
    let mut atoms = [1, 2, 5]
        .into_iter()
        .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
        .collect::<Vec<_>>();
    let nucleotide_start = atoms.len();
    for (offset, name) in ["O3'", "C1'", "N1", "C2", "N3", "C4", "C5", "C6"]
        .into_iter()
        .enumerate()
    {
        let mut atom = test_atom(
            nucleotide_start + offset + 1,
            name,
            "N",
            1,
            vec3(offset as f32 * 0.2, 2.0, 0.0),
        );
        atom.residue = "C".to_string();
        atoms.push(atom);
    }
    let molecule = Molecule {
        atoms,
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = representation_summary_json(&molecule, &options);
    assert!(summary.contains(r#""name":"cartoon""#));
    assert!(
        summary.contains(r#""selected_visuals":["polymer-trace","nucleotide-ring","polymer-gap"]"#)
    );
    assert!(
        summary.contains(r#""realized_visuals":["polymer-trace","nucleotide-ring","polymer-gap"]"#)
    );
}

#[test]
fn polymer_gap_visual_emits_molstar_fixed_count_dashed_cylinders() {
    let molecule = Molecule {
        atoms: [1, 2, 5]
            .into_iter()
            .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
            .collect(),
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        visuals: vec!["polymer-gap".to_string()],
        radial_segments: 3,
        ..MeshOptions::default()
    };

    let gaps = build_render_objects(&molecule, &options)
        .into_iter()
        .filter_map(|object| {
            if let RenderObject::FixedCountDashedCylinder {
                start,
                end,
                radius,
                length_scale,
                segment_count,
            } = object
            {
                Some((start, end, radius, length_scale, segment_count))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    assert_eq!(gaps.len(), 2);

    for (gap, expected_start, expected_end) in [
        (&gaps[0], vec3(2.0, 0.0, 0.0), vec3(5.0, 0.0, 0.0)),
        (&gaps[1], vec3(5.0, 0.0, 0.0), vec3(2.0, 0.0, 0.0)),
    ] {
        assert_vec3_close(gap.0, expected_start, 0.000_001);
        assert_vec3_close(gap.1, expected_end, 0.000_001);
        assert!((gap.2 - 0.20).abs() < 0.000_001);
        assert_eq!(gap.3, 0.5);
        assert_eq!(gap.4, 10);
    }

    let summary = render_object_summary_json(&molecule, &options);
    assert_eq!(summary.matches(r#""visual":"polymer-gap""#).count(), 2);
    assert!(summary.contains(r#""geometry_type":"dashed-cylinder","visual":"polymer-gap","representation":"cartoon","secondary_type":"gap","chain":"A","residue_start":2,"residue_end":5,"group_id":0"#));
    assert!(summary.contains(r#""geometry_type":"dashed-cylinder","visual":"polymer-gap","representation":"cartoon","secondary_type":"gap","chain":"A","residue_start":5,"residue_end":2,"group_id":1"#));

    let mesh = build_mesh(
        &molecule,
        &MeshOptions {
            representation: Representation::Backbone,
            center: false,
            assembly: None,
            visuals: vec!["polymer-gap".to_string()],
            radial_segments: 3,
            ..MeshOptions::default()
        },
    );
    assert_eq!(mesh.vertices.len(), 180);
    assert_eq!(mesh.normals.len(), 180);
    assert_eq!(mesh.faces.len(), 80);
    assert_eq!(mesh.face_groups.len(), 80);
    assert_eq!(mesh.group_count, 2);
    assert!(faces_have_valid_indices(&mesh));
}

#[test]
fn polymer_gap_radius_ignores_type_symbol_for_uniform_size_theme() {
    let mut atoms = [1, 2, 5]
        .into_iter()
        .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
        .collect::<Vec<_>>();
    atoms[1].element = "C".to_string();
    atoms[1].type_symbol = "H".to_string();
    let molecule = Molecule {
        atoms,
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        visuals: vec!["polymer-gap".to_string()],
        ..MeshOptions::default()
    };

    let gaps = build_render_objects(&molecule, &options)
        .into_iter()
        .filter_map(|object| {
            if let RenderObject::FixedCountDashedCylinder { radius, .. } = object {
                Some(radius)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    assert_eq!(gaps.len(), 2);
    assert!((gaps[0] - 0.20).abs() <= 0.000_001);
    assert!((gaps[1] - 0.20).abs() <= 0.000_001);
}

#[test]
fn molstar_representation_summary_selects_branched_ball_and_stick_for_carbohydrates() {
    let molecule = Molecule {
        atoms: vec![
            carbohydrate_atom(1, "C1", "A", 1, "MAN", vec3(0.0, 0.0, 0.0)),
            carbohydrate_atom(2, "C2", "A", 1, "MAN", vec3(1.0, 0.0, 0.0)),
            carbohydrate_atom(3, "C3", "A", 1, "MAN", vec3(1.5, 1.0, 0.0)),
            carbohydrate_atom(4, "C4", "A", 1, "MAN", vec3(1.0, 2.0, 0.0)),
            carbohydrate_atom(5, "C5", "A", 1, "MAN", vec3(0.0, 2.0, 0.0)),
            carbohydrate_atom(6, "O5", "A", 1, "MAN", vec3(-0.5, 1.0, 0.0)),
        ],
        bonds: carbohydrate_bonds(&[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]),
        bond_metadata: carbohydrate_bond_metadata(6),
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = representation_summary_json(&molecule, &options);
    assert!(summary.contains(r#""name":"molstar""#));
    assert!(summary.contains(r#""selected_visuals":["element-sphere","intra-bond","inter-bond"]"#));
    assert!(summary.contains(r#""realized_visuals":["element-sphere","intra-bond"]"#));

    let render_summary = render_object_summary_json(&molecule, &options);
    assert!(render_summary.contains(
        r#""visual":"element-sphere","representation":"molstar","secondary_type":"branched""#
    ));
    assert!(render_summary.contains(
        r#""visual":"intra-bond","representation":"molstar","secondary_type":"branched""#
    ));
    assert!(!render_summary.contains(r#""visual":"carbohydrate-symbol""#));

    let cartoon_options = MeshOptions {
        representation: Representation::Cartoon,
        ..options
    };
    let cartoon_summary = representation_summary_json(&molecule, &cartoon_options);
    assert!(cartoon_summary.contains(r#""selected_visuals":[]"#));
}

#[test]
fn molstar_representation_summary_realizes_ligand_branched_and_ion_components() {
    let atoms = vec![
        test_atom(1, "CA", "A", 1, vec3(0.0, 0.0, 0.0)),
        test_atom(2, "CA", "A", 2, vec3(1.0, 0.0, 0.0)),
        het_atom(3, "C1", "L", 1, "LIG", vec3(3.0, 0.0, 0.0)),
        het_atom(4, "O1", "L", 1, "LIG", vec3(4.0, 0.0, 0.0)),
        carbohydrate_atom(5, "C1", "B", 1, "MAN", vec3(0.0, 3.0, 0.0)),
        carbohydrate_atom(6, "C2", "B", 1, "MAN", vec3(1.0, 3.0, 0.0)),
        het_atom(7, "ZN", "I", 1, "ZN", vec3(6.0, 0.0, 0.0)),
    ];
    let molecule = Molecule {
        atoms,
        bonds: vec![Bond { a: 2, b: 3 }, Bond { a: 4, b: 5 }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let representation = representation_summary_json(&molecule, &options);
    assert!(representation.contains(
        r#""selected_visuals":["polymer-trace","element-sphere","intra-bond","inter-bond"]"#
    ));
    assert!(representation
        .contains(r#""realized_visuals":["polymer-trace","element-sphere","intra-bond"]"#));

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""visual":"element-sphere","representation":"molstar","secondary_type":"ligand","chain":"L""#));
    assert!(summary.contains(
        r#""visual":"intra-bond","representation":"molstar","secondary_type":"ligand","chain":"L""#
    ));
    assert!(summary.contains(r#""visual":"element-sphere","representation":"molstar","secondary_type":"branched","chain":"B""#));
    assert!(summary.contains(r#""visual":"intra-bond","representation":"molstar","secondary_type":"branched","chain":"B""#));
    assert!(summary.contains(
        r#""visual":"element-sphere","representation":"molstar","secondary_type":"ion","chain":"I""#
    ));
}

#[test]
fn molstar_component_intra_bond_objects_follow_unit_adjacency_slot_order() {
    let molecule = Molecule {
        atoms: vec![
            het_atom(1, "C1", "L", 1, "LIG", vec3(0.0, 0.0, 0.0)),
            het_atom(2, "C2", "L", 2, "LIG", vec3(1.0, 0.0, 0.0)),
            het_atom(3, "C3", "L", 3, "LIG", vec3(0.0, 1.0, 0.0)),
        ],
        bonds: vec![Bond { a: 2, b: 0 }, Bond { a: 0, b: 1 }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = render_object_summary_json(&molecule, &options);
    assert_eq!(summary.matches(r#""visual":"intra-bond""#).count(), 4);

    let group0 = summary
        .find(r#""geometry_type":"cylinder","visual":"intra-bond","representation":"molstar","secondary_type":"ligand","chain":"L","residue_start":1,"residue_end":3,"group_id":0"#)
        .expect("missing Mol* adjacency slot 0 half-link 1->3");
    let group1 = summary
        .find(r#""geometry_type":"cylinder","visual":"intra-bond","representation":"molstar","secondary_type":"ligand","chain":"L","residue_start":1,"residue_end":2,"group_id":1"#)
        .expect("missing Mol* adjacency slot 1 half-link 1->2");
    let group2 = summary
        .find(r#""geometry_type":"cylinder","visual":"intra-bond","representation":"molstar","secondary_type":"ligand","chain":"L","residue_start":2,"residue_end":1,"group_id":2"#)
        .expect("missing Mol* adjacency slot 2 half-link 2->1");
    let group3 = summary
        .find(r#""geometry_type":"cylinder","visual":"intra-bond","representation":"molstar","secondary_type":"ligand","chain":"L","residue_start":3,"residue_end":1,"group_id":3"#)
        .expect("missing Mol* adjacency slot 3 half-link 3->1");

    assert!(group0 < group1);
    assert!(group1 < group2);
    assert!(group2 < group3);
}

#[test]
fn molstar_component_intra_bond_slots_follow_branching_edgebuilder_offsets() {
    let molecule = Molecule {
        atoms: vec![
            het_atom(1, "OH2", "A", 1, "1PE", vec3(0.0, 0.0, 0.0)),
            het_atom(2, "C12", "A", 2, "1PE", vec3(1.0, 0.0, 0.0)),
            het_atom(3, "C22", "A", 3, "1PE", vec3(2.0, 0.0, 0.0)),
            het_atom(4, "OH3", "A", 4, "1PE", vec3(3.0, 0.0, 0.0)),
            het_atom(5, "C13", "A", 5, "1PE", vec3(4.0, 0.0, 0.0)),
            het_atom(6, "C23", "A", 6, "1PE", vec3(5.0, 0.0, 0.0)),
            het_atom(7, "OH4", "A", 7, "1PE", vec3(6.0, 0.0, 0.0)),
            het_atom(8, "C14", "A", 8, "1PE", vec3(7.0, 0.0, 0.0)),
            het_atom(9, "C24", "A", 9, "1PE", vec3(8.0, 0.0, 0.0)),
            het_atom(10, "OH5", "A", 10, "1PE", vec3(9.0, 0.0, 0.0)),
            het_atom(11, "C15", "A", 11, "1PE", vec3(10.0, 0.0, 0.0)),
            het_atom(12, "C25", "A", 12, "1PE", vec3(11.0, 0.0, 0.0)),
            het_atom(13, "OH6", "A", 13, "1PE", vec3(12.0, 0.0, 0.0)),
        ],
        bonds: vec![
            Bond { a: 0, b: 1 },
            Bond { a: 1, b: 2 },
            Bond { a: 2, b: 3 },
            Bond { a: 3, b: 5 },
            Bond { a: 4, b: 5 },
            Bond { a: 4, b: 6 },
            Bond { a: 6, b: 8 },
            Bond { a: 7, b: 8 },
            Bond { a: 7, b: 9 },
            Bond { a: 9, b: 11 },
            Bond { a: 10, b: 11 },
            Bond { a: 10, b: 12 },
        ],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = render_object_summary_json(&molecule, &options);
    assert_eq!(summary.matches(r#""visual":"intra-bond""#).count(), 24);
    assert!(summary.contains(r#""residue_start":4,"residue_end":6,"group_id":6"#));
    assert!(summary.contains(r#""residue_start":5,"residue_end":6,"group_id":7"#));
    assert!(summary.contains(r#""residue_start":5,"residue_end":7,"group_id":8"#));
    assert!(summary.contains(r#""residue_start":6,"residue_end":4,"group_id":9"#));
}

#[test]
fn molstar_ion_only_component_realizes_ball_and_stick_sphere_without_bonds() {
    let molecule = Molecule {
        atoms: vec![het_atom(1, "ZN", "I", 1, "ZN", vec3(6.0, 0.0, 0.0))],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let representation = representation_summary_json(&molecule, &options);
    assert!(representation
        .contains(r#""selected_visuals":["element-sphere","intra-bond","inter-bond"]"#));
    assert!(representation.contains(r#""realized_visuals":["element-sphere"]"#));

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(
        r#""visual":"element-sphere","representation":"molstar","secondary_type":"ion","chain":"I""#
    ));
    assert!(!summary.contains(r#""secondary_type":"ligand""#));
    assert!(!summary.contains(r#""secondary_type":"branched""#));
    assert!(!summary.contains(r#""visual":"intra-bond""#));

    let mesh = build_mesh(&molecule, &options);
    assert_eq!(mesh.vertices.len(), 162);
    assert_eq!(mesh.faces.len(), 320);
    validate_mesh_for_export(&mesh).unwrap();
}

#[test]
fn molstar_water_only_component_realizes_ball_and_stick_spheres_and_bonds() {
    let molecule = Molecule {
        atoms: vec![
            het_atom(1, "O", "W", 1, "HOH", vec3(0.0, 0.0, 0.0)),
            het_atom(2, "H1", "W", 1, "HOH", vec3(0.957, 0.0, 0.0)),
            het_atom(3, "H2", "W", 1, "HOH", vec3(-0.239, 0.927, 0.0)),
        ],
        bonds: vec![Bond { a: 0, b: 1 }, Bond { a: 0, b: 2 }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        sphere_detail: 1,
        ..MeshOptions::default()
    };

    let representation = representation_summary_json(&molecule, &options);
    assert!(representation
        .contains(r#""selected_visuals":["element-sphere","intra-bond","inter-bond"]"#));
    assert!(representation.contains(r#""realized_visuals":["element-sphere","intra-bond"]"#));

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(
        r#""visual":"element-sphere","representation":"molstar","secondary_type":"water","chain":"W""#
    ));
    assert!(!summary.contains(r#""secondary_type":"ligand""#));
    assert!(!summary.contains(r#""secondary_type":"ion""#));
    assert_eq!(summary.matches(r#""visual":"intra-bond""#).count(), 4);

    let mesh = build_mesh(&molecule, &options);
    assert_eq!(mesh.vertices.len(), 422);
    assert_eq!(mesh.faces.len(), 528);
    validate_mesh_for_export(&mesh).unwrap();
}

#[test]
fn carbohydrate_symbol_mesh_visual_uses_molstar_symbol_defaults() {
    let mut molecule = Molecule {
        atoms: vec![
            carbohydrate_atom(1, "C1", "A", 1, "GLC", vec3(0.0, 0.0, 0.0)),
            carbohydrate_atom(2, "C2", "A", 1, "GLC", vec3(1.0, 0.0, 0.0)),
            carbohydrate_atom(3, "C3", "A", 1, "GLC", vec3(1.5, 1.0, 0.0)),
            carbohydrate_atom(4, "C4", "A", 1, "GLC", vec3(1.0, 2.0, 0.0)),
            carbohydrate_atom(5, "C5", "A", 1, "GLC", vec3(0.0, 2.0, 0.0)),
            carbohydrate_atom(6, "O5", "A", 1, "GLC", vec3(-0.5, 1.0, 0.0)),
            carbohydrate_atom(7, "C1", "A", 2, "NAG", vec3(4.0, 0.0, 0.0)),
            carbohydrate_atom(8, "C2", "A", 2, "NAG", vec3(5.0, 0.0, 0.0)),
            carbohydrate_atom(9, "C3", "A", 2, "NAG", vec3(5.5, 1.0, 0.0)),
            carbohydrate_atom(10, "C4", "A", 2, "NAG", vec3(5.0, 2.0, 0.0)),
            carbohydrate_atom(11, "C5", "A", 2, "NAG", vec3(4.0, 2.0, 0.0)),
            carbohydrate_atom(12, "O5", "A", 2, "NAG", vec3(3.5, 1.0, 0.0)),
        ],
        bonds: carbohydrate_bonds(&[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10),
            (10, 11),
            (11, 6),
        ]),
        bond_metadata: carbohydrate_bond_metadata(12),
        ..Molecule::default()
    };
    molecule.refresh_topology_metadata();
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        visuals: vec!["carbohydrate-symbol".to_string()],
        ..MeshOptions::default()
    };

    let representation = representation_summary_json(&molecule, &options);
    assert!(representation.contains(r#""realized_visuals":["carbohydrate-symbol"]"#));

    let summary = render_object_summary_json(&molecule, &options);
    assert_eq!(
        summary.matches(r#""visual":"carbohydrate-symbol""#).count(),
        2
    );
    assert!(summary.contains(r#""visual":"carbohydrate-symbol","representation":"molstar","secondary_type":"carbohydrate","chain":"A","residue_start":1,"residue_end":1,"group_id":0"#));
    assert!(summary.contains(r#""visual":"carbohydrate-symbol","representation":"molstar","secondary_type":"carbohydrate","chain":"A","residue_start":2,"residue_end":2,"group_id":2"#));
    assert!(summary.contains(r#""drawCount":60,"uVertexCount":12"#));
    assert!(summary.contains(r#""drawCount":36,"uVertexCount":24"#));

    let symbols = build_render_objects(&molecule, &options)
        .into_iter()
        .filter_map(|object| match object {
            RenderObject::CarbohydrateSymbol {
                center,
                shape,
                part,
                ..
            } => Some((center, shape, part)),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(symbols.len(), 2);
    assert_vec3_close(symbols[0].0, vec3(0.5, 1.0, 0.0), 0.000_001);
    assert_vec3_close(symbols[1].0, vec3(4.5, 1.0, 0.0), 0.000_001);
    assert_eq!(symbols[0].2, crate::mesh::CarbohydrateSymbolPart::Whole);
    assert_eq!(symbols[1].2, crate::mesh::CarbohydrateSymbolPart::Whole);
}

#[test]
fn ribbon_representation_summary_selects_putty_tube_and_gap_visuals() {
    let molecule = Molecule {
        atoms: [1, 2, 5]
            .into_iter()
            .map(|seq| test_atom(seq as usize, "CA", "A", seq, vec3(seq as f32, 0.0, 0.0)))
            .collect(),
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Ribbon,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = representation_summary_json(&molecule, &options);
    assert!(summary.contains(r#""name":"ribbon""#));
    assert!(summary.contains(r#""selected_visuals":["polymer-tube","polymer-gap"]"#));
    assert!(summary.contains(r#""realized_visuals":["polymer-tube","polymer-gap"]"#));
    assert!(!summary.contains("polymer-trace"));
    assert!(!summary.contains("nucleotide-ring"));
}

#[test]
fn spacefill_representation_summary_selects_element_sphere_visual() {
    let molecule = Molecule {
        atoms: vec![test_atom(1, "O", "A", 1, vec3(0.0, 0.0, 0.0))],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Spacefill,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = representation_summary_json(&molecule, &options);
    assert!(summary.contains(r#""name":"spacefill""#));
    assert!(summary.contains(r#""selected_visuals":["element-sphere"]"#));
    assert!(summary.contains(r#""realized_visuals":["element-sphere"]"#));
}

#[test]
fn ball_and_stick_representation_summary_selects_inter_bond_for_multiple_symmetry_groups() {
    let molecule = Molecule {
        atoms: vec![
            test_atom(1, "C1", "A", 1, vec3(0.0, 0.0, 0.0)),
            test_atom(2, "O1", "A", 1, vec3(1.0, 0.0, 0.0)),
            test_atom(3, "C1", "B", 1, vec3(3.0, 0.0, 0.0)),
            test_atom(4, "O1", "B", 1, vec3(4.0, 0.0, 0.0)),
        ],
        bonds: vec![Bond { a: 0, b: 1 }, Bond { a: 2, b: 3 }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::BallAndStick,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    assert!(molecule.atomic_structure().symmetry_groups.len() > 1);
    let summary = representation_summary_json(&molecule, &options);
    assert!(summary.contains(r#""name":"ball-and-stick""#));
    assert!(summary.contains(r#""selected_visuals":["element-sphere","intra-bond","inter-bond"]"#));
    assert!(summary.contains(r#""realized_visuals":["element-sphere","intra-bond"]"#));
}

#[test]
fn ball_and_stick_representation_summary_selects_default_inter_bond_visual() {
    let molecule = Molecule {
        atoms: vec![
            test_atom(1, "C1", "A", 1, vec3(0.0, 0.0, 0.0)),
            test_atom(2, "O1", "A", 1, vec3(1.0, 0.0, 0.0)),
        ],
        bonds: vec![Bond { a: 0, b: 1 }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::BallAndStick,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    assert_eq!(molecule.atomic_structure().symmetry_groups.len(), 1);
    let summary = representation_summary_json(&molecule, &options);
    assert!(summary.contains(r#""selected_visuals":["element-sphere","intra-bond","inter-bond"]"#));
    assert!(summary.contains(r#""realized_visuals":["element-sphere","intra-bond"]"#));
}

#[test]
fn molstar_summary_routes_sheet_to_sheet_geometry() {
    let molecule = Molecule {
        atoms: (1..=6)
            .map(|seq| {
                test_atom(
                    seq as usize,
                    "CA",
                    "A",
                    seq,
                    vec3(seq as f32, (seq % 2) as f32 * 0.25, 0.0),
                )
            })
            .collect(),
        helices: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 1,
            start_insertion_code: String::new(),
            end: 3,
            end_insertion_code: String::new(),
        }],
        sheets: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 4,
            start_insertion_code: String::new(),
            end: 6,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""geometry_type":"tube","visual":"polymer-trace","representation":"molstar","secondary_type":"helix","chain":"A","residue_start":1,"residue_end":1,"group_id":0"#));
    assert!(summary.contains(r#""geometry_type":"sheet","visual":"polymer-trace","representation":"molstar","secondary_type":"sheet","chain":"A","residue_start":4,"residue_end":4,"group_id":3"#));
    assert!(summary.contains(r#""u_group_count":6"#));
    assert!(summary.contains(r#""uGroupCount":1"#));
}

#[test]
fn molstar_radial_segments_two_routes_polymer_trace_to_ribbon_geometry() {
    let molecule = Molecule {
        atoms: (1..=6)
            .map(|seq| {
                test_atom(
                    seq as usize,
                    "CA",
                    "A",
                    seq,
                    vec3(seq as f32, (seq % 2) as f32 * 0.25, 0.0),
                )
            })
            .collect(),
        helices: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 1,
            start_insertion_code: String::new(),
            end: 3,
            end_insertion_code: String::new(),
        }],
        sheets: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 4,
            start_insertion_code: String::new(),
            end: 6,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        radial_segments: 2,
        ..MeshOptions::default()
    };

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""geometry_type":"ribbon","visual":"polymer-trace","representation":"molstar","secondary_type":"helix","chain":"A","residue_start":1,"residue_end":1,"group_id":0"#));
    assert!(summary.contains(r#""geometry_type":"ribbon","visual":"polymer-trace","representation":"molstar","secondary_type":"sheet","chain":"A","residue_start":4,"residue_end":4,"group_id":3"#));

    let mesh = build_mesh(&molecule, &options);
    assert!(!mesh.vertices.is_empty());
    assert_eq!(mesh.vertices.len(), mesh.normals.len());
    assert!(faces_have_valid_indices(&mesh));
}

#[test]
fn molstar_sheet_mesh_and_ply_faces_use_valid_finite_indices() {
    let molecule = sheet_test_molecule();
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };
    let mesh = build_mesh(&molecule, &options);

    assert!(!mesh.vertices.is_empty());
    assert!(!mesh.faces.is_empty());
    assert_eq!(mesh.vertices.len(), mesh.normals.len());
    assert_eq!(mesh.faces.len(), mesh.face_groups.len());
    assert!(mesh.group_count >= 1);
    assert!(mesh.vertices.iter().all(|v| vec3_is_finite(*v)));
    assert!(mesh.normals.iter().all(|v| vec3_is_finite(*v)));
    assert!(faces_have_valid_indices(&mesh));

    let pdb = b"SHEET    1 S1 1 SER A   1  THR A   4  0\nATOM      1  CA  SER A   1       0.000   0.000   0.000  1.00 10.00           C\nATOM      2  CA  THR A   2       1.100   0.250   0.200  1.00 10.00           C\nATOM      3  CA  SER A   3       2.200  -0.200   0.100  1.00 10.00           C\nATOM      4  CA  THR A   4       3.300   0.150   0.000  1.00 10.00           C\nEND\n";
    let ply = String::from_utf8(
        convert_to_ply(
            pdb,
            br#"{"format":"pdb","representation":"molstar","center":false,"infer-bonds":false}"#,
        )
        .unwrap(),
    )
    .unwrap();

    assert!(ply_header_count(&ply, "element vertex ") > 0);
    assert!(ply_header_count(&ply, "element face ") > 0);
    assert!(ply.contains("comment molfig_group_count "));
    assert!(ply_vertices_are_finite(&ply));
    assert!(ply_faces_have_valid_indices(&ply));
}

#[test]
fn molstar_ligand_overlay_includes_molstar_connected_whole_residues() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nLIG non-polymer\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A 1 0.000 0.000 0.000\nATOM 2 C CA ALA A 2 1.000 0.000 0.000\nHETATM 3 C C1 LIG L 1 2.000 0.000 0.000\nHETATM 4 O O1 LIG L 1 3.000 0.000 0.000\nHETATM 5 O O HOH W 1 4.000 0.000 0.000\nHETATM 6 Na NA NA I 1 5.000 0.000 0.000\n#\nloop_\n_molstar_bond_site.atom_id_1\n_molstar_bond_site.atom_id_2\n_molstar_bond_site.value_order\n_molstar_bond_site.type_id\n3 4 sing covale\n3 5 sing covale\n3 6 sing covale\n2 3 sing covale\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = render_object_summary_json(&molecule, &options);
    assert_eq!(summary.matches(r#""secondary_type":"ligand""#).count(), 13);
    assert_eq!(summary.matches(r#""visual":"element-sphere""#).count(), 7);
    assert_eq!(summary.matches(r#""visual":"intra-bond""#).count(), 8);
    assert!(summary.contains(r#""secondary_type":"ligand","chain":"L""#));
    assert!(summary.contains(r#""secondary_type":"ligand","chain":"W""#));
    assert!(summary.contains(r#""secondary_type":"ligand","chain":"I""#));
    assert!(summary.contains(r#""secondary_type":"ligand","chain":"A""#));
    assert!(summary.contains(r#""secondary_type":"water","chain":"W""#));
    assert!(summary.contains(r#""secondary_type":"ion","chain":"I""#));
}

#[test]
fn molstar_default_cartoon_selects_nucleotide_ring_only() {
    let mut atoms = Vec::new();
    for seq in 1..=2 {
        let x = seq as f32 * 2.0;
        for (name, offset) in [
            ("P", vec3(x, 0.0, 0.0)),
            ("O3'", vec3(x + 0.2, 0.0, 0.0)),
            ("C4'", vec3(x + 0.2, 0.5, 0.0)),
            ("C3'", vec3(x + 0.3, 0.4, 0.2)),
            ("C1'", vec3(x + 0.1, 0.7, 0.1)),
            ("N1", vec3(x, 1.0, 0.0)),
            ("C2", vec3(x + 0.4, 1.2, 0.0)),
            ("N3", vec3(x + 0.8, 1.0, 0.0)),
            ("C4", vec3(x + 0.8, 0.6, 0.0)),
            ("C5", vec3(x + 0.4, 0.4, 0.0)),
            ("C6", vec3(x, 0.6, 0.0)),
        ] {
            let mut atom = test_atom(atoms.len() + 1, name, "R", seq, offset);
            atom.residue = "A".to_string();
            atoms.push(atom);
        }
    }
    let molecule = Molecule {
        atoms,
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""geometry_type":"nucleotide-ring""#));
    assert!(!summary.contains(r#""geometry_type":"nucleotide-block""#));
    assert!(!summary.contains(r#""geometry_type":"direction-wedge""#));

    let mesh = build_mesh(&molecule, &options);
    assert!(!mesh.faces.is_empty());
    assert_eq!(mesh.faces.len(), mesh.face_groups.len());
    assert!(faces_have_valid_indices(&mesh));
    assert!(mesh.vertices.iter().all(|v| vec3_is_finite(*v)));
}

#[test]
fn backbone_representation_uses_molstar_backbone_cylinder_and_sphere_visuals() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A 1 0.000 0.000 0.000\nATOM 2 C CA ALA A 2 3.800 0.000 0.000\nATOM 3 C CA ALA A 3 7.600 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let options = MeshOptions {
        representation: Representation::Backbone,
        center: false,
        assembly: None,
        infer_bonds: false,
        ..MeshOptions::default()
    };

    let representation = representation_summary_json(&molecule, &options);
    assert!(representation
        .contains(r#""selected_visuals":["polymer-backbone-cylinder","polymer-backbone-sphere"]"#));
    assert!(representation
        .contains(r#""realized_visuals":["polymer-backbone-cylinder","polymer-backbone-sphere"]"#));

    let summary = render_object_summary_json(&molecule, &options);
    assert_eq!(
        summary
            .matches(r#""visual":"polymer-backbone-cylinder""#)
            .count(),
        4
    );
    assert_eq!(
        summary
            .matches(r#""visual":"polymer-backbone-sphere""#)
            .count(),
        3
    );

    let objects = build_render_objects(&molecule, &options);
    let cylinders = objects
        .iter()
        .filter_map(|object| match object {
            RenderObject::Cylinder { start, end, radius } => Some((*start, *end, *radius)),
            _ => None,
        })
        .collect::<Vec<_>>();
    let spheres = objects
        .iter()
        .filter_map(|object| match object {
            RenderObject::Sphere { center, radius } => Some((*center, *radius)),
            _ => None,
        })
        .collect::<Vec<_>>();

    assert_eq!(cylinders.len(), 4);
    assert_eq!(spheres.len(), 3);
    assert_vec3_close(cylinders[0].0, vec3(0.0, 0.0, 0.0), 0.000_001);
    assert_vec3_close(cylinders[0].1, vec3(1.9, 0.0, 0.0), 0.000_001);
    assert_vec3_close(cylinders[1].0, vec3(3.8, 0.0, 0.0), 0.000_001);
    assert_vec3_close(cylinders[1].1, vec3(1.9, 0.0, 0.0), 0.000_001);
    assert!((cylinders[0].2 - 0.3).abs() <= 0.000_001);
    assert_vec3_close(spheres[2].0, vec3(7.6, 0.0, 0.0), 0.000_001);
    assert!((spheres[2].1 - 0.3).abs() <= 0.000_001);
}

#[test]
fn molstar_preset_uses_physical_ball_and_stick_for_ion_and_water() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nHOH non-polymer\nNA ion\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 O O HOH W 1 1 0.000 0.000 0.000\nHETATM 2 Na NA NA I 2 1 3.000 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        infer_bonds: false,
        sphere_detail: 1,
        ..MeshOptions::default()
    };

    let representation = representation_summary_json(&molecule, &options);
    assert!(representation
        .contains(r#""selected_visuals":["element-sphere","intra-bond","inter-bond"]"#));
    assert!(representation.contains(r#""realized_visuals":["element-sphere"]"#));

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""visual":"element-sphere","representation":"molstar","secondary_type":"water","chain":"W","residue_start":1,"residue_end":1,"group_id":0"#));
    assert!(summary.contains(r#""visual":"element-sphere","representation":"molstar","secondary_type":"ion","chain":"I","residue_start":1,"residue_end":1,"group_id":0"#));
    assert!(!summary.contains(r#""representation":"spacefill""#));

    let spheres = build_render_objects(&molecule, &options)
        .into_iter()
        .filter_map(|object| match object {
            RenderObject::Sphere { center, radius } => Some((center, radius)),
            _ => None,
        })
        .collect::<Vec<_>>();

    assert_eq!(spheres.len(), 2);
    assert_vec3_close(spheres[0].0, vec3(3.0, 0.0, 0.0), 0.000_001);
    assert!((spheres[0].1 - 2.27 * 0.15).abs() <= 0.000_001);
    assert_vec3_close(spheres[1].0, vec3(0.0, 0.0, 0.0), 0.000_001);
    assert!((spheres[1].1 - 1.52 * 0.15).abs() <= 0.000_001);
}

#[test]
fn molstar_component_ball_and_stick_uses_physical_bond_radius() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nLIG non-polymer\n#\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\nLIG C1 O1 sing\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 LIG L 1 1 0.000 0.000 0.000\nHETATM 2 O O1 LIG L 1 1 1.400 0.000 0.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        infer_bonds: false,
        sphere_detail: 1,
        ..MeshOptions::default()
    };

    let objects = build_render_objects(&molecule, &options);
    let spheres = objects
        .iter()
        .filter_map(|object| match object {
            RenderObject::Sphere { radius, .. } => Some(*radius),
            _ => None,
        })
        .collect::<Vec<_>>();
    let cylinders = objects
        .iter()
        .filter_map(|object| match object {
            RenderObject::LinkCylinderWithSegments { radius, .. } => Some(*radius),
            _ => None,
        })
        .collect::<Vec<_>>();

    assert_eq!(spheres.len(), 2);
    assert!((spheres[0] - 1.70 * 0.15).abs() <= 0.000_001);
    assert!((spheres[1] - 1.52 * 0.15).abs() <= 0.000_001);
    assert_eq!(cylinders.len(), 2);
    for radius in cylinders {
        assert!((radius - 1.52 * 0.15 * (2.0 / 3.0)).abs() <= 0.000_001);
    }
}

#[test]
fn molstar_component_ball_and_stick_uses_type_symbol_for_physical_size_theme() {
    let mut atom_a = het_atom(1, "C1", "L", 1, "LIG", vec3(0.0, 0.0, 0.0));
    let mut atom_b = het_atom(2, "O1", "L", 1, "LIG", vec3(1.4, 0.0, 0.0));
    atom_a.element = "C".to_string();
    atom_a.type_symbol = "H".to_string();
    atom_b.element = "O".to_string();
    atom_b.type_symbol = "O".to_string();
    let molecule = Molecule {
        atoms: vec![atom_a, atom_b],
        bonds: vec![Bond { a: 0, b: 1 }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        infer_bonds: false,
        sphere_detail: 1,
        ..MeshOptions::default()
    };

    let objects = build_render_objects(&molecule, &options);
    let spheres = objects
        .iter()
        .filter_map(|object| match object {
            RenderObject::Sphere { radius, .. } => Some(*radius),
            _ => None,
        })
        .collect::<Vec<_>>();
    let cylinders = objects
        .iter()
        .filter_map(|object| match object {
            RenderObject::LinkCylinderWithSegments { radius, .. } => Some(*radius),
            _ => None,
        })
        .collect::<Vec<_>>();

    assert_eq!(spheres.len(), 2);
    assert!((spheres[0] - 1.10 * 0.15).abs() <= 0.000_001);
    assert!((spheres[1] - 1.52 * 0.15).abs() <= 0.000_001);
    assert_eq!(cylinders.len(), 2);
    for radius in cylinders {
        assert!((radius - 1.10 * 0.15 * (2.0 / 3.0)).abs() <= 0.000_001);
    }
}

#[test]
fn explicit_cartoon_visuals_realize_nucleotide_block_and_direction_wedge() {
    let named_atoms = [
        ("O3'", vec3(0.0, 0.0, 0.0)),
        ("C4'", vec3(0.2, 0.0, 0.0)),
        ("C3'", vec3(0.9, 0.0, 0.0)),
        ("C1'", vec3(0.2, 0.8, 0.0)),
        ("N1", vec3(0.2, 1.3, 0.0)),
        ("C2", vec3(0.8, 1.5, 0.0)),
        ("N3", vec3(1.4, 1.2, 0.0)),
        ("C4", vec3(1.4, 0.5, 0.0)),
        ("C5", vec3(0.8, 0.2, 0.0)),
        ("C6", vec3(0.2, 0.5, 0.0)),
    ];
    let molecule = nucleotide_molecule("C", &named_atoms);
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        visuals: vec![
            "nucleotide-block".to_string(),
            "direction-wedge".to_string(),
        ],
        ..MeshOptions::default()
    };

    let representation = representation_summary_json(&molecule, &options);
    assert!(representation.contains(r#""selected_visuals":["nucleotide-block","direction-wedge"]"#));
    assert!(representation.contains(r#""realized_visuals":["nucleotide-block","direction-wedge"]"#));

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""geometry_type":"nucleotide-block","visual":"nucleotide-block""#));
    assert!(summary.contains(r#""geometry_type":"direction-wedge","visual":"direction-wedge""#));
    assert!(!summary.contains(r#""geometry_type":"nucleotide-ring""#));

    let mesh = build_mesh(&molecule, &options);
    assert!(!mesh.faces.is_empty());
    assert_eq!(mesh.faces.len(), mesh.face_groups.len());
    assert!(faces_have_valid_indices(&mesh));
    assert!(mesh.vertices.iter().all(|v| vec3_is_finite(*v)));
}

#[test]
fn nucleotide_ring_and_block_use_resolved_quality_detail_and_radial_segments() {
    let Some(molstar_nucleotide_ring) =
        read_molstar_source("mol-repr/structure/visual/nucleotide-ring-mesh.ts")
    else {
        eprintln!("skipping pinned Mol* nucleotide ring source audit; artifacts is absent");
        return;
    };
    let Some(molstar_nucleotide_block) =
        read_molstar_source("mol-repr/structure/visual/nucleotide-block-mesh.ts")
    else {
        eprintln!("skipping pinned Mol* nucleotide block source audit; artifacts is absent");
        return;
    };
    assert!(molstar_nucleotide_ring
        .contains("radialSegments: PD.Numeric(16, { min: 2, max: 56, step: 2 }"));
    assert!(molstar_nucleotide_ring.contains("detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }"));
    assert!(molstar_nucleotide_block
        .contains("radialSegments: PD.Numeric(16, { min: 2, max: 56, step: 2 }"));

    let named_atoms = [
        ("O3'", vec3(0.0, 0.0, 0.0)),
        ("C4'", vec3(0.2, 0.0, 0.0)),
        ("C3'", vec3(0.9, 0.0, 0.0)),
        ("C1'", vec3(0.2, 0.8, 0.0)),
        ("N1", vec3(0.2, 1.3, 0.0)),
        ("C2", vec3(0.8, 1.5, 0.0)),
        ("N3", vec3(1.4, 1.2, 0.0)),
        ("C4", vec3(1.4, 0.5, 0.0)),
        ("C5", vec3(0.8, 0.2, 0.0)),
        ("C6", vec3(0.2, 0.5, 0.0)),
    ];
    let molecule = nucleotide_molecule("C", &named_atoms);
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        quality: Some(VisualQuality::Medium),
        visuals: vec![
            "nucleotide-ring".to_string(),
            "nucleotide-block".to_string(),
        ],
        ..MeshOptions::default()
    };

    let objects = build_render_objects(&molecule, &options);
    let ring_options = objects
        .iter()
        .find_map(|object| match object {
            RenderObject::NucleotideRing {
                detail,
                radial_segments,
                ..
            } => Some((*detail, *radial_segments)),
            _ => None,
        })
        .expect("nucleotide ring");
    let block_radial = objects
        .iter()
        .find_map(|object| match object {
            RenderObject::NucleotideBlock {
                radial_segments, ..
            } => Some(*radial_segments),
            _ => None,
        })
        .expect("nucleotide block");

    assert_eq!(ring_options, (1, 12));
    assert_eq!(block_radial, 12);

    let summary = render_object_summary_json(&molecule, &options);
    eprintln!("{summary}");
    assert!(summary.contains(r#""geometry_type":"nucleotide-ring","visual":"nucleotide-ring""#));
    assert!(summary.contains(r#""value_cell":{"group_id":0,"draw_count":124"#));
    assert!(summary.contains(r#""valueCell":{"drawCount":372,"uVertexCount":128"#));
    assert!(summary.contains(r#""geometry_type":"nucleotide-block","visual":"nucleotide-block""#));
    assert!(summary.contains(r#""value_cell":{"group_id":0,"draw_count":48"#));
    assert!(summary.contains(r#""valueCell":{"drawCount":144,"uVertexCount":75"#));
}

#[test]
fn nucleic_trace_orientation_uses_molstar_polymer_atom_roles() {
    let cif = b"data_demo\nloop_\n_chem_comp.id\n_chem_comp.type\nRCX 'rna linking'\nDCX 'dna linking'\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 O \"O3'\" RCX A 1 0.000 0.000 0.000\nATOM 2 C \"C4'\" RCX A 1 1.000 0.000 0.000\nATOM 3 C \"C3'\" RCX A 1 1.000 2.000 0.000\nATOM 4 C \"C1'\" RCX A 1 4.000 0.000 0.000\nATOM 5 O \"O3'\" DCX A 2 8.000 0.000 0.000\nATOM 6 C \"C4'\" DCX A 2 8.000 5.000 0.000\nATOM 7 C \"C3'\" DCX A 2 8.000 0.000 0.000\nATOM 8 C \"C1'\" DCX A 2 8.000 0.000 3.000\n#\n";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = molecule.atomic_structure();
    let derived = &structure.model.hierarchy.derived.residue;

    assert_eq!(
        derived.polymer_type,
        vec![PolymerType::Rna, PolymerType::Dna]
    );
    assert_eq!(derived.trace_element_index, vec![Some(0), Some(4)]);
    assert_eq!(derived.direction_from_element_index, vec![Some(1), Some(6)]);
    assert_eq!(derived.direction_to_element_index, vec![Some(2), Some(7)]);

    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        infer_bonds: false,
        visuals: vec!["direction-wedge".to_string()],
        ..MeshOptions::default()
    };
    let wedges = build_render_objects(&molecule, &options)
        .into_iter()
        .filter_map(|object| match object {
            RenderObject::DirectionWedge {
                center,
                tangent,
                up,
                ..
            } => Some((center, tangent, up)),
            _ => None,
        })
        .collect::<Vec<_>>();

    assert_eq!(wedges.len(), 2);
    assert_vec3_close(wedges[0].0, vec3(0.0, 0.0, 0.0), 0.000_001);
    assert_vec3_close(wedges[0].1, vec3(1.0, 0.0, 0.0), 0.000_001);
    assert_vec3_close(wedges[0].2, vec3(0.0, 0.788_205_44, 0.615_412_2), 0.000_001);
    assert_vec3_close(wedges[1].0, vec3(8.0, 0.0, 0.0), 0.000_001);
    assert_vec3_close(wedges[1].1, vec3(1.0, 0.0, 0.0), 0.000_001);
    assert_vec3_close(wedges[1].2, vec3(0.0, -0.994_029, -0.109_116_76), 0.000_001);
}

#[test]
fn direction_wedge_omits_molstar_terminal_beta_sheet_residue() {
    let molecule = sheet_test_molecule();
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        visuals: vec!["direction-wedge".to_string()],
        ..MeshOptions::default()
    };

    let objects = build_render_objects(&molecule, &options);
    let wedge_count = objects
        .iter()
        .filter(|object| matches!(object, RenderObject::DirectionWedge { .. }))
        .count();

    assert_eq!(wedge_count, 3);
}

#[test]
fn mesh_options_parse_explicit_visuals_array() {
    let options = MeshOptions::from_json(
        br#"{"representation":"cartoon","visuals":["nucleotide-block","direction-wedge"]}"#,
    )
    .unwrap();
    assert_eq!(options.representation, Representation::Cartoon);
    assert_eq!(
        options.visuals,
        vec![
            "nucleotide-block".to_string(),
            "direction-wedge".to_string()
        ]
    );
}

#[test]
fn surface_and_volume_exclusion_is_documented_and_enforced() {
    let Some(checklist) = read_internal_doc("molstar-parity-checklist.md") else {
        eprintln!("skipping internal surface/volume documentation audit; dev-docs is absent");
        return;
    };
    for snippet in [
        "N/A: gaussian density volume generation is not required",
        "N/A: marching cubes table and algorithm are not required",
        "N/A: molecular surface mesh generation is not required",
        "N/A: gaussian surface mesh generation is not required",
        "N/A: color smoothing metadata is only relevant to excluded",
    ] {
        assert!(
            checklist.contains(snippet),
            "surface/volume checklist should document: {snippet}"
        );
    }

    let Some(api_contract) = read_internal_doc("api-contract.md") else {
        eprintln!("skipping internal API contract audit; dev-docs is absent");
        return;
    };
    for snippet in [
        "Molecular surface and gaussian surface visuals are intentionally excluded",
        "Surface and volume representations such as",
        "enables the Mol* gaussian density",
        "IHM coarse gaussian rows remain supported as coarse units",
    ] {
        assert!(
            api_contract.contains(snippet),
            "API contract should document: {snippet}"
        );
    }

    let readme = include_str!("../../../README.md");
    for snippet in [
        "Molecular surface, gaussian surface, gaussian",
        "density/volume, marching-cubes, and surface color-smoothing metadata",
        "they are not converted into Mol* gaussian",
    ] {
        assert!(
            readme.contains(snippet),
            "README should document: {snippet}"
        );
    }

    for representation in [
        "molecular-surface",
        "gaussian-surface",
        "gaussian-volume",
        "volume",
    ] {
        let json = format!(r#"{{"representation":"{representation}"}}"#);
        assert_eq!(
            MeshOptions::from_json(json.as_bytes()).unwrap_err(),
            format!("unsupported representation: {representation}")
        );
    }

    let molecule = Molecule {
        atoms: vec![test_atom(1, "CA", "A", 1, vec3(0.0, 0.0, 0.0))],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        visuals: vec![
            "molecular-surface".to_string(),
            "gaussian-surface".to_string(),
            "gaussian-volume".to_string(),
            "volume".to_string(),
        ],
        ..MeshOptions::default()
    };

    let summary = representation_summary_json(&molecule, &options);
    assert!(summary.contains(r#""selected_visuals":[]"#));
    assert!(summary.contains(r#""realized_visuals":[]"#));
    let render_summary = render_object_summary_json(&molecule, &options);
    for visual in [
        "molecular-surface",
        "gaussian-surface",
        "gaussian-volume",
        "volume",
    ] {
        assert!(!render_summary.contains(visual));
    }
}

#[test]
fn nucleotide_ring_uses_molstar_named_atom_fallbacks() {
    let named_atoms = [
        ("O3'", vec3(0.0, 0.0, 0.0)),
        ("N1", vec3(0.0, 1.0, 0.0)),
        ("C2", vec3(0.4, 1.2, 0.0)),
        ("N3", vec3(0.8, 1.0, 0.0)),
        ("C4", vec3(0.8, 0.6, 0.0)),
        ("N5", vec3(0.45, 0.35, 0.0)),
        ("C6", vec3(0.0, 0.6, 0.0)),
        ("C7", vec3(0.95, 0.25, 0.0)),
        ("C8", vec3(0.45, 0.0, 0.0)),
        ("N9", vec3(0.1, 0.25, 0.0)),
    ];
    let molecule = nucleotide_molecule("A", &named_atoms);
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let ring = build_render_objects(&molecule, &options)
        .into_iter()
        .find_map(|object| match object {
            RenderObject::NucleotideRing {
                base: Some(base), ..
            } => Some(base),
            _ => None,
        })
        .expect("nucleotide ring base");

    match ring {
        crate::mesh::NucleotideRingBase::Purine { c5, n7, .. } => {
            assert_eq!(c5, vec3(0.45, 0.35, 0.0));
            assert_eq!(n7, vec3(0.95, 0.25, 0.0));
        }
        _ => panic!("expected purine base"),
    }

    let mesh = build_mesh(&molecule, &options);
    assert!(!mesh.faces.is_empty());
    validate_mesh_for_export(&mesh).unwrap();
}

#[test]
fn nucleotide_ring_uses_pyrimidine_c1_fallback() {
    let named_atoms = [
        ("O3'", vec3(0.0, 0.0, 0.0)),
        ("C1", vec3(-0.1, 0.9, 0.0)),
        ("C2", vec3(0.4, 1.2, 0.0)),
        ("N3", vec3(0.8, 1.0, 0.0)),
        ("C4", vec3(0.8, 0.6, 0.0)),
        ("C5", vec3(0.4, 0.4, 0.0)),
        ("C6", vec3(0.0, 0.6, 0.0)),
    ];
    let molecule = nucleotide_molecule("C", &named_atoms);
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let ring = build_render_objects(&molecule, &options)
        .into_iter()
        .find_map(|object| match object {
            RenderObject::NucleotideRing {
                base: Some(base), ..
            } => Some(base),
            _ => None,
        })
        .expect("nucleotide ring base");

    match ring {
        crate::mesh::NucleotideRingBase::Pyrimidine { n1, .. } => {
            assert_eq!(n1, vec3(-0.1, 0.9, 0.0));
        }
        _ => panic!("expected pyrimidine base"),
    }
}

#[test]
fn nucleotide_ring_missing_base_atoms_uses_molstar_connector_not_annulus() {
    let Some(molstar_nucleotide_ring) =
        read_molstar_source("mol-repr/structure/visual/nucleotide-ring-mesh.ts")
    else {
        eprintln!("skipping pinned Mol* nucleotide ring source audit; artifacts is absent");
        return;
    };
    assert!(molstar_nucleotide_ring.contains("if (idx.N1 !== -1 && idx.trace !== -1)"));
    assert!(molstar_nucleotide_ring.contains("if (hasPyrimidineIndices(idx))"));

    let molecule = nucleotide_molecule(
        "C",
        &[
            ("O3'", vec3(0.0, 0.0, 0.0)),
            ("N1", vec3(0.0, 1.0, 0.0)),
            ("C2", vec3(0.4, 1.2, 0.0)),
        ],
    );
    let options = MeshOptions {
        representation: Representation::Cartoon,
        center: false,
        assembly: None,
        sphere_detail: 0,
        radial_segments: 16,
        visuals: vec!["nucleotide-ring".to_string()],
        ..MeshOptions::default()
    };

    let ring = build_render_objects(&molecule, &options)
        .into_iter()
        .find_map(|object| match object {
            RenderObject::NucleotideRing {
                base: Some(base), ..
            } => Some(base),
            _ => None,
        })
        .expect("nucleotide connector");

    match ring {
        crate::mesh::NucleotideRingBase::PyrimidineConnector { trace, n1 } => {
            assert_eq!(trace, vec3(0.0, 0.0, 0.0));
            assert_eq!(n1, vec3(0.0, 1.0, 0.0));
        }
        _ => panic!("expected pyrimidine connector without ring faces"),
    }

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""geometry_type":"nucleotide-ring","visual":"nucleotide-ring""#));
    assert!(
        summary.contains(r#""value_cell":{"group_id":0,"draw_count":52"#),
        "Mol* emits the connector cylinder plus anchor sphere, not a generic annulus: {summary}"
    );
    let mesh = build_mesh(&molecule, &options);
    validate_mesh_for_export(&mesh).unwrap();
}

#[test]
fn nucleotide_base_ring_reference_tests_cover_rna_and_dna_names() {
    let purine_atoms = [
        ("O3'", vec3(0.0, 0.0, 0.0)),
        ("N1", vec3(0.0, 1.0, 0.0)),
        ("C2", vec3(0.4, 1.2, 0.0)),
        ("N3", vec3(0.8, 1.0, 0.0)),
        ("C4", vec3(0.8, 0.6, 0.0)),
        ("C5", vec3(0.45, 0.35, 0.0)),
        ("C6", vec3(0.0, 0.6, 0.0)),
        ("N7", vec3(0.95, 0.25, 0.0)),
        ("C8", vec3(0.45, 0.0, 0.0)),
        ("N9", vec3(0.1, 0.25, 0.0)),
    ];
    let pyrimidine_atoms = [
        ("O3'", vec3(0.0, 0.0, 0.0)),
        ("N1", vec3(0.0, 0.6, 0.0)),
        ("C2", vec3(0.4, 1.2, 0.0)),
        ("N3", vec3(0.8, 1.0, 0.0)),
        ("C4", vec3(0.8, 0.6, 0.0)),
        ("C5", vec3(0.4, 0.4, 0.0)),
        ("C6", vec3(0.0, 0.6, 0.0)),
    ];
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    for residue in ["A", "G", "DA", "DG"] {
        let molecule = nucleotide_molecule(residue, &purine_atoms);
        let ring = build_render_objects(&molecule, &options)
            .into_iter()
            .find_map(|object| match object {
                RenderObject::NucleotideRing {
                    base: Some(base), ..
                } => Some(base),
                _ => None,
            })
            .unwrap_or_else(|| panic!("missing purine ring for {residue}"));
        assert!(
            matches!(ring, crate::mesh::NucleotideRingBase::Purine { .. }),
            "{residue}"
        );
    }

    for residue in ["C", "U", "DC", "DT"] {
        let molecule = nucleotide_molecule(residue, &pyrimidine_atoms);
        let ring = build_render_objects(&molecule, &options)
            .into_iter()
            .find_map(|object| match object {
                RenderObject::NucleotideRing {
                    base: Some(base), ..
                } => Some(base),
                _ => None,
            })
            .unwrap_or_else(|| panic!("missing pyrimidine ring for {residue}"));
        assert!(
            matches!(ring, crate::mesh::NucleotideRingBase::Pyrimidine { .. }),
            "{residue}"
        );
    }
}

#[test]
fn nucleotide_base_geometry_uses_molstar_ring_atoms_only() {
    let named_atoms = [
        ("O3'", vec3(0.0, 0.0, 0.0)),
        ("O6", vec3(-10.0, -10.0, 8.0)),
        ("N6", vec3(10.0, -5.0, -4.0)),
        ("N1", vec3(0.0, 1.0, 0.0)),
        ("C2", vec3(0.4, 1.2, 0.0)),
        ("N3", vec3(0.8, 1.0, 0.0)),
        ("C4", vec3(0.8, 0.6, 0.0)),
        ("C5", vec3(0.45, 0.35, 0.0)),
        ("C6", vec3(0.0, 0.6, 0.0)),
        ("N7", vec3(0.95, 0.25, 0.0)),
        ("C8", vec3(0.45, 0.0, 0.0)),
        ("N9", vec3(0.1, 0.25, 0.0)),
    ];
    let molecule = nucleotide_molecule("A", &named_atoms);
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let (center, normal) = build_render_objects(&molecule, &options)
        .into_iter()
        .find_map(|object| match object {
            RenderObject::NucleotideRing { center, normal, .. } => Some((center, normal)),
            _ => None,
        })
        .expect("nucleotide ring geometry");

    let expected_center = vec3(3.95 / 9.0, 5.25 / 9.0, 0.0);
    assert!(
        center.distance(expected_center) < 0.000_001,
        "exocyclic atoms must not move the base center"
    );
    assert_eq!(
        normal,
        vec3(0.0, 0.0, -1.0),
        "normal follows Mol* triangleNormal(N1, C4, C5)"
    );
}

#[test]
fn nucleic_acid_fixture_covers_rna_and_dna_base_geometry() {
    let fixtures = [
        (
            include_bytes!("../../tests/fixtures/pdb/nucleic-acid-rna-dna.pdb").as_slice(),
            InputFormat::Pdb,
        ),
        (
            include_bytes!("../../tests/fixtures/cif/nucleic-acid-rna-dna.cif").as_slice(),
            InputFormat::Cif,
        ),
    ];
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        infer_bonds: false,
        ..MeshOptions::default()
    };

    for (fixture, format) in fixtures {
        let molecule = parse_molecule_with_options(
            fixture,
            &MeshOptions {
                format,
                assembly: None,
                infer_bonds: false,
                ..MeshOptions::default()
            },
        )
        .unwrap();
        let structure = molecule.atomic_structure();

        assert_eq!(molecule.atoms.len(), 23);
        assert_eq!(
            structure
                .model
                .hierarchy
                .derived
                .molecule_type_count(MoleculeType::Rna),
            1
        );
        assert_eq!(
            structure
                .model
                .hierarchy
                .derived
                .molecule_type_count(MoleculeType::Dna),
            1
        );
        assert_eq!(
            structure
                .units
                .iter()
                .map(|unit| unit.props.nucleotide_elements.len())
                .sum::<usize>(),
            2
        );

        let rings = build_render_objects(&molecule, &options)
            .into_iter()
            .filter_map(|object| match object {
                RenderObject::NucleotideRing {
                    center,
                    normal,
                    base: Some(base),
                    ..
                } => Some((center, normal, base)),
                _ => None,
            })
            .collect::<Vec<_>>();

        assert_eq!(rings.len(), 2);
        assert!(matches!(
            rings[0].2,
            crate::mesh::NucleotideRingBase::Purine { .. }
        ));
        assert!(rings[0].0.distance(vec3(3.95 / 9.0, 5.25 / 9.0, 0.0)) < 0.000_001);
        assert_eq!(rings[0].1, vec3(0.0, 0.0, -1.0));
        assert!(matches!(
            rings[1].2,
            crate::mesh::NucleotideRingBase::Pyrimidine { .. }
        ));
        assert!(rings[1].0.distance(vec3(4.4, 0.8, 0.0)) < 0.000_001);
        assert_eq!(rings[1].1, vec3(0.0, 0.0, -1.0));

        let representation = representation_summary_json(&molecule, &options);
        assert!(representation.contains(r#""nucleotide-ring""#));
        let summary = render_object_summary_json(&molecule, &options);
        assert_eq!(
            summary
                .matches(r#""geometry_type":"nucleotide-ring""#)
                .count(),
            2
        );
        assert!(summary.contains(r#""visual":"nucleotide-ring""#));

        let mesh = build_mesh(&molecule, &options);
        assert!(!mesh.faces.is_empty());
        assert_eq!(mesh.faces.len(), mesh.face_groups.len());
        validate_mesh_for_export(&mesh).unwrap();
    }
}

#[test]
fn rna_fixture_matches_molstar_nucleotide_reference_summary() {
    let expected =
        include_str!("../../tests/expected/nucleotide-rna-dna-reference-summary.json").trim();
    let molecule = parse_molecule_with_options(
        include_bytes!("../../tests/fixtures/cif/nucleic-acid-rna-dna.cif"),
        &MeshOptions {
            format: InputFormat::Cif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let actual = nucleotide_reference_summary_json(&molecule);
    assert_eq!(actual, expected);
    assert!(actual.contains(r#""polymer_type":"rna""#));
    assert!(actual.contains(r#""residue":"A""#));
}

#[test]
fn dna_fixture_matches_molstar_nucleotide_reference_summary() {
    let expected =
        include_str!("../../tests/expected/nucleotide-rna-dna-reference-summary.json").trim();
    let molecule = parse_molecule_with_options(
        include_bytes!("../../tests/fixtures/bcif/nucleic-acid-rna-dna.bcif"),
        &MeshOptions {
            format: InputFormat::BinaryCif,
            assembly: None,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let actual = nucleotide_reference_summary_json(&molecule);
    assert_eq!(actual, expected);
    assert!(actual.contains(r#""polymer_type":"dna""#));
    assert!(actual.contains(r#""residue":"DT""#));
}

#[test]
fn mixed_protein_nucleic_acid_fixture_derives_both_polymer_families() {
    let fixtures = [
        (
            include_bytes!("../../tests/fixtures/pdb/mixed-protein-nucleic.pdb").as_slice(),
            InputFormat::Pdb,
        ),
        (
            include_bytes!("../../tests/fixtures/cif/mixed-protein-nucleic.cif").as_slice(),
            InputFormat::Cif,
        ),
    ];
    let render_options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        infer_bonds: false,
        ..MeshOptions::default()
    };

    for (fixture, format) in fixtures {
        let molecule = parse_molecule_with_options(
            fixture,
            &MeshOptions {
                format,
                assembly: None,
                infer_bonds: false,
                ..MeshOptions::default()
            },
        )
        .unwrap();
        let structure = molecule.atomic_structure();
        let derived = &structure.model.hierarchy.derived;

        assert_eq!(molecule.atoms.len(), 13);
        assert_eq!(molecule.bonds.len(), 2);
        assert_eq!(
            derived.residue.molecule_type,
            vec![MoleculeType::Protein, MoleculeType::Dna]
        );
        assert_eq!(
            derived.residue.polymer_type,
            vec![PolymerType::PeptideL, PolymerType::Dna]
        );
        assert_eq!(derived.molecule_type_count(MoleculeType::Protein), 1);
        assert_eq!(derived.molecule_type_count(MoleculeType::Dna), 1);
        assert_eq!(&derived.atom.is_protein[0..3], [true, true, true]);
        assert!(derived.atom.is_nucleotide[3..]
            .iter()
            .all(|is_nucleotide| *is_nucleotide));
        assert_eq!(
            structure
                .units
                .iter()
                .map(|unit| unit.props.protein_elements.len())
                .sum::<usize>(),
            1
        );
        assert_eq!(
            structure
                .units
                .iter()
                .map(|unit| unit.props.nucleotide_elements.len())
                .sum::<usize>(),
            1
        );

        let rings = build_render_objects(&molecule, &render_options)
            .into_iter()
            .filter_map(|object| match object {
                RenderObject::NucleotideRing {
                    base: Some(base), ..
                } => Some(base),
                _ => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(rings.len(), 1);
        assert!(matches!(
            rings[0],
            crate::mesh::NucleotideRingBase::Pyrimidine { .. }
        ));

        let mesh = build_mesh(&molecule, &render_options);
        assert!(!mesh.faces.is_empty());
        validate_mesh_for_export(&mesh).unwrap();
    }
}

#[test]
fn zero_face_render_objects_do_not_create_phantom_groups() {
    let mut first = test_atom(1, "C1", "A", 1, vec3(0.0, 0.0, 0.0));
    first.residue = "LIG".to_string();
    let mut second = test_atom(2, "C2", "A", 1, vec3(0.0, 0.0, 0.0));
    second.residue = "LIG".to_string();
    let molecule = Molecule {
        atoms: vec![first, second],
        bonds: vec![Bond { a: 0, b: 1 }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::BallAndStick,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let mesh = build_mesh(&molecule, &options);
    assert_eq!(mesh.group_count, 2);
    assert!(mesh
        .face_groups
        .iter()
        .all(|group| *group < mesh.group_count));
    validate_mesh_for_export(&mesh).unwrap();
}

#[test]
fn export_mesh_validation_rejects_nan_and_bad_indices() {
    let mut mesh = grouped_triangle_mesh();
    mesh.vertices[0].x = f32::NAN;
    assert!(validate_mesh_for_export(&mesh)
        .unwrap_err()
        .contains("NaN or infinity"));

    let mut mesh = grouped_triangle_mesh();
    mesh.normals.pop();
    assert!(validate_mesh_for_export(&mesh)
        .unwrap_err()
        .contains("normal count"));

    let mut mesh = grouped_triangle_mesh();
    mesh.faces[0].a = mesh.vertices.len();
    assert!(validate_mesh_for_export(&mesh)
        .unwrap_err()
        .contains("out-of-range"));

    let mut mesh = grouped_triangle_mesh();
    mesh.vertex_groups.pop();
    assert!(validate_mesh_for_export(&mesh)
        .unwrap_err()
        .contains("vertex group count"));

    let mut mesh = grouped_triangle_mesh();
    mesh.face_groups.pop();
    assert!(validate_mesh_for_export(&mesh)
        .unwrap_err()
        .contains("face group count"));

    let mut mesh = grouped_triangle_mesh();
    mesh.face_materials
        .push(crate::model::MeshMaterial::opaque(0x1b9e77));
    assert!(validate_mesh_for_export(&mesh)
        .unwrap_err()
        .contains("face material count"));

    let mut mesh = grouped_triangle_mesh();
    mesh.sections = vec![crate::model::MeshSection {
        key: "broken".to_string(),
        vertex_start: 0,
        vertex_end: mesh.vertices.len() - 1,
        face_start: 0,
        face_end: mesh.faces.len(),
    }];
    assert!(validate_mesh_for_export(&mesh)
        .unwrap_err()
        .contains("sections do not cover"));
}

#[test]
fn render_object_summary_includes_atom_bond_semantics_and_escapes_chain() {
    let molecule = Molecule {
        atoms: vec![
            test_atom(1, "C1", "A\"B", 10, vec3(0.0, 0.0, 0.0)),
            test_atom(2, "O1", "A\"B", 11, vec3(1.0, 0.0, 0.0)),
        ],
        bonds: vec![Bond { a: 0, b: 1 }],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::BallAndStick,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""geometry_type":"sphere","visual":"element-sphere","representation":"ball-and-stick","secondary_type":"atom","chain":"A\"B","residue_start":10,"residue_end":10,"group_id":0"#));
    assert!(summary.contains(r#""geometry_type":"sphere","visual":"element-sphere","representation":"ball-and-stick","secondary_type":"atom","chain":"A\"B","residue_start":11,"residue_end":11,"group_id":1"#));
    assert!(summary.contains(r#""geometry_type":"cylinder","visual":"intra-bond","representation":"ball-and-stick","secondary_type":"bond","chain":"A\"B","residue_start":10,"residue_end":11,"group_id":0"#));
    assert!(summary.contains(r#""value_cell":{"group_id":0"#));
}

#[test]
fn spacefill_summary_uses_molstar_element_sphere_visual_name() {
    let molecule = Molecule {
        atoms: vec![test_atom(1, "O", "A", 10, vec3(0.0, 0.0, 0.0))],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Spacefill,
        center: false,
        assembly: None,
        ..MeshOptions::default()
    };

    let summary = render_object_summary_json(&molecule, &options);
    assert!(summary.contains(r#""geometry_type":"sphere","visual":"element-sphere","representation":"spacefill","secondary_type":"atom","chain":"A","residue_start":10,"residue_end":10,"group_id":0"#));
}

#[test]
fn static_exports_preserve_face_group_metadata() {
    let mesh = grouped_triangle_mesh();

    let obj = export_obj(&mesh);
    assert!(obj.starts_with("mtllib molfig.mtl\n"));
    assert!(!obj.lines().any(|line| line.starts_with('#')));
    assert!(!obj.lines().any(|line| line.starts_with("o ")));
    assert_eq!(obj_group_sequence(&obj), vec![0, 1]);
    assert!(obj.contains("\nv 0 0 0\n"));
    assert!(obj.contains("\nv 1.235 0 0\n"));
    assert!(obj.contains("\nvn 0 0 1\n"));
    assert!(obj.contains("\ng molfig_group_0\nf 1//1 2//2 3//3\n"));
    assert!(obj.contains("\ng molfig_group_1\nf 1//1 3//3 4//4\n"));

    let ply = export_ply(&mesh);
    assert_eq!(
            ply,
            "ply\nformat ascii 1.0\ncomment Exported by molfig\ncomment molfig_group_count 2\ncomment molfig_face_group_property molfig_group\nelement vertex 4\nproperty float x\nproperty float y\nproperty float z\nelement face 2\nproperty list uchar int vertex_indices\nproperty int molfig_group\nend_header\n0.00000 0.00000 0.00000\n1.23456 0.00000 0.00000\n1.00000 1.00000 0.00000\n0.00000 1.00000 0.00000\n3 0 1 2 0\n3 0 2 3 1\n"
        );
    assert!(ply_faces_have_valid_indices(&ply));
    assert_eq!(ply_face_groups(&ply), vec![0, 1]);

    let stl = export_stl(&mesh);
    assert!(stl.starts_with(b"Exported from Mol* 5.9.0"));
    assert_eq!(u32::from_le_bytes(stl[80..84].try_into().unwrap()), 6);
    assert_eq!(stl.len(), 84 + 50 * 6);
    assert_eq!(stl_f32(&stl, 84), 0.0);
    assert_eq!(stl_f32(&stl, 88), 0.0);
    assert_eq!(stl_f32(&stl, 92), 1.0);
    assert!(stl[84..132].iter().any(|byte| *byte != 0));
    assert_eq!(u16::from_le_bytes(stl[132..134].try_into().unwrap()), 0);
    assert!(stl[134..184].iter().all(|byte| *byte == 0));
    assert_eq!(u16::from_le_bytes(stl[182..184].try_into().unwrap()), 0);
    assert!(stl[184..234].iter().all(|byte| *byte == 0));
    assert_eq!(stl_f32(&stl, 234), 0.0);
    assert_eq!(stl_f32(&stl, 238), 0.0);
    assert_eq!(stl_f32(&stl, 242), 1.0);
    assert!(stl[234..282].iter().any(|byte| *byte != 0));
    assert_eq!(u16::from_le_bytes(stl[282..284].try_into().unwrap()), 0);
    assert!(stl[284..334].iter().all(|byte| *byte == 0));
    assert!(stl[334..384].iter().all(|byte| *byte == 0));
}

#[test]
fn api_obj_and_ply_exports_preserve_representable_operator_metadata() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 1 0.0 0.0 0.0\n#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n1 1,2 A\n#\nloop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\n1 1 0 0 0 0 1 0 0 0 0 1 0\n2 1 0 0 5 0 1 0 0 0 0 1 0\n#\n";
    let options =
        br#"{"format":"cif","representation":"spacefill","assembly":"1","center":false,"sphere-detail":1}"#;

    let obj = String::from_utf8(convert_to_obj(cif, options).unwrap()).unwrap();
    assert!(obj.starts_with("mtllib molfig.mtl\n# molfig_operator_metadata "));
    assert!(obj.contains(
        r##"# molfig_operator_metadata {"assembly_id":"1","operator_count":2,"operators":[{"name":"ASM_1","instance_id":"ASM-1","assembly_id":"1","oper_id":1,"oper_list_ids":["1"],"is_identity":true},{"name":"ASM_2","instance_id":"ASM-2","assembly_id":"1","oper_id":2,"oper_list_ids":["2"],"is_identity":false}]}"##
    ));
    assert!(obj.contains("\ng molfig_group_0\n"));

    let ply = String::from_utf8(convert_to_ply(cif, options).unwrap()).unwrap();
    assert!(ply.contains(
        r#"comment molfig_operator_metadata {"assembly_id":"1","operator_count":2,"operators":[{"name":"ASM_1","instance_id":"ASM-1","assembly_id":"1","oper_id":1,"oper_list_ids":["1"],"is_identity":true},{"name":"ASM_2","instance_id":"ASM-2","assembly_id":"1","oper_id":2,"oper_list_ids":["2"],"is_identity":false}]}"#
    ));
    assert!(ply.contains("comment molfig_face_group_property molfig_group\n"));
    assert!(ply_faces_have_valid_indices(&ply));

    let asymmetric_ply = String::from_utf8(
        convert_to_ply(
            cif,
            br#"{"format":"cif","representation":"spacefill","assembly":"asymmetric-unit","center":false,"sphere-detail":1}"#,
        )
        .unwrap(),
    )
    .unwrap();
    assert!(!asymmetric_ply.contains("molfig_operator_metadata"));
}

#[test]
fn api_obj_export_options_can_match_molstar_reference_header_shape() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 1 0.0 0.0 0.0\n#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n1 1 A\n#\nloop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\n1 1 0 0 0 0 1 0 0 0 0 1 0\n#\n";
    let options = br#"{"format":"cif","representation":"spacefill","assembly":"1","center":false,"sphere-detail":1,"obj-basename":"9R1O.mtl","operator-metadata":false,"obj-groups":false}"#;

    let obj = String::from_utf8(convert_to_obj(cif, options).unwrap()).unwrap();
    assert!(obj.starts_with("mtllib 9R1O.mtl\nv "));
    assert!(!obj.lines().any(|line| line.starts_with('#')));
    assert!(!obj.lines().any(|line| line.starts_with("g ")));

    let ply = String::from_utf8(convert_to_ply(cif, options).unwrap()).unwrap();
    assert!(!ply.contains("molfig_operator_metadata"));
    assert!(ply.contains("comment molfig_face_group_property molfig_group\n"));
}

#[test]
fn api_obj_export_emits_molstar_default_materials_from_semantic_objects() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A A 1 1 0.0 0.0 0.0\nATOM 2 O O HOH B B 2 1 4.0 0.0 0.0\n#\n";
    let options =
        br#"{"format":"cif","representation":"spacefill","assembly":"asymmetric-unit","center":false,"sphere-detail":1,"obj-groups":false}"#;

    let parsed_options = MeshOptions::from_json(options).unwrap();
    let molecule = parse_molecule_with_options(cif, &parsed_options).unwrap();
    let mesh = build_mesh(&molecule, &parsed_options);
    assert_eq!(mesh.face_materials.len(), mesh.faces.len());
    assert_eq!(mesh.sections.len(), 1);
    assert_eq!(mesh.sections[0].key, "element-sphere");
    assert_eq!(mesh.sections[0].vertex_start, 0);
    assert_eq!(mesh.sections[0].vertex_end, mesh.vertices.len());
    assert_eq!(mesh.sections[0].face_start, 0);
    assert_eq!(mesh.sections[0].face_end, mesh.faces.len());
    assert!(mesh
        .face_materials
        .iter()
        .any(|material| material.color == 0x1b9e77 && material.alpha_tenths == 10));
    assert!(mesh
        .face_materials
        .iter()
        .any(|material| material.color == 0xff2618 && material.alpha_tenths == 10));

    let obj = String::from_utf8(convert_to_obj(cif, options).unwrap()).unwrap();
    let switches = obj
        .lines()
        .filter(|line| line.starts_with("usemtl "))
        .collect::<Vec<_>>();
    assert_eq!(switches, vec!["usemtl 0x1b9e771", "usemtl 0xff26181"]);
    assert!(obj.contains("\nusemtl 0x1b9e771\nf "));
    assert!(obj.contains("\nusemtl 0xff26181\nf "));

    let mtl = String::from_utf8(convert_to_mtl(cif, options).unwrap()).unwrap();
    if let Some(reference_mtl) = read_repo_file_if_present("package/examples/data/9R1O.mtl") {
        assert_eq!(mtl, reference_mtl);
    }
    assert_eq!(mtl.len(), 181);
    assert_eq!(
        format!("{:016x}", stable_test_hash64(mtl.as_bytes())),
        "6820621c995f3add"
    );

    let material_map = String::from_utf8(maquette_material_map(obj.as_bytes()).unwrap()).unwrap();
    assert_eq!(
        material_map,
        r##"{"0x1b9e771":"#1b9e77","0xff26181":"#ff2618"}"##
    );
}

#[test]
fn api_exports_center_on_visible_renderable_sphere_box_like_molstar_geo_export() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 H H ALA A 1 1 0.0 0.0 0.0\nATOM 2 C C ALA A 1 2 8.0 0.0 0.0\n#\n";
    let raw_options =
        br#"{"format":"cif","representation":"spacefill","assembly":"asymmetric-unit","center":false,"sphere-detail":1}"#;
    let centered_options =
        br#"{"format":"cif","representation":"spacefill","assembly":"asymmetric-unit","center":true,"sphere-detail":1}"#;
    let options = MeshOptions::from_json(raw_options).unwrap();
    let molecule = parse_molecule_with_options(cif, &options).unwrap();
    let (raw_mesh, visible_sphere) = build_mesh_with_visible_bounding_sphere(&molecule, &options);
    validate_mesh_for_export(&raw_mesh).unwrap();

    let sphere = visible_sphere.expect("visible renderable sphere");
    let export_center = if sphere.extrema.len() >= 14 {
        let (min, max) = sphere.extrema.iter().fold(
            (
                vec3(f32::INFINITY, f32::INFINITY, f32::INFINITY),
                vec3(f32::NEG_INFINITY, f32::NEG_INFINITY, f32::NEG_INFINITY),
            ),
            |(mut min, mut max), vertex| {
                min.x = min.x.min(vertex.x);
                min.y = min.y.min(vertex.y);
                min.z = min.z.min(vertex.z);
                max.x = max.x.max(vertex.x);
                max.y = max.y.max(vertex.y);
                max.z = max.z.max(vertex.z);
                (min, max)
            },
        );
        (min + max) * 0.5
    } else {
        sphere.center
    };
    assert!(export_center.distance(Vec3::default()) > 0.1);
    let to_array = |v: Vec3| [v.x, v.y, v.z];
    let expected_first_vertex = raw_mesh.vertices[0] - export_center;

    let raw_obj = String::from_utf8(convert_to_obj(cif, raw_options).unwrap()).unwrap();
    assert_f32_array_close(
        obj_vectors(&raw_obj, "v ")[0],
        to_array(raw_mesh.vertices[0]),
        0.0015,
        "raw OBJ first vertex",
    );

    let obj = String::from_utf8(convert_to_obj(cif, centered_options).unwrap()).unwrap();
    assert_f32_array_close(
        obj_vectors(&obj, "v ")[0],
        to_array(expected_first_vertex),
        0.0015,
        "centered OBJ first vertex",
    );

    let ply = String::from_utf8(convert_to_ply(cif, centered_options).unwrap()).unwrap();
    assert_f32_array_close(
        ply_vertex_samples(&ply)[0].vertex,
        to_array(expected_first_vertex),
        0.00002,
        "centered PLY first vertex",
    );

    let stl = convert_to_stl(cif, centered_options).unwrap();
    let first_stl_vertex = [stl_f32(&stl, 96), stl_f32(&stl, 100), stl_f32(&stl, 104)];
    let expected_first_stl_vertex = raw_mesh.vertices[raw_mesh.faces[0].a] - export_center;
    assert_f32_array_close(
        first_stl_vertex,
        to_array(expected_first_stl_vertex),
        0.00002,
        "centered STL first vertex",
    );
}

#[test]
fn obj_faces_follow_molstar_draw_count_index_order() {
    let mesh = Mesh {
        vertices: vec![
            vec3(0.0, 0.0, 0.0),
            vec3(1.0, 0.0, 0.0),
            vec3(2.0, 0.0, 0.0),
            vec3(3.0, 0.0, 0.0),
            vec3(4.0, 0.0, 0.0),
        ],
        normals: vec![
            vec3(0.0, 0.0, 1.0),
            vec3(0.0, 0.0, 1.0),
            vec3(0.0, 0.0, 1.0),
            vec3(0.0, 0.0, 1.0),
            vec3(0.0, 0.0, 1.0),
        ],
        faces: vec![
            Face { a: 3, b: 1, c: 4 },
            Face { a: 2, b: 4, c: 0 },
            Face { a: 4, b: 3, c: 2 },
        ],
        vertex_groups: vec![7, 7, 7, 3, 3],
        face_groups: vec![7, 7, 3],
        face_materials: Vec::new(),
        sections: Vec::new(),
        group_count: 8,
    };

    let obj = export_obj(&mesh);

    assert_eq!(obj_group_sequence(&obj), vec![7, 3]);
    assert_eq!(
        obj_face_lines(&obj),
        vec!["f 4//4 2//2 5//5", "f 3//3 5//5 1//1", "f 5//5 4//4 3//3",]
    );
}

#[test]
fn static_exports_default_missing_face_groups_to_zero() {
    let mut mesh = grouped_triangle_mesh();
    mesh.face_groups.clear();
    mesh.group_count = 0;

    let obj = export_obj(&mesh);
    assert!(obj.starts_with("mtllib molfig.mtl\n"));
    assert_eq!(obj_group_sequence(&obj), vec![0]);

    let ply = export_ply(&mesh);
    assert!(ply.contains("3 0 1 2 0\n"));
    assert!(ply.contains("3 0 2 3 0\n"));

    let stl = export_stl(&mesh);
    assert!(stl.starts_with(b"Exported from Mol* 5.9.0"));
    assert_eq!(u32::from_le_bytes(stl[80..84].try_into().unwrap()), 6);
    assert_eq!(u16::from_le_bytes(stl[132..134].try_into().unwrap()), 0);
    assert_eq!(u16::from_le_bytes(stl[182..184].try_into().unwrap()), 0);
    assert_eq!(u16::from_le_bytes(stl[282..284].try_into().unwrap()), 0);
}

#[test]
fn exports_obj() {
    let obj =
        String::from_utf8(convert_to_obj(PDB, br#"{"format":"pdb","sphere-detail":1}"#).unwrap())
            .unwrap();
    assert!(obj.contains("\nv "));
    assert!(obj.contains("\nf "));
}

#[test]
fn exports_stl_header_count() {
    let stl = convert_to_stl(PDB, br#"{"format":"pdb","sphere-detail":1}"#).unwrap();
    assert!(stl.starts_with(b"Exported from Mol* 5.9.0"));
    assert!(stl.len() > 84);
    let count = u32::from_le_bytes(stl[80..84].try_into().unwrap());
    assert!(count > 0);
}

#[test]
fn static_exports_are_deterministic_for_same_input() {
    let options =
        br#"{"format":"pdb","representation":"molstar","sphere-detail":1,"center":false}"#;

    assert_eq!(
        convert_to_obj(PDB, options).unwrap(),
        convert_to_obj(PDB, options).unwrap()
    );
    assert_eq!(
        convert_to_stl(PDB, options).unwrap(),
        convert_to_stl(PDB, options).unwrap()
    );
    assert_eq!(
        convert_to_ply(PDB, options).unwrap(),
        convert_to_ply(PDB, options).unwrap()
    );
}

#[test]
fn visible_renderable_sphere_padding_uses_type_symbol_physical_size() {
    let mut atom_a = test_atom(1, "CA", "A", 1, vec3(0.0, 0.0, 0.0));
    let mut atom_b = test_atom(2, "CB", "A", 1, vec3(2.0, 0.0, 0.0));
    atom_a.element = "C".to_string();
    atom_a.type_symbol = "H".to_string();
    atom_b.element = "C".to_string();
    atom_b.type_symbol = "H".to_string();
    let molecule = Molecule {
        atoms: vec![atom_a, atom_b],
        ..Molecule::default()
    };
    let options = MeshOptions {
        representation: Representation::Spacefill,
        center: false,
        assembly: None,
        sphere_detail: 1,
        ..MeshOptions::default()
    };

    let (_, visible_sphere) = build_mesh_with_visible_bounding_sphere(&molecule, &options);
    let sphere = visible_sphere.expect("visible renderable sphere");

    assert!(
        (sphere.radius - 3.736_083).abs() <= 0.000_001,
        "visible sphere radius should use type_symbol vdW padding, got {}",
        sphere.radius
    );
}

#[test]
fn visible_renderable_sphere_padding_tracks_custom_ribbon_tube_radius() {
    let pdb = b"ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C\nATOM      2  CA  GLY A   2       4.000   0.000   0.000  1.00 10.00           C\nEND\n";
    let molecule = parse_molecule(pdb, InputFormat::Pdb).unwrap();
    let default_options = MeshOptions {
        representation: Representation::Ribbon,
        center: false,
        assembly: None,
        sphere_detail: 1,
        ..MeshOptions::default()
    };
    let custom_options = MeshOptions {
        ribbon_radius: 1.0,
        ..default_options.clone()
    };
    let cartoon_options = MeshOptions {
        representation: Representation::Cartoon,
        ribbon_radius: 1.0,
        ..default_options.clone()
    };

    let (_, default_sphere) = build_mesh_with_visible_bounding_sphere(&molecule, &default_options);
    let (_, custom_sphere) = build_mesh_with_visible_bounding_sphere(&molecule, &custom_options);
    let (_, cartoon_sphere) = build_mesh_with_visible_bounding_sphere(&molecule, &cartoon_options);
    let default_radius = default_sphere.expect("default ribbon sphere").radius;
    let custom_radius = custom_sphere.expect("custom ribbon sphere").radius;
    let cartoon_radius = cartoon_sphere.expect("cartoon sphere").radius;

    assert!(
        custom_radius > default_radius + 0.7,
        "custom ribbon tube visible sphere should track the larger custom radius, default={default_radius} custom={custom_radius}"
    );
    assert!(
        (cartoon_radius - default_radius).abs() <= 0.000_01,
        "cartoon visible sphere should keep Mol* uniform sizeFactor padding despite ribbon-radius, default={default_radius} cartoon={cartoon_radius}"
    );
}

#[test]
fn export_diff_reports_clear_pass_and_failure_details() {
    let pass = crate::diff_text("v 0 0 0\n", "v 0 0 0\n", "obj");
    assert!(pass.passed);
    assert_eq!(pass.message, "PASS obj: text match (8 bytes)");

    let fail = crate::diff_text("v 0 0 0\nf 1 2 3\n", "v 0 0 0\nf 1 3 2\n", "obj");
    assert!(!fail.passed);
    assert!(fail
        .message
        .contains("FAIL obj: first difference at line 2"));
    assert!(fail.message.contains(r#"reference="f 1 2 3""#));
    assert!(fail.message.contains(r#"generated="f 1 3 2""#));
    assert!(fail
        .details
        .iter()
        .any(|detail| detail == &("first_line".to_string(), "2".to_string())));

    let binary = crate::diff_bytes(&[0, 1, 2], &[0, 1, 3, 4], "stl");
    assert!(!binary.passed);
    assert!(binary.message.contains("first difference at byte 2"));
    assert!(binary.message.contains("reference=0x02"));
    assert!(binary.message.contains("generated=0x03"));
    assert!(!binary.message.contains("stl_context="));
    assert!(binary
        .details
        .iter()
        .any(|detail| detail == &("first_byte".to_string(), "2".to_string())));

    let mut reference_stl = vec![0u8; 84 + 50];
    let mut generated_stl = vec![0u8; 84 + 50];
    reference_stl[80..84].copy_from_slice(&1u32.to_le_bytes());
    generated_stl[80..84].copy_from_slice(&1u32.to_le_bytes());
    reference_stl[84..88].copy_from_slice(&1.0f32.to_le_bytes());
    generated_stl[84..88].copy_from_slice(&2.0f32.to_le_bytes());
    write_stl_vec3(&mut reference_stl, 96, [0.0, 0.0, 0.0]);
    write_stl_vec3(&mut reference_stl, 108, [1.0, 0.0, 0.0]);
    write_stl_vec3(&mut reference_stl, 120, [0.0, 1.0, 0.0]);
    write_stl_vec3(&mut generated_stl, 96, [1.0, 2.0, 3.0]);
    write_stl_vec3(&mut generated_stl, 108, [2.0, 2.0, 3.0]);
    write_stl_vec3(&mut generated_stl, 120, [1.0, 3.0, 3.0]);

    let stl = crate::diff_bytes(&reference_stl, &generated_stl, "stl");
    assert!(!stl.passed);
    assert!(stl.message.contains("first difference at byte 86"));
    assert!(stl.message.contains("stl_context=facet 0 normal.x byte 2"));
    assert!(stl.message.contains("reference_f32=1"));
    assert!(stl.message.contains("generated_f32=2"));
    assert!(stl
        .message
        .contains("generated_minus_reference={normal:[1,0,0],vertices:[[1,2,3],[1,2,3],[1,2,3]],vertex_centroid:[1,2,3],vertex_residuals:[[0,0,0],[0,0,0],[0,0,0]]}"));
    assert!(stl.message.contains("stl_delta_scan={facet_count:1"));
    assert!(stl.message.contains("nonzero_vertex_delta_facets:1"));
    assert!(stl.message.contains("nonzero_vertex_residual_facets:0"));
    assert!(stl.message.contains("center_like_vertex_delta_facets:1"));
    assert!(stl.message.contains("shape_like_vertex_delta_facets:0"));
    assert!(stl
        .message
        .contains("center_like_residual_ratio_threshold:0.1"));
    assert!(stl
        .message
        .contains("max_vertex_delta_abs:{facet:0,vertex:0,max_component_abs:3,delta:[1,2,3]}"));
    assert!(stl
        .message
        .contains("max_vertex_centroid_abs:{facet:0,max_component_abs:3,delta:[1,2,3]}"));
    assert!(stl
        .message
        .contains("max_vertex_residual_abs:{facet:0,vertex:0,max_component_abs:0,delta:[0,0,0]}"));
    assert!(stl
        .message
        .contains("max_vertex_residual_to_centroid_ratio:{facet:0,ratio:0,centroid_delta:[1,2,3],max_residual_delta:[0,0,0]}"));
    assert!(stl
        .message
        .contains("max_normal_abs:{facet:0,max_component_abs:1,delta:[1,0,0]}"));
    assert!(stl.message.contains("center_fit:{"));
    assert!(stl.message.contains("real_vertex_delta_facets:1"));
    assert!(stl.message.contains("median_vertex_centroid_delta:[1,2,3]"));
    assert!(stl.message.contains("mean_vertex_centroid_delta:[1,2,3]"));
    assert!(stl.message.contains(
        "max_vertex_delta_after_median_center_abs:{facet:0,vertex:0,max_component_abs:0,delta:[0,0,0]}"
    ));
    assert!(stl.message.contains(
        "max_vertex_delta_after_mean_center_abs:{facet:0,vertex:0,max_component_abs:0,delta:[0,0,0]}"
    ));
    assert!(stl
        .details
        .iter()
        .any(|(key, value)| key == "stl_delta_scan"
            && value.contains("facet_count:1")
            && value.contains("center_like_vertex_delta_facets:1")
            && value.contains("center_fit:{")
            && value.contains("max_vertex_residual_to_centroid_ratio:")));
}

fn write_stl_vec3(stl: &mut [u8], offset: usize, values: [f32; 3]) {
    for (axis, value) in values.iter().enumerate() {
        stl[offset + axis * 4..offset + axis * 4 + 4].copy_from_slice(&value.to_le_bytes());
    }
}

#[test]
fn performance_baseline_artifact_covers_large_lookup_mesh_and_wasm_memory() {
    let baseline = include_str!("../../tests/expected/performance-baselines.json");
    assert!(baseline.contains(r#""name": "9r1o-asymmetric-lookup3d""#));
    assert!(baseline.contains(r#""query_count": 64"#));
    assert!(baseline.contains(r#""name": "9r1o-asymmetric-mesh-generation""#));
    assert!(baseline.contains(r#""vertex_count": 74694"#));
    assert!(baseline.contains(r#""face_count": 111504"#));
    assert!(baseline.contains(r#""group_count": 362"#));
    assert!(baseline.contains(r#""debug_max_elapsed_ms": 600000"#));
    assert!(baseline.contains(r#""name": "checked-in-wasm-memory""#));
    assert!(baseline.contains(r#""byte_len": 826019"#));
    assert!(baseline.contains(r#""initial_pages": 18"#));
}

#[test]
fn checked_in_wasm_memory_usage_matches_baseline() {
    let wasm = include_bytes!("../../../package/molfig.wasm");
    assert_eq!(wasm.len(), 826_019);
    assert!(wasm.len() <= 900_000);

    let memory = parse_wasm_memory_summary(wasm).unwrap();
    assert_eq!(
        memory,
        WasmMemorySummary {
            initial_pages: 18,
            maximum_pages: None,
            exported_memory: true,
        }
    );
}

#[test]
#[ignore = "performance baseline; run explicitly when updating large-structure lookup3d numbers"]
fn large_structure_lookup3d_performance_baseline() {
    let pdb = include_bytes!("../../../package/examples/data/9R1O.pdb");
    let options = MeshOptions {
        format: InputFormat::Pdb,
        representation: Representation::Molstar,
        assembly: None,
        sphere_detail: 1,
        ..MeshOptions::default()
    };
    let molecule = parse_molecule_with_options(pdb, &options).unwrap();
    assert_eq!(molecule.atoms.len(), 2870);
    let structure = molecule.atomic_structure();
    let queries = molecule
        .atoms
        .iter()
        .take(64)
        .map(|atom| atom.position)
        .collect::<Vec<_>>();
    assert_eq!(queries.len(), 64);

    let start = std::time::Instant::now();
    let mut total_hits = 0usize;
    for point in queries {
        let hits = structure.lookup3d.find(point, 0.25);
        assert!(!hits.is_empty());
        total_hits += hits.len();
    }
    let elapsed = start.elapsed();

    assert!(total_hits >= 64, "total_hits={total_hits}");
    assert!(
        elapsed <= std::time::Duration::from_millis(1000),
        "lookup3d baseline exceeded: elapsed_ms={}",
        elapsed.as_millis()
    );
}

#[test]
#[ignore = "performance baseline; run explicitly when updating large-structure mesh generation numbers"]
fn large_structure_mesh_generation_performance_baseline() {
    let pdb = include_bytes!("../../../package/examples/data/9R1O.pdb");
    let options = MeshOptions {
        format: InputFormat::Pdb,
        representation: Representation::Molstar,
        assembly: None,
        sphere_detail: 1,
        ..MeshOptions::default()
    };
    let molecule = parse_molecule_with_options(pdb, &options).unwrap();

    let start = std::time::Instant::now();
    let mesh = build_mesh(&molecule, &options);
    let elapsed = start.elapsed();

    validate_mesh_for_export(&mesh).unwrap();
    assert_eq!(mesh.vertices.len(), 74_694);
    assert_eq!(mesh.faces.len(), 111_504);
    assert_eq!(mesh.group_count, 362);
    let max_elapsed = if cfg!(debug_assertions) {
        std::time::Duration::from_millis(600_000)
    } else {
        std::time::Duration::from_millis(5000)
    };
    assert!(
        elapsed <= max_elapsed,
        "mesh generation baseline exceeded: elapsed_ms={}",
        elapsed.as_millis()
    );
}

#[test]
fn topology_derived_bond_sets_iterate_in_canonical_order() {
    let mut molecule = Molecule::default();
    molecule.derived_aromatic_bonds.extend([5, 1, 3, 1]);
    molecule.derived_resonance_bonds.extend([4, 2, 4, 0]);

    assert_eq!(
        molecule
            .derived_aromatic_bonds
            .iter()
            .copied()
            .collect::<Vec<_>>(),
        vec![1, 3, 5]
    );
    assert_eq!(
        molecule
            .derived_resonance_bonds
            .iter()
            .copied()
            .collect::<Vec<_>>(),
        vec![0, 2, 4]
    );
}

#[test]
fn structure_unit_order_is_stable_across_repeated_parse() {
    let cif = b"data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA GLY A 1 0.000 0.000 0.000\nATOM 2 C CA GLY B 1 1.000 0.000 0.000\n#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n1 2 B\n1 1,3 A\n#\nloop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\n1 1 0 0 0 0 1 0 0 0 0 1 0\n2 1 0 0 20 0 1 0 0 0 0 1 0\n3 1 0 0 30 0 1 0 0 0 0 1 0\n#\n";
    let signature = || {
        parse_molecule_with_options(
            cif,
            &MeshOptions {
                format: InputFormat::Cif,
                assembly: Some("1".to_string()),
                infer_bonds: false,
                ..MeshOptions::default()
            },
        )
        .unwrap()
        .atomic_structure()
        .units
        .into_iter()
        .map(|unit| {
            (
                unit.id,
                unit.invariant_id,
                unit.chain_indices,
                unit.elements,
                unit.operator.instance_id,
            )
        })
        .collect::<Vec<_>>()
    };

    let first = signature();
    for _ in 0..5 {
        assert_eq!(signature(), first);
    }
}

#[test]
fn inter_unit_bond_order_is_stable_across_repeated_builds() {
    let assembly = Assembly {
        id: "1".to_string(),
        details: String::new(),
        oligomeric_details: String::new(),
        oligomeric_count: None,
        asym_ids: vec!["A".to_string()],
        transforms: Vec::new(),
        generators: vec![AssemblyGenerator::from_transforms(
            "1",
            vec!["A".to_string()],
            0,
            vec![
                Transform::identity(),
                Transform {
                    m: [
                        [1.0, 0.0, 0.0, 10.0],
                        [0.0, 1.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0, 0.0],
                    ],
                },
            ],
            vec![vec!["1".to_string()], vec!["2".to_string()]],
        )],
    };
    let molecule = Molecule {
        atoms: vec![
            test_atom(1, "C1", "A", 1, vec3(0.0, 0.0, 0.0)),
            test_atom(2, "C2", "A", 1, vec3(1.0, 0.0, 0.0)),
            test_atom(3, "C3", "A", 1, vec3(2.0, 0.0, 0.0)),
        ],
        bonds: vec![Bond { a: 0, b: 1 }, Bond { a: 1, b: 2 }],
        bond_metadata: vec![
            BondMetadata {
                source: BondSource::IndexPair,
                order: 1,
                flags: BondFlags::COVALENT,
                key: 101,
                distance: None,
                operator_a: 1,
                operator_b: 2,
                struct_conn: None,
            },
            BondMetadata {
                source: BondSource::IndexPair,
                order: 1,
                flags: BondFlags::COVALENT,
                key: 102,
                distance: None,
                operator_a: 1,
                operator_b: 2,
                struct_conn: None,
            },
        ],
        selected_assembly: Some(assembly),
        ..Molecule::default()
    };
    let signature = || {
        molecule
            .atomic_structure()
            .inter_unit_bonds
            .into_iter()
            .map(|bond| {
                (
                    bond.unit_a,
                    bond.index_a,
                    bond.unit_b,
                    bond.index_b,
                    bond.source_bond,
                    bond.key,
                )
            })
            .collect::<Vec<_>>()
    };

    let first = signature();
    assert_eq!(first.len(), 2);
    for _ in 0..5 {
        assert_eq!(signature(), first);
    }
}

#[test]
fn inter_unit_bond_graph_construction_canonicalizes_input_bond_order() {
    let assembly = Assembly {
        id: "1".to_string(),
        details: String::new(),
        oligomeric_details: String::new(),
        oligomeric_count: None,
        asym_ids: vec!["A".to_string()],
        transforms: Vec::new(),
        generators: vec![AssemblyGenerator::from_transforms(
            "1",
            vec!["A".to_string()],
            0,
            vec![
                Transform::identity(),
                Transform {
                    m: [
                        [1.0, 0.0, 0.0, 10.0],
                        [0.0, 1.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0, 0.0],
                    ],
                },
                Transform {
                    m: [
                        [1.0, 0.0, 0.0, 20.0],
                        [0.0, 1.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0, 0.0],
                    ],
                },
            ],
            vec![
                vec!["1".to_string()],
                vec!["2".to_string()],
                vec!["3".to_string()],
            ],
        )],
    };
    let atoms = vec![
        test_atom(1, "C1", "A", 1, vec3(0.0, 0.0, 0.0)),
        test_atom(2, "C2", "A", 1, vec3(1.0, 0.0, 0.0)),
        test_atom(3, "C3", "A", 1, vec3(2.0, 0.0, 0.0)),
    ];
    let metadata = [
        BondMetadata {
            source: BondSource::IndexPair,
            order: 1,
            flags: BondFlags::COVALENT,
            key: 9,
            distance: None,
            operator_a: 1,
            operator_b: 3,
            struct_conn: None,
        },
        BondMetadata {
            source: BondSource::IndexPair,
            order: 2,
            flags: BondFlags::COVALENT.union(BondFlags::AROMATIC),
            key: 7,
            distance: None,
            operator_a: 1,
            operator_b: 2,
            struct_conn: None,
        },
    ];
    let build = |bonds: Vec<Bond>, metadata: Vec<BondMetadata>| {
        Molecule {
            atoms: atoms.clone(),
            bonds,
            bond_metadata: metadata,
            selected_assembly: Some(assembly.clone()),
            ..Molecule::default()
        }
        .atomic_structure()
        .inter_unit_bond_graph
        .edges
    };

    let sorted = build(
        vec![Bond { a: 0, b: 1 }, Bond { a: 0, b: 2 }],
        vec![metadata[1].clone(), metadata[0].clone()],
    );
    let shuffled = build(
        vec![Bond { a: 0, b: 2 }, Bond { a: 0, b: 1 }],
        vec![metadata[0].clone(), metadata[1].clone()],
    );

    assert_eq!(shuffled, sorted);
    assert_eq!(
        sorted
            .iter()
            .map(|edge| (
                edge.unit_a,
                edge.index_a,
                edge.unit_b,
                edge.index_b,
                edge.props.key
            ))
            .collect::<Vec<_>>(),
        vec![
            (0, 0, 1, 1, 7),
            (0, 0, 2, 2, 9),
            (1, 1, 0, 0, 7),
            (2, 2, 0, 0, 9),
        ]
    );
}

#[test]
fn carbohydrate_links_are_detected_from_intra_unit_ring_connectivity() {
    let mut molecule = Molecule {
        atoms: vec![
            carbohydrate_atom(1, "C1", "A", 1, "GLC", vec3(0.0, 0.0, 0.0)),
            carbohydrate_atom(2, "C2", "A", 1, "GLC", vec3(1.0, 0.0, 0.0)),
            carbohydrate_atom(3, "C3", "A", 1, "GLC", vec3(1.5, 1.0, 0.0)),
            carbohydrate_atom(4, "C4", "A", 1, "GLC", vec3(1.0, 2.0, 0.0)),
            carbohydrate_atom(5, "C5", "A", 1, "GLC", vec3(0.0, 2.0, 0.0)),
            carbohydrate_atom(6, "O5", "A", 1, "GLC", vec3(-0.5, 1.0, 0.0)),
            carbohydrate_atom(7, "C1", "A", 2, "GLC", vec3(4.0, 0.0, 0.0)),
            carbohydrate_atom(8, "C2", "A", 2, "GLC", vec3(5.0, 0.0, 0.0)),
            carbohydrate_atom(9, "C3", "A", 2, "GLC", vec3(5.5, 1.0, 0.0)),
            carbohydrate_atom(10, "C4", "A", 2, "GLC", vec3(5.0, 2.0, 0.0)),
            carbohydrate_atom(11, "C5", "A", 2, "GLC", vec3(4.0, 2.0, 0.0)),
            carbohydrate_atom(12, "O5", "A", 2, "GLC", vec3(3.5, 1.0, 0.0)),
            carbohydrate_atom(13, "O4", "A", 2, "GLC", vec3(2.5, 1.0, 0.0)),
        ],
        bonds: carbohydrate_bonds(&[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10),
            (10, 11),
            (11, 6),
            (0, 12),
            (12, 9),
        ]),
        ..Molecule::default()
    };
    molecule.atoms[0].element = "H".to_string();
    molecule.atoms[0].type_symbol = "C".to_string();
    molecule.refresh_topology_metadata();

    let structure = molecule.atomic_structure();
    let carbohydrates = &structure.carbohydrates;

    assert_eq!(carbohydrates.elements.len(), 2);
    assert!(carbohydrates.partial_elements.is_empty());
    assert_eq!(
        carbohydrates
            .links
            .iter()
            .map(|link| (link.carbohydrate_index_a, link.carbohydrate_index_b))
            .collect::<Vec<_>>(),
        vec![(0, 1), (1, 0)]
    );
    assert_eq!(structure.carbohydrate_element_indices(0, 0), &[0]);
    assert_eq!(structure.carbohydrate_link_indices(0, 0), &[0]);
    assert_eq!(structure.carbohydrate_link_indices(0, 6), &[1]);
    assert!(
        carbohydrates.elements[0]
            .geometry
            .center
            .distance(vec3(0.5, 1.0, 0.0))
            < 0.000_001
    );
    assert!(
        carbohydrates.elements[1]
            .geometry
            .center
            .distance(vec3(4.5, 1.0, 0.0))
            < 0.000_001
    );
    assert!((carbohydrates.elements[0].geometry.normal.z.abs() - 1.0).abs() < 0.000_001);
    assert_eq!(
        carbohydrates.elements[0].geometry.direction,
        vec3(1.0, 0.0, 0.0)
    );
    assert_eq!(
        carbohydrates.elements[1].geometry.direction,
        vec3(-1.0, 0.0, 0.0)
    );
}

#[test]
fn carbohydrate_link_cylinder_visual_uses_molstar_directed_half_links() {
    let mut molecule = Molecule {
        atoms: vec![
            carbohydrate_atom(1, "C1", "A", 1, "GLC", vec3(0.0, 0.0, 0.0)),
            carbohydrate_atom(2, "C2", "A", 1, "GLC", vec3(1.0, 0.0, 0.0)),
            carbohydrate_atom(3, "C3", "A", 1, "GLC", vec3(1.5, 1.0, 0.0)),
            carbohydrate_atom(4, "C4", "A", 1, "GLC", vec3(1.0, 2.0, 0.0)),
            carbohydrate_atom(5, "C5", "A", 1, "GLC", vec3(0.0, 2.0, 0.0)),
            carbohydrate_atom(6, "O5", "A", 1, "GLC", vec3(-0.5, 1.0, 0.0)),
            carbohydrate_atom(7, "C1", "A", 2, "GLC", vec3(4.0, 0.0, 0.0)),
            carbohydrate_atom(8, "C2", "A", 2, "GLC", vec3(5.0, 0.0, 0.0)),
            carbohydrate_atom(9, "C3", "A", 2, "GLC", vec3(5.5, 1.0, 0.0)),
            carbohydrate_atom(10, "C4", "A", 2, "GLC", vec3(5.0, 2.0, 0.0)),
            carbohydrate_atom(11, "C5", "A", 2, "GLC", vec3(4.0, 2.0, 0.0)),
            carbohydrate_atom(12, "O5", "A", 2, "GLC", vec3(3.5, 1.0, 0.0)),
            carbohydrate_atom(13, "O4", "A", 2, "GLC", vec3(2.5, 1.0, 0.0)),
        ],
        bonds: carbohydrate_bonds(&[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10),
            (10, 11),
            (11, 6),
            (0, 12),
            (12, 9),
        ]),
        ..Molecule::default()
    };
    molecule.refresh_topology_metadata();
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        visuals: vec!["carbohydrate-link".to_string()],
        ..MeshOptions::default()
    };

    let representation = representation_summary_json(&molecule, &options);
    assert!(representation.contains(r#""realized_visuals":["carbohydrate-link"]"#));

    let summary = render_object_summary_json(&molecule, &options);
    assert_eq!(
        summary.matches(r#""visual":"carbohydrate-link""#).count(),
        2
    );
    assert!(summary.contains(r#""visual":"carbohydrate-link","representation":"molstar","secondary_type":"carbohydrate","chain":"A","residue_start":1,"residue_end":2,"group_id":0"#));
    assert!(summary.contains(r#""visual":"carbohydrate-link","representation":"molstar","secondary_type":"carbohydrate","chain":"A","residue_start":2,"residue_end":1,"group_id":1"#));
    assert_eq!(
        summary
            .matches(r#""drawCount":96,"uVertexCount":34"#)
            .count(),
        2
    );

    let link_cylinders = build_render_objects(&molecule, &options)
        .into_iter()
        .filter_map(|object| match object {
            RenderObject::LinkCylinder { start, end, radius } => Some((start, end, radius)),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(link_cylinders.len(), 2);
    assert_vec3_close(link_cylinders[0].0, vec3(0.5, 1.0, 0.0), 0.000_001);
    assert_vec3_close(link_cylinders[0].1, vec3(4.5, 1.0, 0.0), 0.000_001);
    assert_vec3_close(link_cylinders[1].0, vec3(4.5, 1.0, 0.0), 0.000_001);
    assert_vec3_close(link_cylinders[1].1, vec3(0.5, 1.0, 0.0), 0.000_001);
    assert!((link_cylinders[0].2 - 0.51).abs() < 0.000_001);
    assert!((link_cylinders[1].2 - 0.51).abs() < 0.000_001);
}

#[test]
fn carbohydrate_terminal_links_are_detected_from_inter_unit_covalent_bonds() {
    let mut molecule = Molecule {
        atoms: vec![
            carbohydrate_atom(1, "C1", "A", 1, "GLC", vec3(0.0, 0.0, 0.0)),
            carbohydrate_atom(2, "C2", "A", 1, "GLC", vec3(1.0, 0.0, 0.0)),
            carbohydrate_atom(3, "C3", "A", 1, "GLC", vec3(1.5, 1.0, 0.0)),
            carbohydrate_atom(4, "C4", "A", 1, "GLC", vec3(1.0, 2.0, 0.0)),
            carbohydrate_atom(5, "C5", "A", 1, "GLC", vec3(0.0, 2.0, 0.0)),
            carbohydrate_atom(6, "O5", "A", 1, "GLC", vec3(-0.5, 1.0, 0.0)),
            test_atom(7, "ND2", "B", 1, vec3(2.5, 0.0, 0.0)),
        ],
        bonds: carbohydrate_bonds(&[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6)]),
        ..Molecule::default()
    };
    molecule.atoms[0].element = "H".to_string();
    molecule.atoms[0].type_symbol = "C".to_string();
    molecule.refresh_topology_metadata();

    let structure = molecule.atomic_structure();
    let carbohydrates = &structure.carbohydrates;

    assert_eq!(carbohydrates.elements.len(), 1);
    assert!(carbohydrates.links.is_empty());
    assert_eq!(carbohydrates.terminal_links.len(), 2);
    assert!(carbohydrates
        .terminal_links
        .iter()
        .any(|link| link.carbohydrate_index == 0
            && link.element_unit_id == 1
            && link.element_index == 0
            && link.from_carbohydrate));
    assert!(carbohydrates
        .terminal_links
        .iter()
        .any(|link| link.carbohydrate_index == 0
            && link.element_unit_id == 1
            && link.element_index == 0
            && !link.from_carbohydrate));
    assert_eq!(structure.carbohydrate_terminal_link_indices(0, 0), &[0]);
    assert_eq!(structure.carbohydrate_terminal_link_indices(1, 6), &[1]);
    assert_eq!(molecule.carbohydrates().terminal_links.len(), 2);
    assert!(
        carbohydrates.elements[0]
            .geometry
            .center
            .distance(vec3(0.5, 1.0, 0.0))
            < 0.000_001
    );
    assert!(
        carbohydrates.elements[0].geometry.direction.distance(vec3(
            2.0 / 5.0_f32.sqrt(),
            -1.0 / 5.0_f32.sqrt(),
            0.0
        )) < 0.000_001
    );

    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        visuals: vec!["carbohydrate-terminal-link".to_string()],
        ..MeshOptions::default()
    };
    let representation = representation_summary_json(&molecule, &options);
    assert!(representation.contains(r#""realized_visuals":["carbohydrate-terminal-link"]"#));

    let summary = render_object_summary_json(&molecule, &options);
    assert_eq!(
        summary
            .matches(r#""visual":"carbohydrate-terminal-link""#)
            .count(),
        2
    );
    assert!(summary.contains(r#""visual":"carbohydrate-terminal-link","representation":"molstar","secondary_type":"carbohydrate","chain":"A","residue_start":1,"residue_end":1,"group_id":0"#));
    assert!(summary.contains(r#""visual":"carbohydrate-terminal-link","representation":"molstar","secondary_type":"carbohydrate","chain":"B","residue_start":1,"residue_end":1,"group_id":1"#));
    assert_eq!(
        summary
            .matches(r#""drawCount":96,"uVertexCount":34"#)
            .count(),
        2
    );

    let terminal_cylinders = build_render_objects(&molecule, &options)
        .into_iter()
        .filter_map(|object| match object {
            RenderObject::LinkCylinder { start, end, radius } => Some((start, end, radius)),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(terminal_cylinders.len(), 2);
    assert_vec3_close(terminal_cylinders[0].0, vec3(0.5, 1.0, 0.0), 0.000_001);
    assert_vec3_close(terminal_cylinders[0].1, vec3(2.5, 0.0, 0.0), 0.000_001);
    assert_vec3_close(terminal_cylinders[1].0, vec3(2.5, 0.0, 0.0), 0.000_001);
    assert_vec3_close(terminal_cylinders[1].1, vec3(0.5, 1.0, 0.0), 0.000_001);
    assert!((terminal_cylinders[0].2 - 0.34).abs() < 0.000_001);
    assert!((terminal_cylinders[1].2 - 0.31).abs() < 0.000_001);
}

#[test]
fn mesh_group_order_is_stable_across_repeated_builds() {
    let molecule = parse_molecule(PDB, InputFormat::Pdb).unwrap();
    let options = MeshOptions {
        representation: Representation::BallAndStick,
        center: false,
        assembly: None,
        sphere_detail: 1,
        ..MeshOptions::default()
    };
    let signature = || {
        let mesh = build_mesh(&molecule, &options);
        (mesh.group_count, mesh.face_groups)
    };

    let first = signature();
    assert!(first.0 > 0);
    assert!(!first.1.is_empty());
    for _ in 0..5 {
        assert_eq!(signature(), first);
    }
}

fn vec3_is_finite(v: Vec3) -> bool {
    v.x.is_finite() && v.y.is_finite() && v.z.is_finite()
}

fn mesh_bounds(mesh: &Mesh) -> (Vec3, Vec3) {
    let first = mesh.vertices[0];
    mesh.vertices
        .iter()
        .copied()
        .fold((first, first), |(min, max), vertex| {
            (min.min(vertex), max.max(vertex))
        })
}

fn ply_header_count(ply: &str, prefix: &str) -> usize {
    ply.lines()
        .find_map(|line| line.strip_prefix(prefix))
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0)
}

fn ply_header_comment_usize(ply: &str, prefix: &str) -> usize {
    ply.lines()
        .take_while(|line| *line != "end_header")
        .find_map(|line| line.strip_prefix(prefix))
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0)
}

fn stable_test_hash64(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

fn contract_value<'a>(contract: &'a str, key: &str) -> &'a str {
    contract_optional_value(contract, key).unwrap_or_else(|| panic!("missing contract key {key}"))
}

fn contract_optional_value<'a>(contract: &'a str, key: &str) -> Option<&'a str> {
    let prefix = format!("{key}=");
    contract.lines().find_map(|line| line.strip_prefix(&prefix))
}

fn contract_usize(contract: &str, key: &str) -> usize {
    contract_value(contract, key)
        .parse::<usize>()
        .unwrap_or_else(|_| panic!("invalid usize contract key {key}"))
}

fn contract_f32(contract: &str, key: &str) -> f32 {
    contract_value(contract, key)
        .parse::<f32>()
        .unwrap_or_else(|_| panic!("invalid f32 contract key {key}"))
}

fn contract_usize_array4(contract: &str, key: &str) -> [usize; 4] {
    let values = contract_value(contract, key)
        .split(',')
        .map(|value| {
            value
                .parse::<usize>()
                .unwrap_or_else(|_| panic!("invalid usize array contract key {key}"))
        })
        .collect::<Vec<_>>();
    values
        .try_into()
        .unwrap_or_else(|_| panic!("expected four values for contract key {key}"))
}

fn contract_f32_array(contract: &str, key: &str) -> [f32; 3] {
    let values = contract_value(contract, key)
        .split(',')
        .map(|value| {
            value
                .parse::<f32>()
                .unwrap_or_else(|_| panic!("invalid f32 array contract key {key}"))
        })
        .collect::<Vec<_>>();
    values
        .try_into()
        .unwrap_or_else(|_| panic!("expected three values for contract key {key}"))
}

fn assert_f32_array_close(actual: [f32; 3], expected: [f32; 3], tolerance: f32, label: &str) {
    for axis in 0..3 {
        assert!(
            (actual[axis] - expected[axis]).abs() <= tolerance,
            "{label}[{axis}] differs: actual={} expected={} tolerance={}",
            actual[axis],
            expected[axis],
            tolerance
        );
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ReferenceManifestEntry<'a> {
    contract_path: &'a str,
    formats: Option<Vec<&'a str>>,
}

#[derive(Debug)]
struct ParityEvidenceEntry {
    domain: String,
    unit_tests: Vec<String>,
    fixtures: Vec<String>,
    reference_comparisons: Vec<String>,
    notes: String,
}

#[derive(Debug)]
struct ParitySectionEvidence {
    domains: Vec<String>,
    checked_items: usize,
}

#[derive(Debug)]
struct ParityMeshReferenceEvidence {
    domains: Vec<String>,
    checked_items: usize,
    notes: String,
}

impl ReferenceManifestEntry<'_> {
    fn includes_format(&self, format: &str) -> bool {
        self.formats
            .as_ref()
            .map(|formats| formats.contains(&format))
            .unwrap_or(true)
    }
}

fn reference_manifest_entries(manifest: &str) -> Vec<&str> {
    reference_manifest_entries_with_formats(manifest)
        .into_iter()
        .map(|entry| entry.contract_path)
        .collect()
}

fn reference_manifest_entries_for_format<'a>(manifest: &'a str, format: &str) -> Vec<&'a str> {
    reference_manifest_entries_for_formats(manifest, &[format])
}

fn reference_manifest_entries_for_formats<'a>(manifest: &'a str, formats: &[&str]) -> Vec<&'a str> {
    validate_reference_manifest_formats(formats, "reference_manifest_entries_for_formats");
    reference_manifest_entries_with_formats(manifest)
        .into_iter()
        .filter(|entry| formats.iter().all(|format| entry.includes_format(format)))
        .map(|entry| entry.contract_path)
        .collect()
}

fn reference_manifest_entries_with_formats(manifest: &str) -> Vec<ReferenceManifestEntry<'_>> {
    manifest
        .lines()
        .filter_map(reference_manifest_entry)
        .collect()
}

fn parity_evidence_entries() -> Vec<ParityEvidenceEntry> {
    let manifest = read_manifest_file("tests/expected/molstar-parity-evidence.tsv");
    let mut lines = manifest
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.trim_start().starts_with('#'));
    let header = lines
        .next()
        .expect("missing Mol* parity evidence manifest header");
    assert_eq!(
        header, "domain\tunit_tests\tfixtures\treference_comparisons\tnotes",
        "unexpected Mol* parity evidence manifest header"
    );

    lines
        .enumerate()
        .map(|(index, line)| {
            let fields = line.split('\t').collect::<Vec<_>>();
            assert_eq!(
                fields.len(),
                5,
                "invalid Mol* parity evidence row {}: expected 5 tab-separated fields",
                index + 2
            );
            ParityEvidenceEntry {
                domain: fields[0].to_string(),
                unit_tests: split_evidence_list(fields[1]),
                fixtures: split_evidence_list(fields[2]),
                reference_comparisons: split_evidence_list(fields[3]),
                notes: fields[4].to_string(),
            }
        })
        .collect()
}

fn parity_section_evidence_domains() -> std::collections::BTreeMap<String, ParitySectionEvidence> {
    let manifest = read_manifest_file("tests/expected/molstar-parity-section-evidence.tsv");
    let mut lines = manifest
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.trim_start().starts_with('#'));
    let header = lines
        .next()
        .expect("missing Mol* parity section evidence manifest header");
    assert_eq!(
        header, "section\tdomains\tchecked_items",
        "unexpected Mol* parity section evidence manifest header"
    );

    let mut out = std::collections::BTreeMap::new();
    for (index, line) in lines.enumerate() {
        let fields = line.split('\t').collect::<Vec<_>>();
        assert_eq!(
            fields.len(),
            3,
            "invalid Mol* parity section evidence row {}: expected 3 tab-separated fields",
            index + 2
        );
        let checked_items = fields[2].parse::<usize>().unwrap_or_else(|_| {
            panic!(
                "invalid checked_items in Mol* parity section evidence row {}",
                index + 2
            )
        });
        let previous = out.insert(
            fields[0].to_string(),
            ParitySectionEvidence {
                domains: split_evidence_list(fields[1]),
                checked_items,
            },
        );
        assert!(
            previous.is_none(),
            "duplicate Mol* parity section evidence row for {}",
            fields[0]
        );
    }
    out
}

fn parity_mesh_reference_evidence(
) -> std::collections::BTreeMap<String, ParityMeshReferenceEvidence> {
    let manifest = read_manifest_file("tests/expected/molstar-parity-mesh-reference-evidence.tsv");
    let mut lines = manifest
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.trim_start().starts_with('#'));
    let header = lines
        .next()
        .expect("missing Mol* mesh reference evidence manifest header");
    assert_eq!(
        header, "mesh_section\tmolstar_reference_domains\tchecked_items\tnotes",
        "unexpected Mol* mesh reference evidence manifest header"
    );

    let mut out = std::collections::BTreeMap::new();
    for (index, line) in lines.enumerate() {
        let fields = line.split('\t').collect::<Vec<_>>();
        assert_eq!(
            fields.len(),
            4,
            "invalid Mol* mesh reference evidence row {}: expected 4 tab-separated fields",
            index + 2
        );
        let checked_items = fields[2].parse::<usize>().unwrap_or_else(|_| {
            panic!(
                "invalid checked_items in Mol* mesh reference evidence row {}",
                index + 2
            )
        });
        let previous = out.insert(
            fields[0].to_string(),
            ParityMeshReferenceEvidence {
                domains: split_evidence_list(fields[1]),
                checked_items,
                notes: fields[3].to_string(),
            },
        );
        assert!(
            previous.is_none(),
            "duplicate Mol* mesh reference evidence row for {}",
            fields[0]
        );
    }
    out
}

fn checklist_checked_items_by_section(
    checklist: &str,
) -> std::collections::BTreeMap<String, usize> {
    let mut out = std::collections::BTreeMap::new();
    let mut current_section = "";
    for line in checklist.lines() {
        if let Some(section) = line.strip_prefix("## ") {
            current_section = section.trim();
        } else if line.trim_start().starts_with("- [x]") {
            *out.entry(current_section.to_string()).or_default() += 1;
        }
    }
    out
}

fn split_evidence_list(value: &str) -> Vec<String> {
    value
        .split('|')
        .filter_map(|item| {
            let item = item.trim();
            (!item.is_empty() && item != "-").then(|| item.to_string())
        })
        .collect()
}

fn parity_evidence_test_sources() -> String {
    [
        "src/tests/mod.rs",
        "src/tests/molstar_model_parity.rs",
        "src/parser/tests.rs",
        "src/mesh/geometry.rs",
        "src/export/tests.rs",
        "tests/molfig_diff_cli.rs",
    ]
    .into_iter()
    .map(read_manifest_file)
    .collect::<Vec<_>>()
    .join("\n")
}

fn assert_manifest_token_exists(
    repo: &std::path::Path,
    test_sources: &str,
    domain: &str,
    field: &str,
    token: &str,
) {
    if looks_like_repo_path(token) {
        let path = repo.join(token);
        assert!(
            path.exists(),
            "{field} listed for {domain} does not exist: {token}"
        );
    } else {
        assert!(
            test_sources.contains(&format!("fn {token}(")),
            "{field} listed for {domain} is neither an existing path nor a known test: {token}"
        );
    }
}

fn looks_like_repo_path(token: &str) -> bool {
    token.contains('/')
        || token.ends_with(".wasm")
        || token.ends_with(".typ")
        || token.ends_with(".md")
        || token.ends_with(".json")
        || token.ends_with(".tsv")
        || token.ends_with(".txt")
        || token.ends_with(".contract")
}

fn reference_manifest_entry(line: &str) -> Option<ReferenceManifestEntry<'_>> {
    let line = line
        .split_once('#')
        .map_or(line, |(before, _)| before)
        .trim();
    if line.is_empty() {
        return None;
    }

    let mut bare_path = None;
    let mut contract_path = None;
    let mut formats = None;
    for token in line.split_whitespace() {
        if let Some((key, value)) = token.split_once('=') {
            assert!(
                !value.is_empty(),
                "manifest field '{key}' requires a value: {line}"
            );
            match key {
                "contract" | "path" => {
                    assert!(
                        contract_path.is_none() && bare_path.is_none(),
                        "duplicate reference manifest contract path: {line}"
                    );
                    contract_path = Some(value);
                }
                "formats" => {
                    assert!(
                        formats.is_none(),
                        "duplicate formats= manifest field: {line}"
                    );
                    let parsed = value.split(',').map(str::trim).collect::<Vec<_>>();
                    validate_reference_manifest_formats(&parsed, "formats");
                    formats = Some(parsed);
                }
                "tag" | "tags" | "note" => {}
                _ => panic!("unknown reference manifest field '{key}': {line}"),
            }
        } else {
            assert!(
                bare_path.is_none() && contract_path.is_none(),
                "duplicate reference manifest contract path: {line}"
            );
            bare_path = Some(token);
        }
    }

    let contract_path = contract_path
        .or(bare_path)
        .unwrap_or_else(|| panic!("reference manifest entry is missing a contract path: {line}"));
    Some(ReferenceManifestEntry {
        contract_path,
        formats,
    })
}

fn validate_reference_manifest_formats(formats: &[&str], label: &str) {
    assert!(!formats.is_empty(), "{label}: expected at least one format");
    let mut seen = Vec::new();
    for format in formats {
        assert!(
            matches!(*format, "json" | "obj" | "stl"),
            "{label}: unsupported format '{format}'"
        );
        assert!(
            !seen.contains(format),
            "{label}: duplicate format '{format}'"
        );
        seen.push(*format);
    }
}

fn read_manifest_file(relative_path: &str) -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(relative_path);
    std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()))
}

fn read_manifest_file_bytes(relative_path: &str) -> Vec<u8> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(relative_path);
    std::fs::read(&path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()))
}

fn faces_have_valid_indices(mesh: &Mesh) -> bool {
    mesh.faces.iter().all(|face| {
        face.a < mesh.vertices.len() && face.b < mesh.vertices.len() && face.c < mesh.vertices.len()
    })
}

fn ply_vertices_are_finite(ply: &str) -> bool {
    let vertex_count = ply_header_count(ply, "element vertex ");
    let mut lines = ply.lines().skip_while(|line| *line != "end_header").skip(1);
    for _ in 0..vertex_count {
        let Some(line) = lines.next() else {
            return false;
        };
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() != 3 {
            return false;
        }
        for value in fields {
            let Ok(number) = value.parse::<f32>() else {
                return false;
            };
            if !number.is_finite() {
                return false;
            }
        }
    }
    true
}

fn ply_faces_have_valid_indices(ply: &str) -> bool {
    let vertex_count = ply_header_count(ply, "element vertex ");
    let face_count = ply_header_count(ply, "element face ");
    let mut lines = ply.lines().skip_while(|line| *line != "end_header").skip(1);
    for _ in 0..vertex_count {
        if lines.next().is_none() {
            return false;
        }
    }

    let mut observed_faces = 0usize;
    for line in lines {
        let fields: Vec<&str> = line.split_whitespace().collect();
        let Ok(index_count) = fields.first().unwrap_or(&"").parse::<usize>() else {
            return false;
        };
        if index_count != 3 || fields.len() < 1 + index_count {
            return false;
        }
        for value in &fields[1..1 + index_count] {
            let Ok(index) = value.parse::<usize>() else {
                return false;
            };
            if index >= vertex_count {
                return false;
            }
        }
        observed_faces += 1;
    }
    observed_faces == face_count
}

fn ply_vertex_bounds(ply: &str) -> ([f32; 3], [f32; 3]) {
    let vertex_count = ply_header_count(ply, "element vertex ");
    let mut min = [f32::INFINITY; 3];
    let mut max = [f32::NEG_INFINITY; 3];
    for line in ply
        .lines()
        .skip_while(|line| *line != "end_header")
        .skip(1)
        .take(vertex_count)
    {
        for (axis, value) in line.split_whitespace().take(3).enumerate() {
            let value = value.parse::<f32>().unwrap();
            min[axis] = min[axis].min(value);
            max[axis] = max[axis].max(value);
        }
    }
    (min, max)
}

fn contract_triplet(values: [f32; 3]) -> String {
    format!("{:.4},{:.4},{:.4}", values[0], values[1], values[2])
}

fn ply_face_groups(ply: &str) -> Vec<usize> {
    let vertex_count = ply_header_count(ply, "element vertex ");
    ply.lines()
        .skip_while(|line| *line != "end_header")
        .skip(1 + vertex_count)
        .filter_map(|line| line.split_whitespace().last()?.parse::<usize>().ok())
        .collect()
}

fn obj_group_sequence(obj: &str) -> Vec<usize> {
    obj.lines()
        .filter_map(|line| line.strip_prefix("g molfig_group_"))
        .filter_map(|group| group.parse::<usize>().ok())
        .collect()
}

fn obj_face_lines(obj: &str) -> Vec<&str> {
    obj.lines().filter(|line| line.starts_with("f ")).collect()
}

fn obj_face_indices(obj: &str) -> Vec<[usize; 3]> {
    obj_face_lines(obj)
        .into_iter()
        .map(|line| {
            let fields = line.split_whitespace().skip(1).collect::<Vec<_>>();
            assert_eq!(fields.len(), 3, "expected triangular OBJ face: {line}");
            [
                obj_face_vertex_index(fields[0]),
                obj_face_vertex_index(fields[1]),
                obj_face_vertex_index(fields[2]),
            ]
        })
        .collect()
}

fn obj_face_vertex_index(value: &str) -> usize {
    value
        .split('/')
        .next()
        .expect("missing OBJ vertex index")
        .parse::<usize>()
        .expect("invalid OBJ vertex index")
        - 1
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct ObjStats {
    vertex_count: usize,
    normal_count: usize,
    face_count: usize,
    min: [f32; 3],
    max: [f32; 3],
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct StlStats {
    byte_len: usize,
    facet_count: usize,
    min: [f32; 3],
    max: [f32; 3],
}

fn obj_stats(obj: &str) -> ObjStats {
    let mut stats = ObjStats {
        vertex_count: 0,
        normal_count: 0,
        face_count: 0,
        min: [f32::INFINITY; 3],
        max: [f32::NEG_INFINITY; 3],
    };
    for line in obj.lines() {
        if let Some(rest) = line.strip_prefix("v ") {
            stats.vertex_count += 1;
            for (axis, value) in rest.split_whitespace().take(3).enumerate() {
                let value = value.parse::<f32>().unwrap();
                stats.min[axis] = stats.min[axis].min(value);
                stats.max[axis] = stats.max[axis].max(value);
            }
        } else if line.starts_with("vn ") {
            stats.normal_count += 1;
        } else if line.starts_with("f ") {
            stats.face_count += 1;
        }
    }
    stats
}

fn stl_stats(stl: &[u8]) -> StlStats {
    assert!(
        stl.len() >= 84,
        "STL must contain an 80-byte header and count"
    );
    let facet_count = u32::from_le_bytes(stl[80..84].try_into().unwrap()) as usize;
    assert_eq!(stl.len(), 84 + facet_count * 50);
    let mut stats = StlStats {
        byte_len: stl.len(),
        facet_count,
        min: [f32::INFINITY; 3],
        max: [f32::NEG_INFINITY; 3],
    };
    for facet in 0..facet_count {
        let base = 84 + facet * 50 + 12;
        for vertex in 0..3 {
            for axis in 0..3 {
                let offset = base + (vertex * 3 + axis) * 4;
                let value = stl_f32(stl, offset);
                stats.min[axis] = stats.min[axis].min(value);
                stats.max[axis] = stats.max[axis].max(value);
            }
        }
    }
    stats
}

fn nine_r1o_default_assembly_live_export_diff_json() -> Option<String> {
    const FIXTURE: &str = "package/examples/data/9R1O.pdb";
    const OPTIONS: &str = r#"{"format":"pdb","representation":"molstar","sphere-detail":1}"#;
    let fixture = include_bytes!("../../../package/examples/data/9R1O.pdb");
    let reference_obj = read_repo_file_if_present("package/examples/data/9R1O.obj")?;
    let reference_stl = read_repo_bytes_if_present("package/examples/data/9R1O.stl")?;
    let live_obj = String::from_utf8(convert_to_obj(fixture, OPTIONS.as_bytes()).unwrap())
        .expect("live OBJ must be UTF-8");
    let live_stl = convert_to_stl(fixture, OPTIONS.as_bytes()).unwrap();

    let live_obj_stats = obj_stats(&live_obj);
    let reference_obj_stats = obj_stats(&reference_obj);
    let live_stl_stats = stl_stats(&live_stl);
    let reference_stl_stats = stl_stats(&reference_stl);
    let live_obj_vertices = obj_vectors(&live_obj, "v ");
    let live_obj_normals = obj_vectors(&live_obj, "vn ");
    let reference_obj_vertices = obj_vectors(&reference_obj, "v ");
    let reference_obj_normals = obj_vectors(&reference_obj, "vn ");

    Some(format!(
        concat!(
            "{{\n",
            "  \"schema\": 1,\n",
            "  \"name\": \"9r1o-default-assembly-live-export-diff\",\n",
            "  \"molstar_reference_commit\": \"{}\",\n",
            "  \"fixture\": \"{}\",\n",
            "  \"options\": {},\n",
            "  \"molfig_live\": {{\n",
            "    \"source\": \"live molfig export\",\n",
            "    \"obj\": {},\n",
            "    \"stl\": {}\n",
            "  }},\n",
            "  \"molstar_reference\": {{\n",
            "    \"source\": \"pinned Mol* geo-export\",\n",
            "    \"obj\": {},\n",
            "    \"stl\": {}\n",
            "  }},\n",
            "  \"diff\": {{\n",
            "    \"obj_counts_delta\": {{\"vertex\": {}, \"normal\": {}, \"face\": {}}},\n",
            "    \"obj_bounds_abs_delta\": {{\"min\": {}, \"max\": {}}},\n",
            "    \"obj_sample_abs_delta\": [{}],\n",
            "    \"stl_counts_delta\": {{\"facet\": {}}},\n",
            "    \"stl_bounds_abs_delta\": {{\"min\": {}, \"max\": {}}},\n",
            "    \"stl_sample_abs_delta\": [{}],\n",
            "    \"stl_sample_signed_delta\": [{}],\n",
            "    \"stl_first_facet_signed_delta\": {}\n",
            "  }}\n",
            "}}"
        ),
        MOLSTAR_REFERENCE_COMMIT,
        FIXTURE,
        OPTIONS,
        obj_export_geometry_json(live_obj_stats, &live_obj_vertices, &live_obj_normals),
        stl_export_geometry_json(live_stl_stats, &live_stl),
        obj_export_geometry_json(
            reference_obj_stats,
            &reference_obj_vertices,
            &reference_obj_normals,
        ),
        stl_export_geometry_json(reference_stl_stats, &reference_stl),
        count_delta(
            live_obj_stats.vertex_count,
            reference_obj_stats.vertex_count
        ),
        count_delta(
            live_obj_stats.normal_count,
            reference_obj_stats.normal_count
        ),
        count_delta(live_obj_stats.face_count, reference_obj_stats.face_count),
        f32_triplet_json(abs_delta_array(live_obj_stats.min, reference_obj_stats.min)),
        f32_triplet_json(abs_delta_array(live_obj_stats.max, reference_obj_stats.max)),
        obj_sample_abs_deltas_json(
            &live_obj_vertices,
            &live_obj_normals,
            &reference_obj_vertices,
            &reference_obj_normals,
        ),
        count_delta(live_stl_stats.facet_count, reference_stl_stats.facet_count),
        f32_triplet_json(abs_delta_array(live_stl_stats.min, reference_stl_stats.min)),
        f32_triplet_json(abs_delta_array(live_stl_stats.max, reference_stl_stats.max)),
        stl_sample_abs_deltas_json(
            &live_stl,
            live_stl_stats.facet_count,
            &reference_stl,
            reference_stl_stats.facet_count,
        ),
        stl_sample_signed_deltas_json(
            &live_stl,
            live_stl_stats.facet_count,
            &reference_stl,
            reference_stl_stats.facet_count,
        ),
        stl_facet_signed_delta_json(&live_stl, 0, &reference_stl, 0)
    ))
}

fn nine_r1o_asymmetric_ply_vs_assembly_one_reference_gap_json() -> Option<String> {
    const FIXTURE: &str = "package/examples/data/9R1O.pdb";
    const MOLFIG_OPTIONS: &str = r#"{"format":"pdb","representation":"molstar","assembly":"asymmetric-unit","sphere-detail":1}"#;
    const MOLSTAR_REFERENCE_OPTIONS: &str =
        r#"{"format":"pdb","representation":"molstar","assembly":"1","sphere-detail":1}"#;
    let reference_obj = read_repo_file_if_present("package/examples/data/9R1O.obj")?;
    let reference_stl = read_repo_bytes_if_present("package/examples/data/9R1O.stl")?;
    let ply_contract =
        include_str!("../../tests/expected/ply/9r1o-molstar-asymmetric-unit.ply.contract");

    let reference_obj_stats = obj_stats(&reference_obj);
    let reference_stl_stats = stl_stats(&reference_stl);
    let molfig_vertex_count = contract_usize(ply_contract, "vertex_count");
    let molfig_face_count = contract_usize(ply_contract, "face_count");
    let molfig_bounds_min = contract_f32_array(ply_contract, "bounds_min");
    let molfig_bounds_max = contract_f32_array(ply_contract, "bounds_max");
    let molfig_vertex_samples = contract_vertex_samples(ply_contract);
    let molfig_stl_facet_count = molfig_face_count * 3;
    let molfig_stl_byte_len = 84 + molfig_stl_facet_count * 50;

    assert_eq!(reference_obj_stats.vertex_count, 101_114);
    assert_eq!(
        reference_obj_stats.normal_count,
        reference_obj_stats.vertex_count
    );
    assert_eq!(reference_obj_stats.face_count, 178_864);
    assert_eq!(
        reference_stl_stats.facet_count,
        reference_obj_stats.face_count * 3
    );
    assert_ne!(
        molfig_vertex_count, reference_obj_stats.vertex_count,
        "9R1O cartoon vertex-count parity was reached; update this snapshot and checklist"
    );
    assert_ne!(
        molfig_face_count, reference_obj_stats.face_count,
        "9R1O cartoon face-count parity was reached; update this snapshot and checklist"
    );

    let reference_obj_vertices = obj_vectors(&reference_obj, "v ");
    let reference_obj_normals = obj_vectors(&reference_obj, "vn ");

    Some(format!(
        concat!(
            "{{\n",
            "  \"schema\": 1,\n",
            "  \"name\": \"9r1o-asymmetric-ply-vs-assembly-1-reference-gap\",\n",
            "  \"molstar_reference_commit\": \"{}\",\n",
            "  \"fixture\": \"{}\",\n",
            "  \"comparison_scope\": \"diagnostic only: package-owned asymmetric-unit PLY contract vs pinned Mol* assembly 1 OBJ/STL reference\",\n",
            "  \"molfig_options\": {},\n",
            "  \"molstar_reference_options\": {},\n",
            "  \"tolerances\": {{\"positions_and_normals\": 0.0001, \"bounds\": 0.0001}},\n",
            "  \"molstar_reference\": {{\n",
            "    \"obj\": {},\n",
            "    \"stl\": {}\n",
            "  }},\n",
            "  \"molfig_current\": {{\n",
            "    \"obj\": {{\"vertex_count\":{},\"normal_count\":{},\"face_count\":{},\"bounds\":{}}},\n",
            "    \"ply\": {{\"vertex_count\":{},\"face_count\":{},\"bounds\":{},\"samples\":[{}],\"face_samples\":[{}],\"contract\":\"tests/expected/ply/9r1o-molstar-asymmetric-unit.ply.contract\"}},\n",
            "    \"stl\": {{\"byte_len\":{},\"facet_count\":{},\"bounds\":{}}}\n",
            "  }},\n",
            "  \"parity_gaps\": {{\n",
            "    \"obj_vertex_count_delta\": {},\n",
            "    \"obj_face_count_delta\": {},\n",
            "    \"stl_facet_count_delta\": {},\n",
            "    \"obj_bounds_abs_delta\": {{\"min\": {}, \"max\": {}}},\n",
            "    \"stl_bounds_abs_delta\": {{\"min\": {}, \"max\": {}}},\n",
            "    \"obj_vertex_sample_abs_delta_vs_ply\": [{}],\n",
            "    \"ply_reference\": \"Mol* geo-export has no PLY exporter; molfig PLY is package-owned\"\n",
            "  }},\n",
            "  \"blockers\": [\n",
            "    \"current molfig OBJ/PLY vertex_count=99954 and face_count=166128 do not match Mol* OBJ vertex_count=101114 and face_count=178864\",\n",
            "    \"current molfig STL facet_count=498384 does not match Mol* STL facet_count=536592\",\n",
            "    \"sampled OBJ vertices/normals and STL facet vertices/normals are pinned here as regression data until tube/ribbon/sheet/trace geometry reaches exact parity\"\n",
            "  ]\n",
            "}}"
        ),
        MOLSTAR_REFERENCE_COMMIT,
        FIXTURE,
        MOLFIG_OPTIONS,
        MOLSTAR_REFERENCE_OPTIONS,
        obj_export_geometry_json(
            reference_obj_stats,
            &reference_obj_vertices,
            &reference_obj_normals
        ),
        stl_export_geometry_json(reference_stl_stats, &reference_stl),
        molfig_vertex_count,
        molfig_vertex_count,
        molfig_face_count,
        bounds_json(molfig_bounds_min, molfig_bounds_max),
        molfig_vertex_count,
        molfig_face_count,
        bounds_json(molfig_bounds_min, molfig_bounds_max),
        vertex_samples_json(&molfig_vertex_samples),
        contract_value(ply_contract, "face_samples"),
        molfig_stl_byte_len,
        molfig_stl_facet_count,
        bounds_json(molfig_bounds_min, molfig_bounds_max),
        count_delta(molfig_vertex_count, reference_obj_stats.vertex_count),
        count_delta(molfig_face_count, reference_obj_stats.face_count),
        count_delta(molfig_stl_facet_count, reference_stl_stats.facet_count),
        f32_triplet_json(abs_delta_array(molfig_bounds_min, reference_obj_stats.min)),
        f32_triplet_json(abs_delta_array(molfig_bounds_max, reference_obj_stats.max)),
        f32_triplet_json(abs_delta_array(molfig_bounds_min, reference_stl_stats.min)),
        f32_triplet_json(abs_delta_array(molfig_bounds_max, reference_stl_stats.max)),
        vertex_sample_abs_deltas_json(&molfig_vertex_samples, &reference_obj_vertices)
    ))
}

fn assert_ratio_in_range(
    actual: usize,
    reference: usize,
    min_ratio: f32,
    max_ratio: f32,
    label: &str,
) {
    let ratio = actual as f32 / reference as f32;
    assert!(
        ratio >= min_ratio && ratio <= max_ratio,
        "{label} ratio out of range: actual={actual}, reference={reference}, ratio={ratio}"
    );
}

fn stl_f32(stl: &[u8], offset: usize) -> f32 {
    f32::from_le_bytes(stl[offset..offset + 4].try_into().unwrap())
}

const EXPORT_SAMPLE_FRACTIONS: [(&str, f32); 5] = [
    ("start", 0.0),
    ("quarter", 0.25),
    ("middle", 0.5),
    ("three_quarter", 0.75),
    ("end", 1.0),
];

#[derive(Clone, Copy, Debug, PartialEq)]
struct VertexSample {
    label: &'static str,
    index: usize,
    vertex: [f32; 3],
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct FaceSample {
    label: &'static str,
    index: usize,
    indices: [usize; 3],
    group: usize,
}

fn parse_vec3_fields<'a, I>(mut fields: I) -> [f32; 3]
where
    I: Iterator<Item = &'a str>,
{
    [
        fields
            .next()
            .expect("missing x coordinate")
            .parse::<f32>()
            .expect("invalid x coordinate"),
        fields
            .next()
            .expect("missing y coordinate")
            .parse::<f32>()
            .expect("invalid y coordinate"),
        fields
            .next()
            .expect("missing z coordinate")
            .parse::<f32>()
            .expect("invalid z coordinate"),
    ]
}

fn obj_vectors(obj: &str, prefix: &str) -> Vec<[f32; 3]> {
    obj.lines()
        .filter_map(|line| line.strip_prefix(prefix))
        .map(|rest| parse_vec3_fields(rest.split_whitespace()))
        .collect()
}

fn obj_export_geometry_json(
    stats: ObjStats,
    vertices: &[[f32; 3]],
    normals: &[[f32; 3]],
) -> String {
    assert_eq!(vertices.len(), stats.vertex_count);
    assert_eq!(normals.len(), stats.normal_count);
    format!(
        "{{\"vertex_count\":{},\"normal_count\":{},\"face_count\":{},\"bounds\":{},\"samples\":[{}]}}",
        stats.vertex_count,
        stats.normal_count,
        stats.face_count,
        bounds_json(stats.min, stats.max),
        obj_samples_json(vertices, normals)
    )
}

fn stl_export_geometry_json(stats: StlStats, stl: &[u8]) -> String {
    format!(
        "{{\"byte_len\":{},\"facet_count\":{},\"bounds\":{},\"samples\":[{}]}}",
        stats.byte_len,
        stats.facet_count,
        bounds_json(stats.min, stats.max),
        stl_samples_json(stl, stats.facet_count)
    )
}

fn obj_samples_json(vertices: &[[f32; 3]], normals: &[[f32; 3]]) -> String {
    EXPORT_SAMPLE_FRACTIONS
        .iter()
        .map(|(label, fraction)| {
            let index = sample_index(vertices.len(), *fraction);
            format!(
                "{{\"label\":\"{}\",\"index\":{},\"vertex\":{},\"normal\":{}}}",
                label,
                index,
                f32_triplet_json(vertices[index]),
                f32_triplet_json(normals[index])
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn stl_samples_json(stl: &[u8], facet_count: usize) -> String {
    EXPORT_SAMPLE_FRACTIONS
        .iter()
        .map(|(label, fraction)| {
            let facet_index = sample_index(facet_count, *fraction);
            format!(
                "{{\"label\":\"{}\",\"facet_index\":{},\"normal\":{},\"vertices\":{}}}",
                label,
                facet_index,
                f32_triplet_json(stl_facet_normal(stl, facet_index)),
                stl_facet_vertices_json(stl, facet_index)
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn obj_sample_abs_deltas_json(
    actual_vertices: &[[f32; 3]],
    actual_normals: &[[f32; 3]],
    reference_vertices: &[[f32; 3]],
    reference_normals: &[[f32; 3]],
) -> String {
    EXPORT_SAMPLE_FRACTIONS
        .iter()
        .map(|(label, fraction)| {
            let actual_index = sample_index(actual_vertices.len(), *fraction);
            let reference_index = sample_index(reference_vertices.len(), *fraction);
            format!(
                "{{\"label\":\"{}\",\"molfig_index\":{},\"reference_index\":{},\"vertex_abs_delta\":{},\"normal_abs_delta\":{}}}",
                label,
                actual_index,
                reference_index,
                f32_triplet_json(abs_delta_array(
                    actual_vertices[actual_index],
                    reference_vertices[reference_index],
                )),
                f32_triplet_json(abs_delta_array(
                    actual_normals[actual_index],
                    reference_normals[reference_index],
                ))
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn stl_sample_abs_deltas_json(
    actual_stl: &[u8],
    actual_facet_count: usize,
    reference_stl: &[u8],
    reference_facet_count: usize,
) -> String {
    EXPORT_SAMPLE_FRACTIONS
        .iter()
        .map(|(label, fraction)| {
            let actual_index = sample_index(actual_facet_count, *fraction);
            let reference_index = sample_index(reference_facet_count, *fraction);
            format!(
                "{{\"label\":\"{}\",\"molfig_facet_index\":{},\"reference_facet_index\":{},\"normal_abs_delta\":{},\"vertices_abs_delta\":{}}}",
                label,
                actual_index,
                reference_index,
                f32_triplet_json(abs_delta_array(
                    stl_facet_normal(actual_stl, actual_index),
                    stl_facet_normal(reference_stl, reference_index),
                )),
                stl_facet_vertices_abs_delta_json(
                    actual_stl,
                    actual_index,
                    reference_stl,
                    reference_index,
                )
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn stl_sample_signed_deltas_json(
    actual_stl: &[u8],
    actual_facet_count: usize,
    reference_stl: &[u8],
    reference_facet_count: usize,
) -> String {
    EXPORT_SAMPLE_FRACTIONS
        .iter()
        .map(|(label, fraction)| {
            let actual_index = sample_index(actual_facet_count, *fraction);
            let reference_index = sample_index(reference_facet_count, *fraction);
            format!(
                "{{\"label\":\"{}\",\"molfig_facet_index\":{},\"reference_facet_index\":{},\"normal_signed_delta\":{},\"vertices_signed_delta\":{}}}",
                label,
                actual_index,
                reference_index,
                precise_f32_triplet_json(signed_delta_array(
                    stl_facet_normal(actual_stl, actual_index),
                    stl_facet_normal(reference_stl, reference_index),
                )),
                stl_facet_vertices_signed_delta_json(
                    actual_stl,
                    actual_index,
                    reference_stl,
                    reference_index,
                )
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn stl_facet_signed_delta_json(
    actual_stl: &[u8],
    actual_facet_index: usize,
    reference_stl: &[u8],
    reference_facet_index: usize,
) -> String {
    let vertex_deltas = [
        signed_delta_array(
            stl_facet_vertex(actual_stl, actual_facet_index, 0),
            stl_facet_vertex(reference_stl, reference_facet_index, 0),
        ),
        signed_delta_array(
            stl_facet_vertex(actual_stl, actual_facet_index, 1),
            stl_facet_vertex(reference_stl, reference_facet_index, 1),
        ),
        signed_delta_array(
            stl_facet_vertex(actual_stl, actual_facet_index, 2),
            stl_facet_vertex(reference_stl, reference_facet_index, 2),
        ),
    ];
    let centroid = [
        (vertex_deltas[0][0] + vertex_deltas[1][0] + vertex_deltas[2][0]) / 3.0,
        (vertex_deltas[0][1] + vertex_deltas[1][1] + vertex_deltas[2][1]) / 3.0,
        (vertex_deltas[0][2] + vertex_deltas[1][2] + vertex_deltas[2][2]) / 3.0,
    ];
    let residuals = vertex_deltas
        .iter()
        .map(|delta| {
            precise_f32_triplet_json([
                delta[0] - centroid[0],
                delta[1] - centroid[1],
                delta[2] - centroid[2],
            ])
        })
        .collect::<Vec<_>>()
        .join(",");

    format!(
        "{{\"molfig_facet_index\":{},\"reference_facet_index\":{},\"normal_signed_delta\":{},\"vertices_signed_delta\":[{}],\"vertex_centroid_signed_delta\":{},\"vertex_residual_signed_delta\":[{}]}}",
        actual_facet_index,
        reference_facet_index,
        precise_f32_triplet_json(signed_delta_array(
            stl_facet_normal(actual_stl, actual_facet_index),
            stl_facet_normal(reference_stl, reference_facet_index),
        )),
        vertex_deltas
            .iter()
            .map(|delta| precise_f32_triplet_json(*delta))
            .collect::<Vec<_>>()
            .join(","),
        precise_f32_triplet_json(centroid),
        residuals
    )
}

fn stl_facet_normal(stl: &[u8], facet_index: usize) -> [f32; 3] {
    let base = 84 + facet_index * 50;
    [
        stl_f32(stl, base),
        stl_f32(stl, base + 4),
        stl_f32(stl, base + 8),
    ]
}

fn stl_facet_vertex(stl: &[u8], facet_index: usize, vertex_index: usize) -> [f32; 3] {
    let base = 84 + facet_index * 50 + 12 + vertex_index * 12;
    [
        stl_f32(stl, base),
        stl_f32(stl, base + 4),
        stl_f32(stl, base + 8),
    ]
}

fn stl_facet_vertices_json(stl: &[u8], facet_index: usize) -> String {
    let vertices = (0..3)
        .map(|vertex_index| f32_triplet_json(stl_facet_vertex(stl, facet_index, vertex_index)))
        .collect::<Vec<_>>()
        .join(",");
    format!("[{vertices}]")
}

fn stl_facet_vertices_abs_delta_json(
    actual_stl: &[u8],
    actual_facet_index: usize,
    reference_stl: &[u8],
    reference_facet_index: usize,
) -> String {
    let vertices = (0..3)
        .map(|vertex_index| {
            f32_triplet_json(abs_delta_array(
                stl_facet_vertex(actual_stl, actual_facet_index, vertex_index),
                stl_facet_vertex(reference_stl, reference_facet_index, vertex_index),
            ))
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("[{vertices}]")
}

fn stl_facet_vertices_signed_delta_json(
    actual_stl: &[u8],
    actual_facet_index: usize,
    reference_stl: &[u8],
    reference_facet_index: usize,
) -> String {
    let vertices = (0..3)
        .map(|vertex_index| {
            precise_f32_triplet_json(signed_delta_array(
                stl_facet_vertex(actual_stl, actual_facet_index, vertex_index),
                stl_facet_vertex(reference_stl, reference_facet_index, vertex_index),
            ))
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("[{vertices}]")
}

fn ply_vertex_samples_json(ply: &str) -> String {
    vertex_samples_json(&ply_vertex_samples(ply))
}

fn ply_vertex_samples(ply: &str) -> Vec<VertexSample> {
    let vertex_count = ply_header_count(ply, "element vertex ");
    let vertices = ply
        .lines()
        .skip_while(|line| *line != "end_header")
        .skip(1)
        .take(vertex_count)
        .map(|line| parse_vec3_fields(line.split_whitespace()))
        .collect::<Vec<_>>();
    EXPORT_SAMPLE_FRACTIONS
        .iter()
        .map(|(label, fraction)| {
            let index = sample_index(vertices.len(), *fraction);
            VertexSample {
                label,
                index,
                vertex: vertices[index],
            }
        })
        .collect()
}

fn ply_face_samples_json(ply: &str) -> String {
    face_samples_json(&ply_face_samples(ply))
}

fn ply_face_samples(ply: &str) -> Vec<FaceSample> {
    let vertex_count = ply_header_count(ply, "element vertex ");
    let face_count = ply_header_count(ply, "element face ");
    let faces = ply
        .lines()
        .skip_while(|line| *line != "end_header")
        .skip(1 + vertex_count)
        .take(face_count)
        .map(parse_ply_face_sample)
        .collect::<Vec<_>>();
    EXPORT_SAMPLE_FRACTIONS
        .iter()
        .map(|(label, fraction)| {
            let index = sample_index(faces.len(), *fraction);
            let (indices, group) = faces[index];
            FaceSample {
                label,
                index,
                indices,
                group,
            }
        })
        .collect()
}

fn parse_ply_face_sample(line: &str) -> ([usize; 3], usize) {
    let fields = line
        .split_whitespace()
        .map(|field| field.parse::<usize>().expect("invalid PLY face field"))
        .collect::<Vec<_>>();
    assert_eq!(fields.first(), Some(&3), "expected triangular PLY face");
    assert!(
        fields.len() >= 5,
        "expected PLY face indices followed by molfig_group"
    );
    ([fields[1], fields[2], fields[3]], fields[4])
}

fn vertex_samples_json(samples: &[VertexSample]) -> String {
    samples
        .iter()
        .map(|sample| {
            format!(
                "{{\"label\":\"{}\",\"index\":{},\"vertex\":{}}}",
                sample.label,
                sample.index,
                f32_triplet_json(sample.vertex)
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn face_samples_json(samples: &[FaceSample]) -> String {
    samples
        .iter()
        .map(|sample| {
            format!(
                "{{\"label\":\"{}\",\"index\":{},\"indices\":[{},{},{}],\"group\":{}}}",
                sample.label,
                sample.index,
                sample.indices[0],
                sample.indices[1],
                sample.indices[2],
                sample.group
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn contract_vertex_samples(contract: &str) -> Vec<VertexSample> {
    EXPORT_SAMPLE_FRACTIONS
        .iter()
        .map(|(label, _)| VertexSample {
            label,
            index: contract_usize(contract, &format!("sample_{label}_vertex_index")),
            vertex: contract_f32_array(contract, &format!("sample_{label}_vertex")),
        })
        .collect()
}

fn vertex_sample_abs_deltas_json(
    samples: &[VertexSample],
    reference_vertices: &[[f32; 3]],
) -> String {
    samples
        .iter()
        .zip(EXPORT_SAMPLE_FRACTIONS.iter())
        .map(|(sample, (label, fraction))| {
            assert_eq!(&sample.label, label);
            let reference_index = sample_index(reference_vertices.len(), *fraction);
            format!(
                "{{\"label\":\"{}\",\"molfig_index\":{},\"reference_index\":{},\"abs_delta\":{}}}",
                sample.label,
                sample.index,
                reference_index,
                f32_triplet_json(abs_delta_array(
                    sample.vertex,
                    reference_vertices[reference_index],
                ))
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn normal_sample_abs_deltas_json(
    actual_normals: &[[f32; 3]],
    reference_normals: &[[f32; 3]],
) -> String {
    EXPORT_SAMPLE_FRACTIONS
        .iter()
        .map(|(label, fraction)| {
            let actual_index = sample_index(actual_normals.len(), *fraction);
            let reference_index = sample_index(reference_normals.len(), *fraction);
            format!(
                "{{\"label\":\"{}\",\"molfig_index\":{},\"reference_index\":{},\"abs_delta\":{}}}",
                label,
                actual_index,
                reference_index,
                f32_triplet_json(abs_delta_array(
                    actual_normals[actual_index],
                    reference_normals[reference_index],
                ))
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn sample_index(count: usize, fraction: f32) -> usize {
    assert!(count > 0, "cannot sample an empty geometry array");
    (((count - 1) as f32) * fraction).round() as usize
}

fn bounds_json(min: [f32; 3], max: [f32; 3]) -> String {
    format!(
        "{{\"min\":{},\"max\":{}}}",
        f32_triplet_json(min),
        f32_triplet_json(max)
    )
}

fn f32_triplet_json(values: [f32; 3]) -> String {
    format!(
        "[{},{},{}]",
        test_f32_json(values[0]),
        test_f32_json(values[1]),
        test_f32_json(values[2])
    )
}

fn precise_f32_triplet_json(values: [f32; 3]) -> String {
    format!("[{:.8},{:.8},{:.8}]", values[0], values[1], values[2])
}

fn abs_delta_array(actual: [f32; 3], reference: [f32; 3]) -> [f32; 3] {
    [
        (actual[0] - reference[0]).abs(),
        (actual[1] - reference[1]).abs(),
        (actual[2] - reference[2]).abs(),
    ]
}

fn signed_delta_array(actual: [f32; 3], reference: [f32; 3]) -> [f32; 3] {
    [
        actual[0] - reference[0],
        actual[1] - reference[1],
        actual[2] - reference[2],
    ]
}

fn count_delta(actual: usize, reference: usize) -> i64 {
    actual as i64 - reference as i64
}

fn pretty_f32_array_json(values: [f32; 3]) -> String {
    format!(
        "[{}, {}, {}]",
        test_f32_json(values[0]),
        test_f32_json(values[1]),
        test_f32_json(values[2])
    )
}

fn molstar_reference_name(contract_path: &str) -> String {
    std::path::Path::new(contract_path)
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or(contract_path)
        .trim_end_matches(".reference.contract")
        .to_string()
}

fn molstar_reference_summary_from_contract(contract_path: &str, contract: &str) -> String {
    let obj_path = contract_value(contract, "obj_reference");
    let stl_path = contract_value(contract, "stl_reference");
    let options_json = contract_value(contract, "options");
    let options = MeshOptions::from_json(options_json.as_bytes())
        .unwrap_or_else(|error| panic!("{contract_path}: invalid options JSON: {error}"));
    let input_format = molstar_reference_input_format(contract, &options);
    let representation = molstar_reference_representation(options.representation);
    let assembly = options.assembly.as_deref().unwrap_or("asymmetric-unit");
    let sphere_detail = molstar_reference_sphere_detail(options_json);
    let obj = read_manifest_file(obj_path);
    let stl = read_manifest_file_bytes(stl_path);
    let obj_stats = obj_stats(&obj);
    let stl_stats = stl_stats(&stl);
    let stl_header = std::str::from_utf8(&stl[..80])
        .unwrap()
        .trim_end_matches('\0');

    format!(
        concat!(
            "{{\n",
            "  \"schema\": 1,\n",
            "  \"name\": \"{}\",\n",
            "  \"molstar_reference_commit\": \"{}\",\n",
            "  \"fixture\": \"{}\",\n",
            "  \"options\": {{\n",
            "    \"format\": \"{}\",\n",
            "    \"representation\": \"{}\",\n",
            "    \"assembly\": \"{}\",\n",
            "    \"sphere-detail\": {}\n",
            "  }},\n",
            "  \"source\": {{\n",
            "    \"obj\": \"{}\",\n",
            "    \"stl\": \"{}\"\n",
            "  }},\n",
            "  \"obj\": {{\n",
            "    \"byte_len\": {},\n",
            "    \"fnv1a64\": \"{}\",\n",
            "    \"vertex_count\": {},\n",
            "    \"normal_count\": {},\n",
            "    \"face_count\": {},\n",
            "    \"bounds\": {{\n",
            "      \"min\": {},\n",
            "      \"max\": {}\n",
            "    }}\n",
            "  }},\n",
            "  \"stl\": {{\n",
            "    \"byte_len\": {},\n",
            "    \"fnv1a64\": \"{}\",\n",
            "    \"facet_count\": {},\n",
            "    \"header\": \"{}\",\n",
            "    \"bounds\": {{\n",
            "      \"min\": {},\n",
            "      \"max\": {}\n",
            "    }}\n",
            "  }}\n",
            "}}\n"
        ),
        molstar_reference_name(contract_path),
        contract_value(contract, "molstar_reference_commit"),
        test_json_escape(contract_value(contract, "fixture")),
        input_format,
        representation,
        test_json_escape(assembly),
        sphere_detail,
        obj_path,
        stl_path,
        obj.len(),
        format!("{:016x}", stable_test_hash64(obj.as_bytes())),
        obj_stats.vertex_count,
        obj_stats.normal_count,
        obj_stats.face_count,
        pretty_f32_array_json(obj_stats.min),
        pretty_f32_array_json(obj_stats.max),
        stl_stats.byte_len,
        format!("{:016x}", stable_test_hash64(&stl)),
        stl_stats.facet_count,
        test_json_escape(stl_header),
        pretty_f32_array_json(stl_stats.min),
        pretty_f32_array_json(stl_stats.max)
    )
}

fn molstar_reference_input_format(contract: &str, options: &MeshOptions) -> String {
    if let Some(format) = contract_optional_value(contract, "input_format") {
        return format.to_string();
    }
    match options.format {
        InputFormat::Pdb => "pdb".to_string(),
        InputFormat::Cif => "cif".to_string(),
        InputFormat::BinaryCif => "bcif".to_string(),
        InputFormat::Auto => {
            let fixture = contract_value(contract, "fixture").to_ascii_lowercase();
            if fixture.ends_with(".bcif") {
                "bcif".to_string()
            } else if fixture.ends_with(".cif") || fixture.ends_with(".mmcif") {
                "cif".to_string()
            } else {
                "pdb".to_string()
            }
        }
    }
}

fn molstar_reference_representation(representation: Representation) -> &'static str {
    match representation {
        Representation::Molstar => "molstar",
        Representation::Spacefill => "spacefill",
        Representation::BallAndStick => "ball-and-stick",
        Representation::Cartoon => "cartoon",
        Representation::Ribbon => "ribbon",
        Representation::Backbone => "backbone",
    }
}

fn molstar_reference_sphere_detail(options_json: &str) -> usize {
    contract_json_number(options_json, "sphere-detail")
        .or_else(|| contract_json_number(options_json, "sphere_detail"))
        .unwrap_or(1.0) as usize
}

fn contract_json_number(text: &str, key: &str) -> Option<f32> {
    let quoted_key = format!("\"{key}\"");
    let key_start = text.find(&quoted_key)?;
    let after_key = &text[key_start + quoted_key.len()..];
    let colon = after_key.find(':')?;
    let value = after_key[colon + 1..].trim_start();
    let end = value
        .find(|ch: char| !(ch.is_ascii_digit() || matches!(ch, '.' | '-' | '+' | 'e' | 'E')))
        .unwrap_or(value.len());
    value[..end].parse::<f32>().ok()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct WasmMemorySummary {
    initial_pages: u32,
    maximum_pages: Option<u32>,
    exported_memory: bool,
}

fn parse_wasm_memory_summary(bytes: &[u8]) -> Result<WasmMemorySummary, String> {
    if bytes.len() < 8 || &bytes[..4] != b"\0asm" || bytes[4..8] != [1, 0, 0, 0] {
        return Err("not a WebAssembly MVP binary".to_string());
    }

    let mut offset = 8usize;
    let mut memory = None;
    let mut exported_memory = false;
    while offset < bytes.len() {
        let section_id = read_wasm_byte(bytes, &mut offset)?;
        let section_size = read_wasm_u32(bytes, &mut offset)? as usize;
        let section_end = offset
            .checked_add(section_size)
            .ok_or_else(|| "wasm section size overflow".to_string())?;
        if section_end > bytes.len() {
            return Err("wasm section extends past end of file".to_string());
        }

        match section_id {
            2 => {
                let count = read_wasm_u32(bytes, &mut offset)?;
                for _ in 0..count {
                    read_wasm_name(bytes, &mut offset)?;
                    read_wasm_name(bytes, &mut offset)?;
                    match read_wasm_byte(bytes, &mut offset)? {
                        0 => {
                            read_wasm_u32(bytes, &mut offset)?;
                        }
                        1 => {
                            read_wasm_byte(bytes, &mut offset)?;
                            read_wasm_limits(bytes, &mut offset)?;
                        }
                        2 => {
                            memory = Some(read_wasm_limits(bytes, &mut offset)?);
                        }
                        3 => {
                            read_wasm_byte(bytes, &mut offset)?;
                            read_wasm_byte(bytes, &mut offset)?;
                        }
                        kind => return Err(format!("unsupported wasm import kind {kind}")),
                    }
                }
            }
            5 => {
                let count = read_wasm_u32(bytes, &mut offset)?;
                for _ in 0..count {
                    memory = Some(read_wasm_limits(bytes, &mut offset)?);
                }
            }
            7 => {
                let count = read_wasm_u32(bytes, &mut offset)?;
                for _ in 0..count {
                    let name = read_wasm_name(bytes, &mut offset)?;
                    let kind = read_wasm_byte(bytes, &mut offset)?;
                    read_wasm_u32(bytes, &mut offset)?;
                    if kind == 2 && name == "memory" {
                        exported_memory = true;
                    }
                }
            }
            _ => {}
        }
        offset = section_end;
    }

    let (initial_pages, maximum_pages) =
        memory.ok_or_else(|| "wasm memory section not found".to_string())?;
    Ok(WasmMemorySummary {
        initial_pages,
        maximum_pages,
        exported_memory,
    })
}

fn read_wasm_byte(bytes: &[u8], offset: &mut usize) -> Result<u8, String> {
    let byte = *bytes
        .get(*offset)
        .ok_or_else(|| "unexpected end of wasm".to_string())?;
    *offset += 1;
    Ok(byte)
}

fn read_wasm_u32(bytes: &[u8], offset: &mut usize) -> Result<u32, String> {
    let mut result = 0u32;
    let mut shift = 0;
    loop {
        let byte = read_wasm_byte(bytes, offset)?;
        result |= u32::from(byte & 0x7f) << shift;
        if byte & 0x80 == 0 {
            return Ok(result);
        }
        shift += 7;
        if shift >= 35 {
            return Err("invalid wasm varuint32".to_string());
        }
    }
}

fn read_wasm_name(bytes: &[u8], offset: &mut usize) -> Result<String, String> {
    let len = read_wasm_u32(bytes, offset)? as usize;
    let end = (*offset)
        .checked_add(len)
        .ok_or_else(|| "wasm name length overflow".to_string())?;
    let name = std::str::from_utf8(
        bytes
            .get(*offset..end)
            .ok_or_else(|| "wasm name extends past end of file".to_string())?,
    )
    .map_err(|_| "wasm name is not UTF-8".to_string())?
    .to_string();
    *offset = end;
    Ok(name)
}

fn read_wasm_limits(bytes: &[u8], offset: &mut usize) -> Result<(u32, Option<u32>), String> {
    let flags = read_wasm_byte(bytes, offset)?;
    let min = read_wasm_u32(bytes, offset)?;
    let max = if flags & 0x01 != 0 {
        Some(read_wasm_u32(bytes, offset)?)
    } else {
        None
    };
    Ok((min, max))
}

fn grouped_triangle_mesh() -> Mesh {
    Mesh {
        vertices: vec![
            vec3(0.0, 0.0, 0.0),
            vec3(1.23456, 0.0, 0.0),
            vec3(1.0, 1.0, 0.0),
            vec3(0.0, 1.0, 0.0),
        ],
        normals: vec![
            vec3(0.0, 0.0, 1.0),
            vec3(0.0, 0.0, 1.0),
            vec3(0.0, 0.0, 1.0),
            vec3(0.0, 0.0, 1.0),
        ],
        faces: vec![Face { a: 0, b: 1, c: 2 }, Face { a: 0, b: 2, c: 3 }],
        vertex_groups: vec![0, 0, 1, 1],
        face_groups: vec![0, 1],
        face_materials: Vec::new(),
        sections: Vec::new(),
        group_count: 2,
    }
}

fn sheet_test_molecule() -> Molecule {
    Molecule {
        atoms: [
            vec3(0.0, 0.0, 0.0),
            vec3(1.1, 0.25, 0.2),
            vec3(2.2, -0.2, 0.1),
            vec3(3.3, 0.15, 0.0),
        ]
        .into_iter()
        .enumerate()
        .map(|(index, position)| test_atom(index + 1, "CA", "A", index as i32 + 1, position))
        .collect(),
        sheets: vec![SecondaryRange {
            chain: "A".to_string(),
            start: 1,
            start_insertion_code: String::new(),
            end: 4,
            end_insertion_code: String::new(),
        }],
        ..Molecule::default()
    }
}

fn nucleotide_molecule(residue: &str, atoms: &[(&str, Vec3)]) -> Molecule {
    Molecule {
        atoms: atoms
            .iter()
            .enumerate()
            .map(|(index, (name, position))| {
                let mut atom = test_atom(index + 1, name, "N", 1, *position);
                atom.residue = residue.to_string();
                atom.element = name
                    .chars()
                    .find(|ch| ch.is_ascii_alphabetic())
                    .map(|ch| ch.to_string())
                    .unwrap_or_else(|| "C".to_string());
                atom
            })
            .collect(),
        ..Molecule::default()
    }
}

struct AssemblyReferenceCase {
    case_name: &'static str,
    fixture: &'static str,
    input: &'static str,
    data: &'static [u8],
    format: InputFormat,
    format_name: &'static str,
    assembly_id: &'static str,
    alt_loc: &'static str,
}

fn assembly_fixture_reference_summary_json() -> String {
    let cases = [
        AssemblyReferenceCase {
            case_name: "assembly-altloc-helix-cif",
            fixture: "assembly-altloc-helix",
            input: "tests/fixtures/cif/assembly-altloc-helix.cif",
            data: include_bytes!("../../tests/fixtures/cif/assembly-altloc-helix.cif"),
            format: InputFormat::Cif,
            format_name: "cif",
            assembly_id: "1",
            alt_loc: "A",
        },
        AssemblyReferenceCase {
            case_name: "assembly-altloc-helix-bcif",
            fixture: "assembly-altloc-helix",
            input: "tests/fixtures/bcif/assembly-altloc-helix.bcif",
            data: include_bytes!("../../tests/fixtures/bcif/assembly-altloc-helix.bcif"),
            format: InputFormat::BinaryCif,
            format_name: "bcif",
            assembly_id: "1",
            alt_loc: "A",
        },
        AssemblyReferenceCase {
            case_name: "assembly-altloc-secondary-cif",
            fixture: "assembly-altloc-secondary",
            input: "tests/fixtures/cif/assembly-altloc-secondary.cif",
            data: include_bytes!("../../tests/fixtures/cif/assembly-altloc-secondary.cif"),
            format: InputFormat::Cif,
            format_name: "cif",
            assembly_id: "1",
            alt_loc: "A",
        },
        AssemblyReferenceCase {
            case_name: "assembly-altloc-secondary-bcif",
            fixture: "assembly-altloc-secondary",
            input: "tests/fixtures/bcif/assembly-altloc-secondary.bcif",
            data: include_bytes!("../../tests/fixtures/bcif/assembly-altloc-secondary.bcif"),
            format: InputFormat::BinaryCif,
            format_name: "bcif",
            assembly_id: "1",
            alt_loc: "A",
        },
        AssemblyReferenceCase {
            case_name: "assembly-single-operator-cif",
            fixture: "assembly-operator-matrix",
            input: "tests/fixtures/cif/assembly-operator-matrix.cif",
            data: include_bytes!("../../tests/fixtures/cif/assembly-operator-matrix.cif"),
            format: InputFormat::Cif,
            format_name: "cif",
            assembly_id: "single",
            alt_loc: "",
        },
        AssemblyReferenceCase {
            case_name: "assembly-single-operator-bcif",
            fixture: "assembly-operator-matrix",
            input: "tests/fixtures/bcif/assembly-operator-matrix.bcif",
            data: include_bytes!("../../tests/fixtures/bcif/assembly-operator-matrix.bcif"),
            format: InputFormat::BinaryCif,
            format_name: "bcif",
            assembly_id: "single",
            alt_loc: "",
        },
        AssemblyReferenceCase {
            case_name: "assembly-operator-range-cif",
            fixture: "assembly-operator-matrix",
            input: "tests/fixtures/cif/assembly-operator-matrix.cif",
            data: include_bytes!("../../tests/fixtures/cif/assembly-operator-matrix.cif"),
            format: InputFormat::Cif,
            format_name: "cif",
            assembly_id: "range",
            alt_loc: "",
        },
        AssemblyReferenceCase {
            case_name: "assembly-operator-range-bcif",
            fixture: "assembly-operator-matrix",
            input: "tests/fixtures/bcif/assembly-operator-matrix.bcif",
            data: include_bytes!("../../tests/fixtures/bcif/assembly-operator-matrix.bcif"),
            format: InputFormat::BinaryCif,
            format_name: "bcif",
            assembly_id: "range",
            alt_loc: "",
        },
        AssemblyReferenceCase {
            case_name: "assembly-cartesian-product-cif",
            fixture: "assembly-operator-matrix",
            input: "tests/fixtures/cif/assembly-operator-matrix.cif",
            data: include_bytes!("../../tests/fixtures/cif/assembly-operator-matrix.cif"),
            format: InputFormat::Cif,
            format_name: "cif",
            assembly_id: "cart",
            alt_loc: "",
        },
        AssemblyReferenceCase {
            case_name: "assembly-cartesian-product-bcif",
            fixture: "assembly-operator-matrix",
            input: "tests/fixtures/bcif/assembly-operator-matrix.bcif",
            data: include_bytes!("../../tests/fixtures/bcif/assembly-operator-matrix.bcif"),
            format: InputFormat::BinaryCif,
            format_name: "bcif",
            assembly_id: "cart",
            alt_loc: "",
        },
    ];
    format!(
        "{{\"molstar_reference_commit\":\"{}\",\"molstar_sources\":[\"artifacts/molstar/src/mol-model-formats/structure/property/assembly.ts\",\"artifacts/molstar/src/mol-model/structure/structure/symmetry.ts\",\"artifacts/molstar/src/mol-model/structure/structure/structure.ts\",\"artifacts/molstar/src/mol-model/structure/structure/unit.ts\",\"artifacts/molstar/src/mol-math/geometry/symmetry-operator.ts\"],\"summary_name\":\"assembly-fixture-reference-comparison\",\"source_fields\":[\"createAssemblies\",\"operatorGroupsProvider\",\"parseOperatorList\",\"expandOperators\",\"getAssemblyOperators\",\"StructureSymmetry.buildAssembly\",\"Structure.Builder.addWithOperator\",\"Unit.SymmetryGroup\"],\"cases\":[{}]}}",
        MOLSTAR_REFERENCE_COMMIT,
        cases
            .iter()
            .map(assembly_reference_case_json)
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn assembly_reference_case_json(case: &AssemblyReferenceCase) -> String {
    let molecule = parse_molecule_with_options(
        case.data,
        &MeshOptions {
            format: case.format,
            assembly: Some(case.assembly_id.to_string()),
            alt_loc: case.alt_loc.to_string(),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let assembly = molecule.selected_assembly.as_ref().unwrap();
    let structure = molecule.atomic_structure();
    let (expanded_min, expanded_max) =
        expanded_atom_bounds(&molecule.expanded_for_geometry()).unwrap_or_default();
    let operator_count = assembly
        .generators
        .iter()
        .map(|generator| generator.operators.len())
        .sum::<usize>();

    format!(
        "{{\"case\":\"{}\",\"fixture\":\"{}\",\"input\":\"{}\",\"format\":\"{}\",\"assembly_id\":\"{}\",\"alt_loc\":\"{}\",\"assembly\":{{\"asym_ids\":{},\"generator_count\":{},\"operator_count\":{},\"operators\":[{}]}},\"structure\":{{\"model_count\":{},\"unit_count\":{},\"element_count\":{},\"symmetry_group_count\":{},\"intra_unit_bond_count\":{},\"inter_unit_bond_count\":{},\"coordinate_system\":{},\"boundary\":{},\"expanded_atom_bounds\":{{\"min\":{},\"max\":{}}}}},\"units\":[{}],\"symmetry_groups\":[{}]}}",
        test_json_escape(case.case_name),
        test_json_escape(case.fixture),
        test_json_escape(case.input),
        test_json_escape(case.format_name),
        test_json_escape(case.assembly_id),
        test_json_escape(case.alt_loc),
        crate::json::json_string_array(&assembly.asym_ids),
        assembly.generators.len(),
        operator_count,
        assembly_operator_reference_json(assembly),
        structure.models.len(),
        structure.units.len(),
        structure.element_count,
        structure.symmetry_groups.len(),
        structure.intra_unit_bond_count,
        structure.inter_unit_bonds.len(),
        unit_operator_json(&structure.coordinate_system),
        boundary_json(&structure.boundary),
        test_vec3_json(expanded_min),
        test_vec3_json(expanded_max),
        structure
            .units
            .iter()
            .map(assembly_unit_reference_json)
            .collect::<Vec<_>>()
            .join(","),
        structure
            .symmetry_groups
            .iter()
            .map(assembly_symmetry_group_reference_json)
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn assembly_operator_reference_json(assembly: &Assembly) -> String {
    assembly
        .generators
        .iter()
        .flat_map(|generator| generator.operators.iter())
        .map(|operator| {
            let unit_operator = UnitOperator {
                name: operator.name.clone(),
                instance_id: operator.instance_id.clone(),
                assembly_id: operator.assembly_id.clone(),
                oper_id: operator.oper_id as i32,
                oper_list_ids: operator.oper_list_ids.clone(),
                transform: operator.transform,
                is_identity: operator.transform.is_identity(),
                suffix: if operator.transform.is_identity() {
                    String::new()
                } else {
                    format!("_{}", operator.oper_id)
                },
            };
            format!(
                "{{\"operator\":{},\"probe\":{}}}",
                unit_operator_json(&unit_operator),
                test_vec3_json(operator.transform.apply(vec3(1.0, 2.0, 3.0)))
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn assembly_unit_reference_json(unit: &crate::model::StructureUnit) -> String {
    format!(
        "{{\"id\":{},\"kind\":\"{}\",\"invariant_id\":{},\"chain_group_id\":{},\"chain_indices\":{},\"elements\":{},\"atom_indices\":{},\"residue_indices\":{},\"operator\":{},\"boundary\":{},\"intra_unit_bond_count\":{},\"inter_unit_bond_count\":{}}}",
        unit.id,
        unit_kind_json(unit.kind),
        unit.invariant_id,
        unit.chain_group_id,
        usize_array_json(&unit.chain_indices),
        usize_array_json(&unit.elements),
        usize_array_json(&unit.atom_indices),
        usize_array_json(&unit.residue_indices),
        unit_operator_json(&unit.operator),
        boundary_json(&unit.props.boundary),
        unit.props.intra_unit_bond_count,
        unit.props.inter_unit_bond_count
    )
}

fn assembly_symmetry_group_reference_json(group: &crate::model::UnitSymmetryGroup) -> String {
    format!(
        "{{\"kind\":\"{}\",\"model_id\":{},\"invariant_id\":{},\"elements\":{},\"unit_ids\":{},\"operator_names\":{},\"operator_instance_ids\":{},\"hash_code\":{},\"transform_hash\":{}}}",
        unit_kind_json(group.kind),
        group.model_id,
        group.invariant_id,
        usize_array_json(&group.elements),
        usize_array_json(&group.unit_ids),
        crate::json::json_string_array(&group.operator_names),
        crate::json::json_string_array(&group.operator_instance_ids),
        group.hash_code,
        group.transform_hash
    )
}

fn unit_operator_json(operator: &UnitOperator) -> String {
    format!(
        "{{\"name\":\"{}\",\"instance_id\":\"{}\",\"assembly_id\":\"{}\",\"oper_id\":{},\"oper_list_ids\":{},\"is_identity\":{},\"suffix\":\"{}\",\"transform\":{}}}",
        test_json_escape(&operator.name),
        test_json_escape(&operator.instance_id),
        test_json_escape(&operator.assembly_id),
        operator.oper_id,
        crate::json::json_string_array(&operator.oper_list_ids),
        if operator.is_identity { "true" } else { "false" },
        test_json_escape(&operator.suffix),
        transform_json(operator.transform)
    )
}

fn transform_json(transform: Transform) -> String {
    format!(
        "[{},{},{}]",
        f32_array_json(&transform.m[0]),
        f32_array_json(&transform.m[1]),
        f32_array_json(&transform.m[2])
    )
}

fn boundary_json(boundary: &crate::model::Boundary) -> String {
    format!(
        "{{\"box_min\":{},\"box_max\":{},\"sphere_center\":{},\"sphere_radius\":{}}}",
        test_vec3_json(boundary.box_min),
        test_vec3_json(boundary.box_max),
        test_vec3_json(boundary.sphere.center),
        test_f32_json(boundary.sphere.radius)
    )
}

fn expanded_atom_bounds(molecule: &Molecule) -> Option<(Vec3, Vec3)> {
    let mut atoms = molecule.atoms.iter();
    let first = atoms.next()?.position;
    let mut min = first;
    let mut max = first;
    for atom in atoms {
        min = min.min(atom.position);
        max = max.max(atom.position);
    }
    Some((min, max))
}

fn unit_kind_json(kind: UnitKind) -> &'static str {
    match kind {
        UnitKind::Atomic => "atomic",
        UnitKind::Spheres => "spheres",
        UnitKind::Gaussians => "gaussians",
    }
}

fn usize_array_json(values: &[usize]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn f32_array_json(values: &[f32]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|value| test_f32_json(*value))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn nucleotide_reference_summary_json(molecule: &Molecule) -> String {
    let structure = molecule.atomic_structure();
    let hierarchy = &structure.model.hierarchy;
    let options = MeshOptions {
        representation: Representation::Molstar,
        center: false,
        assembly: None,
        infer_bonds: false,
        ..MeshOptions::default()
    };
    let rings = build_render_objects(molecule, &options)
        .into_iter()
        .filter_map(|object| match object {
            RenderObject::NucleotideRing {
                center,
                normal,
                base: Some(base),
                ..
            } => Some((center, normal, base)),
            _ => None,
        })
        .collect::<Vec<_>>();
    let residues = hierarchy
        .residues
        .iter()
        .enumerate()
        .filter(|(residue_index, _)| {
            matches!(
                hierarchy.derived.residue.molecule_type.get(*residue_index),
                Some(MoleculeType::Rna | MoleculeType::Dna | MoleculeType::Pna)
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(rings.len(), residues.len());
    let base_summaries = residues
        .iter()
        .zip(rings.iter())
        .map(|(&(residue_index, residue), &(center, normal, base))| {
            let chain = &hierarchy.chains[residue.chain_index];
            let molecule_type = hierarchy.derived.residue.molecule_type[residue_index];
            let polymer_type = hierarchy.derived.residue.polymer_type[residue_index];
            let trace_element = hierarchy.derived.residue.trace_element_index[residue_index]
                .expect("nucleotide trace element");
            let (base_kind, molstar_index_names, positions_json) = match base {
                crate::mesh::NucleotideRingBase::PurineConnector { trace, n9 } => (
                    "purine-connector",
                    vec!["trace", "N9"],
                    nucleotide_position_json(&[("trace", trace), ("N9", n9)]),
                ),
                crate::mesh::NucleotideRingBase::Purine {
                    trace,
                    n1,
                    c2,
                    n3,
                    c4,
                    c5,
                    c6,
                    n7,
                    c8,
                    n9,
                } => (
                    "purine",
                    vec!["trace", "N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"],
                    nucleotide_position_json(&[
                        ("trace", trace),
                        ("N1", n1),
                        ("C2", c2),
                        ("N3", n3),
                        ("C4", c4),
                        ("C5", c5),
                        ("C6", c6),
                        ("N7", n7),
                        ("C8", c8),
                        ("N9", n9),
                    ]),
                ),
                crate::mesh::NucleotideRingBase::Pyrimidine {
                    trace,
                    n1,
                    c2,
                    n3,
                    c4,
                    c5,
                    c6,
                } => (
                    "pyrimidine",
                    vec!["trace", "N1", "C2", "N3", "C4", "C5", "C6"],
                    nucleotide_position_json(&[
                        ("trace", trace),
                        ("N1", n1),
                        ("C2", c2),
                        ("N3", n3),
                        ("C4", c4),
                        ("C5", c5),
                        ("C6", c6),
                    ]),
                ),
                crate::mesh::NucleotideRingBase::PyrimidineConnector { trace, n1 } => (
                    "pyrimidine-connector",
                    vec!["trace", "N1"],
                    nucleotide_position_json(&[("trace", trace), ("N1", n1)]),
                ),
            };
            let molstar_indices = nucleotide_atom_refs_json(
                hierarchy,
                residue,
                trace_element,
                &molstar_index_names,
            );
            format!(
                "{{\"residue_index\":{},\"chain\":\"{}\",\"residue\":\"{}\",\"seq\":\"{}\",\"molecule_type\":\"{}\",\"polymer_type\":\"{}\",\"base_kind\":\"{}\",\"trace_element\":{},\"center\":{},\"normal\":{},\"molstar_indices\":{},\"positions\":{}}}",
                residue_index,
                test_json_escape(&chain.id),
                test_json_escape(&residue.comp_id),
                test_json_escape(&residue.label_seq_id),
                molecule_type_json(molecule_type),
                polymer_type_json(polymer_type),
                base_kind,
                trace_element,
                test_vec3_json(center),
                test_vec3_json(normal),
                molstar_indices,
                positions_json
            )
        })
        .collect::<Vec<_>>()
        .join(",");

    format!(
        "{{\"molstar_reference_commit\":\"{}\",\"molstar_sources\":[\"artifacts/molstar/src/mol-repr/structure/visual/util/nucleotide.ts\",\"artifacts/molstar/src/mol-repr/structure/visual/nucleotide-ring-mesh.ts\",\"artifacts/molstar/src/mol-repr/structure/visual/nucleotide-block-mesh.ts\",\"artifacts/molstar/src/mol-repr/structure/visual/polymer-direction-wedge.ts\"],\"fixture\":\"tests/fixtures/cif/nucleic-acid-rna-dna.cif\",\"binary_fixture\":\"tests/fixtures/bcif/nucleic-acid-rna-dna.bcif\",\"summary_name\":\"nucleotide-rna-dna-reference\",\"bases\":[{}]}}",
        MOLSTAR_REFERENCE_COMMIT,
        base_summaries
    )
}

fn nucleotide_atom_refs_json(
    hierarchy: &crate::model::AtomicHierarchy,
    residue: &crate::model::AtomicResidue,
    trace_element: usize,
    names: &[&str],
) -> String {
    format!(
        "[{}]",
        names
            .iter()
            .map(|name| {
                let element = if *name == "trace" {
                    trace_element
                } else {
                    find_atom_on_residue(hierarchy, residue, name)
                        .unwrap_or_else(|| panic!("missing {name} on nucleotide residue"))
                };
                let atom = &hierarchy.atoms[element];
                format!(
                    "{{\"name\":\"{}\",\"element\":{},\"source_atom_index\":{},\"atom_id\":{},\"atom_name\":\"{}\"}}",
                    test_json_escape(name),
                    element,
                    atom.source_index,
                    atom.id,
                    test_json_escape(&atom.name)
                )
            })
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn find_atom_on_residue(
    hierarchy: &crate::model::AtomicHierarchy,
    residue: &crate::model::AtomicResidue,
    name: &str,
) -> Option<usize> {
    (residue.start_atom..residue.end_atom).find(|&element| hierarchy.atoms[element].name == name)
}

fn nucleotide_position_json(values: &[(&str, Vec3)]) -> String {
    format!(
        "{{{}}}",
        values
            .iter()
            .map(|(name, position)| {
                format!(
                    "\"{}\":{}",
                    test_json_escape(name),
                    test_vec3_json(*position)
                )
            })
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn molecule_type_json(value: MoleculeType) -> &'static str {
    match value {
        MoleculeType::Unknown => "unknown",
        MoleculeType::Other => "other",
        MoleculeType::Water => "water",
        MoleculeType::Ion => "ion",
        MoleculeType::Lipid => "lipid",
        MoleculeType::Protein => "protein",
        MoleculeType::Rna => "rna",
        MoleculeType::Dna => "dna",
        MoleculeType::Pna => "pna",
        MoleculeType::Saccharide => "saccharide",
    }
}

fn polymer_type_json(value: crate::model::PolymerType) -> &'static str {
    match value {
        crate::model::PolymerType::None => "none",
        crate::model::PolymerType::PeptideL => "peptide-l",
        crate::model::PolymerType::GammaPeptide => "gamma-peptide",
        crate::model::PolymerType::BetaPeptide => "beta-peptide",
        crate::model::PolymerType::Rna => "rna",
        crate::model::PolymerType::Dna => "dna",
        crate::model::PolymerType::Pna => "pna",
    }
}

struct BondRingReferenceCase {
    case_name: &'static str,
    input: &'static str,
    data: &'static [u8],
    format: InputFormat,
    format_name: &'static str,
    assembly_id: Option<&'static str>,
}

fn bond_ring_reference_summary_json() -> String {
    let cases = [
        BondRingReferenceCase {
            case_name: "ligand-metal-aromatic-cif",
            input: "tests/fixtures/cif/ligand-metal-aromatic.cif",
            data: include_bytes!("../../tests/fixtures/cif/ligand-metal-aromatic.cif"),
            format: InputFormat::Cif,
            format_name: "cif",
            assembly_id: None,
        },
        BondRingReferenceCase {
            case_name: "ligand-metal-aromatic-bcif",
            input: "tests/fixtures/bcif/ligand-metal-aromatic.bcif",
            data: include_bytes!("../../tests/fixtures/bcif/ligand-metal-aromatic.bcif"),
            format: InputFormat::BinaryCif,
            format_name: "bcif",
            assembly_id: None,
        },
        BondRingReferenceCase {
            case_name: "carbohydrate-branched-cif",
            input: "tests/fixtures/cif/carbohydrate-branched.cif",
            data: include_bytes!("../../tests/fixtures/cif/carbohydrate-branched.cif"),
            format: InputFormat::Cif,
            format_name: "cif",
            assembly_id: None,
        },
        BondRingReferenceCase {
            case_name: "covalent-cross-link-cif",
            input: "tests/fixtures/cif/covalent-cross-link.cif",
            data: include_bytes!("../../tests/fixtures/cif/covalent-cross-link.cif"),
            format: InputFormat::Cif,
            format_name: "cif",
            assembly_id: None,
        },
    ];
    format!(
        "{{\"molstar_reference_commit\":\"{}\",\"molstar_sources\":[\"artifacts/molstar/src/mol-model/structure/structure/unit/bonds/intra-compute.ts\",\"artifacts/molstar/src/mol-model/structure/structure/unit/bonds/inter-compute.ts\",\"artifacts/molstar/src/mol-model/structure/structure/unit/bonds/data.ts\",\"artifacts/molstar/src/mol-model-formats/structure/property/bonds/index-pair.ts\",\"artifacts/molstar/src/mol-model/structure/structure/unit/rings/compute.ts\",\"artifacts/molstar/src/mol-model/structure/structure/unit/rings.ts\",\"artifacts/molstar/src/mol-math/graph/int-adjacency-graph.ts\"],\"summary_name\":\"bond-ring-graph-reference-comparison\",\"cases\":[{}]}}",
        MOLSTAR_REFERENCE_COMMIT,
        cases
            .iter()
            .map(bond_ring_reference_case_json)
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn bond_ring_reference_case_json(case: &BondRingReferenceCase) -> String {
    let molecule = parse_molecule_with_options(
        case.data,
        &MeshOptions {
            format: case.format,
            assembly: case.assembly_id.map(str::to_string),
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();
    let structure = molecule.atomic_structure();
    format!(
        "{{\"case\":\"{}\",\"input\":\"{}\",\"format\":\"{}\",\"assembly_id\":{},\"molecule\":{{\"atom_count\":{},\"bond_count\":{},\"bonds\":[{}],\"index_pair\":{},\"rings\":[{}],\"resonance\":{}}},\"structure\":{{\"unit_count\":{},\"intra_unit_bond_count\":{},\"inter_unit_bond_count\":{},\"units\":[{}],\"inter_unit_bonds\":[{}],\"inter_unit_graph\":{}}}}}",
        test_json_escape(case.case_name),
        test_json_escape(case.input),
        test_json_escape(case.format_name),
        opt_str_json(case.assembly_id),
        molecule.atoms.len(),
        molecule.bonds.len(),
        molecule_bonds_json(&molecule),
        index_pair_bonds_json(molecule.index_pair_bonds.as_ref()),
        molecule_rings_json(&molecule),
        resonance_json(&molecule.resonance),
        structure.units.len(),
        structure.intra_unit_bond_count,
        structure.inter_unit_bonds.len(),
        structure
            .units
            .iter()
            .map(structure_unit_bond_json)
            .collect::<Vec<_>>()
            .join(","),
        inter_unit_bonds_json(&structure.inter_unit_bonds),
        inter_unit_bond_graph_json(&structure.inter_unit_bond_graph)
    )
}

fn molecule_bonds_json(molecule: &Molecule) -> String {
    molecule
        .bonds
        .iter()
        .enumerate()
        .map(|(index, bond)| {
            let metadata = molecule
                .bond_metadata
                .get(index)
                .cloned()
                .unwrap_or_default();
            let atom_a = &molecule.atoms[bond.a];
            let atom_b = &molecule.atoms[bond.b];
            format!(
                "{{\"index\":{},\"a\":{},\"b\":{},\"atom_ids\":[{},{}],\"atom_names\":[\"{}\",\"{}\"],\"source\":\"{}\",\"order\":{},\"flags\":{},\"key\":{},\"distance\":{},\"operator_a\":{},\"operator_b\":{},\"struct_conn\":{}}}",
                index,
                bond.a,
                bond.b,
                atom_a.id,
                atom_b.id,
                test_json_escape(&atom_a.name),
                test_json_escape(&atom_b.name),
                bond_source_json(&metadata.source),
                metadata.order,
                metadata.flags.bits,
                metadata.key,
                opt_f32_json(metadata.distance),
                metadata.operator_a,
                metadata.operator_b,
                metadata
                    .struct_conn
                    .as_ref()
                    .map(struct_conn_json)
                    .unwrap_or_else(|| "null".to_string())
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn struct_conn_json(struct_conn: &StructConnMetadata) -> String {
    format!(
        "{{\"id\":\"{}\",\"row_index\":{},\"partner_a_atom_index\":{},\"partner_b_atom_index\":{},\"conn_type_id\":\"{}\",\"value_order\":\"{}\",\"partner_a_symmetry\":\"{}\",\"partner_b_symmetry\":\"{}\"}}",
        test_json_escape(&struct_conn.id),
        struct_conn.row_index,
        struct_conn.partner_a_atom_index,
        struct_conn.partner_b_atom_index,
        test_json_escape(&struct_conn.conn_type_id),
        test_json_escape(&struct_conn.value_order),
        test_json_escape(&struct_conn.partner_a_symmetry),
        test_json_escape(&struct_conn.partner_b_symmetry)
    )
}

fn index_pair_bonds_json(index_pairs: Option<&IndexPairBonds>) -> String {
    let Some(index_pairs) = index_pairs else {
        return "null".to_string();
    };
    format!(
        "{{\"max_distance\":{},\"cacheable\":{},\"has_operators\":{},\"by_same_operator\":[{}],\"graph\":{}}}",
        test_f32_json(index_pairs.max_distance),
        index_pairs.cacheable,
        index_pairs.has_operators,
        index_pairs
            .by_same_operator
            .iter()
            .map(|(operator, slots)| format!(
                "{{\"operator\":{},\"slots\":{}}}",
                operator,
                usize_array_json(slots)
            ))
            .collect::<Vec<_>>()
            .join(","),
        index_pair_graph_json(&index_pairs.bonds)
    )
}

fn index_pair_graph_json(graph: &crate::model::IndexPairGraph) -> String {
    format!(
        "{{\"vertex_count\":{},\"offset\":{},\"a\":{},\"b\":{},\"edge_count\":{},\"key\":{},\"operator_a\":{},\"operator_b\":{},\"order\":{},\"distance\":{},\"flag\":{}}}",
        graph.vertex_count,
        usize_array_json(&graph.offset),
        usize_array_json(&graph.a),
        usize_array_json(&graph.b),
        graph.edge_count,
        i32_array_json(&graph.props.key),
        i32_array_json(&graph.props.operator_a),
        i32_array_json(&graph.props.operator_b),
        i8_array_json(&graph.props.order),
        f32_array_json(&graph.props.distance),
        flags_array_json(&graph.props.flag)
    )
}

fn structure_unit_bond_json(unit: &crate::model::StructureUnit) -> String {
    format!(
        "{{\"id\":{},\"kind\":\"{}\",\"elements\":{},\"atom_indices\":{},\"intra_unit_bond_count\":{},\"inter_unit_bond_count\":{},\"intra_unit_bonds\":{}}}",
        unit.id,
        unit_kind_json(unit.kind),
        usize_array_json(&unit.elements),
        usize_array_json(&unit.atom_indices),
        unit.props.intra_unit_bond_count,
        unit.props.inter_unit_bond_count,
        intra_unit_bonds_json(&unit.props.intra_unit_bonds)
    )
}

fn intra_unit_bonds_json(bonds: &crate::model::IntraUnitBonds) -> String {
    format!(
        "{{\"vertex_count\":{},\"offset\":{},\"a\":{},\"b\":{},\"edge_count\":{},\"key\":{},\"order\":{},\"flags\":{},\"can_remap\":{},\"cacheable\":{}}}",
        bonds.vertex_count,
        usize_array_json(&bonds.offset),
        usize_array_json(&bonds.a),
        usize_array_json(&bonds.b),
        bonds.edge_count,
        i32_array_json(&bonds.props.key),
        i8_array_json(&bonds.props.order),
        flags_array_json(&bonds.props.flags),
        bonds.can_remap,
        bonds.cacheable
    )
}

fn inter_unit_bonds_json(bonds: &[crate::model::InterUnitBond]) -> String {
    bonds
        .iter()
        .map(|bond| {
            format!(
                "{{\"unit_a\":{},\"index_a\":{},\"unit_b\":{},\"index_b\":{},\"source_bond\":{},\"order\":{},\"flags\":{},\"key\":{}}}",
                bond.unit_a,
                bond.index_a,
                bond.unit_b,
                bond.index_b,
                bond.source_bond,
                bond.order,
                bond.flags.bits,
                bond.key
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn inter_unit_bond_graph_json(graph: &crate::model::InterUnitBonds) -> String {
    format!(
        "{{\"edge_count\":{},\"edges\":[{}]}}",
        graph.edge_count,
        graph
            .edges
            .iter()
            .map(|edge| {
                format!(
                    "{{\"unit_a\":{},\"unit_b\":{},\"index_a\":{},\"index_b\":{},\"order\":{},\"flag\":{},\"key\":{}}}",
                    edge.unit_a,
                    edge.unit_b,
                    edge.index_a,
                    edge.index_b,
                    edge.props.order,
                    edge.props.flag.bits,
                    edge.props.key
                )
            })
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn molecule_rings_json(molecule: &Molecule) -> String {
    molecule
        .rings
        .iter()
        .map(|ring| {
            let atom_ids = ring
                .atom_indices
                .iter()
                .filter_map(|&index| molecule.atoms.get(index).map(|atom| atom.id))
                .collect::<Vec<_>>();
            format!(
                "{{\"atom_indices\":{},\"atom_ids\":{},\"aromatic\":{},\"fingerprint\":\"{}\"}}",
                usize_array_json(&ring.atom_indices),
                usize_array_json(&atom_ids),
                ring.aromatic,
                test_json_escape(&ring.fingerprint)
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn resonance_json(resonance: &crate::model::Resonance) -> String {
    format!(
        "{{\"ring_count\":{},\"aromatic_ring_count\":{},\"delocalized_bond_count\":{},\"delocalized_triplets\":[{}],\"element_ring_indices\":{},\"element_aromatic_ring_indices\":{},\"ring_component_index\":{},\"ring_components\":{}}}",
        resonance.ring_count,
        resonance.aromatic_ring_count,
        resonance.delocalized_bond_count,
        triplets_json(&resonance.delocalized_triplets),
        nested_usize_array_json(&resonance.element_ring_indices),
        nested_usize_array_json(&resonance.element_aromatic_ring_indices),
        usize_array_json(&resonance.ring_component_index),
        nested_usize_array_json(&resonance.ring_components)
    )
}

fn triplets_json(triplets: &[[usize; 3]]) -> String {
    triplets
        .iter()
        .map(|triplet| format!("[{},{},{}]", triplet[0], triplet[1], triplet[2]))
        .collect::<Vec<_>>()
        .join(",")
}

fn nested_usize_array_json(values: &[Vec<usize>]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|value| usize_array_json(value))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn flags_array_json(values: &[BondFlags]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|flags| flags.bits.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn i32_array_json(values: &[i32]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn i8_array_json(values: &[i8]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn opt_f32_json(value: Option<f32>) -> String {
    value
        .map(test_f32_json)
        .unwrap_or_else(|| "null".to_string())
}

fn opt_str_json(value: Option<&str>) -> String {
    value
        .map(|value| format!("\"{}\"", test_json_escape(value)))
        .unwrap_or_else(|| "null".to_string())
}

fn bond_source_json(source: &BondSource) -> &'static str {
    match source {
        BondSource::Computed => "computed",
        BondSource::PdbConect => "pdb-conect",
        BondSource::StructConn => "struct-conn",
        BondSource::IndexPair => "index-pair",
        BondSource::ChemComp => "chem-comp",
    }
}

fn test_vec3_json(value: Vec3) -> String {
    format!(
        "[{},{},{}]",
        test_f32_json(value.x),
        test_f32_json(value.y),
        test_f32_json(value.z)
    )
}

fn test_f32_json(value: f32) -> String {
    let value = if value.abs() < 0.000_05 { 0.0 } else { value };
    if (value.round() - value).abs() < 0.000_05 {
        format!("{value:.0}")
    } else {
        let mut out = format!("{value:.4}");
        while out.contains('.') && out.ends_with('0') {
            out.pop();
        }
        if out.ends_with('.') {
            out.pop();
        }
        out
    }
}

fn topology_test_molecule(
    residue: &str,
    atom_names: &[&str],
    bond_pairs: &[(usize, usize)],
) -> Molecule {
    let atoms = atom_names
        .iter()
        .enumerate()
        .map(|(index, name)| {
            let mut atom = test_atom(index + 1, name, "A", 1, vec3(index as f32, 0.0, 0.0));
            atom.residue = residue.to_string();
            atom.auth_residue = residue.to_string();
            atom
        })
        .collect::<Vec<_>>();
    let mut molecule = Molecule {
        atoms,
        bonds: bond_pairs
            .iter()
            .map(|(a, b)| Bond { a: *a, b: *b })
            .collect(),
        bond_metadata: vec![
            BondMetadata {
                source: BondSource::IndexPair,
                order: 1,
                flags: BondFlags::COVALENT,
                key: -1,
                distance: None,
                operator_a: -1,
                operator_b: -1,
                struct_conn: None,
            };
            bond_pairs.len()
        ],
        ..Molecule::default()
    };
    molecule.refresh_topology_metadata();
    molecule
}

fn test_atom(id: usize, name: &str, chain: &str, seq: i32, position: Vec3) -> Atom {
    Atom {
        id,
        source_index: id.saturating_sub(1),
        model_num: 1,
        name: name.to_string(),
        type_symbol: name
            .chars()
            .find(|ch| ch.is_ascii_alphabetic())
            .map(|ch| ch.to_ascii_uppercase().to_string())
            .unwrap_or_else(|| "C".to_string()),
        element: name
            .chars()
            .find(|ch| ch.is_ascii_alphabetic())
            .map(|ch| ch.to_string())
            .unwrap_or_else(|| "C".to_string()),
        chain: chain.to_string(),
        auth_chain: chain.to_string(),
        entity_id: String::new(),
        residue: "ALA".to_string(),
        auth_residue: "ALA".to_string(),
        group_pdb: "ATOM".to_string(),
        residue_seq: seq.to_string(),
        auth_residue_seq: seq.to_string(),
        insertion_code: String::new(),
        alt_id: String::new(),
        auth_name: name.to_string(),
        occupancy: 1.0,
        b_iso: 0.0,
        formal_charge: 0,
        position,
        het: false,
        operator_name: String::new(),
    }
}

fn het_atom(id: usize, name: &str, chain: &str, seq: i32, residue: &str, position: Vec3) -> Atom {
    let mut atom = test_atom(id, name, chain, seq, position);
    atom.residue = residue.to_string();
    atom.auth_residue = residue.to_string();
    atom.group_pdb = "HETATM".to_string();
    atom.het = true;
    atom
}

fn carbohydrate_atom(
    id: usize,
    name: &str,
    chain: &str,
    seq: i32,
    residue: &str,
    position: Vec3,
) -> Atom {
    let mut atom = test_atom(id, name, chain, seq, position);
    atom.residue = residue.to_string();
    atom.auth_residue = residue.to_string();
    atom.group_pdb = "HETATM".to_string();
    atom.het = true;
    atom
}

fn carbohydrate_bonds(bond_pairs: &[(usize, usize)]) -> Vec<Bond> {
    bond_pairs.iter().map(|&(a, b)| Bond { a, b }).collect()
}

fn carbohydrate_bond_metadata(count: usize) -> Vec<BondMetadata> {
    (0..count)
        .map(|_| BondMetadata {
            source: BondSource::Computed,
            flags: BondFlags::COVALENT,
            ..BondMetadata::computed()
        })
        .collect()
}

fn carbohydrate_reference_summary_json(molecule: &Molecule) -> String {
    let structure = molecule.atomic_structure();
    let carbohydrates = &structure.carbohydrates;
    let elements = carbohydrates
        .elements
        .iter()
        .enumerate()
        .map(|(index, element)| {
            let unit = &structure.units[element.unit_id];
            let ring_atoms = element
                .ring_element_indices
                .iter()
                .map(|&unit_index| {
                    let source_atom_index = unit.atom_indices[unit_index];
                    let atom = &molecule.atoms[source_atom_index];
                    format!(
                        "{{\"source_atom_index\":{},\"atom_id\":{},\"atom_name\":\"{}\"}}",
                        source_atom_index,
                        atom.id,
                        crate::json::json_escape(&atom.name)
                    )
                })
                .collect::<Vec<_>>()
                .join(",");
            let shape = get_saccharide_shape(
                element.component.component_type,
                element.ring_element_indices.len(),
            );
            format!(
                "{{\"index\":{},\"unit_id\":{},\"residue_index\":{},\"component\":{{\"abbr\":\"{}\",\"name\":\"{}\",\"type\":\"{:?}\",\"shape\":\"{:?}\",\"color\":{}}},\"ring_index\":{},\"ring_atoms\":[{}],\"alt_id\":\"{}\"}}",
                index,
                element.unit_id,
                element.residue_index,
                crate::json::json_escape(&element.component.abbr),
                crate::json::json_escape(&element.component.name),
                element.component.component_type,
                shape,
                element.component.color,
                element.ring_index,
                ring_atoms,
                crate::json::json_escape(&element.alt_id)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let partial_elements = carbohydrates
        .partial_elements
        .iter()
        .map(|element| {
            format!(
                "{{\"unit_id\":{},\"residue_index\":{},\"component\":\"{}\"}}",
                element.unit_id,
                element.residue_index,
                crate::json::json_escape(&element.component.abbr)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let links = carbohydrates
        .links
        .iter()
        .map(|link| {
            format!(
                "[{},{}]",
                link.carbohydrate_index_a, link.carbohydrate_index_b
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let terminal_links = carbohydrates
        .terminal_links
        .iter()
        .map(|link| {
            format!(
                "{{\"carbohydrate_index\":{},\"element_unit_id\":{},\"element_index\":{},\"from_carbohydrate\":{}}}",
                link.carbohydrate_index,
                link.element_unit_id,
                link.element_index,
                link.from_carbohydrate
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let element_lookups = carbohydrate_lookup_summary(&structure, CarbohydrateLookupKind::Element);
    let link_lookups = carbohydrate_lookup_summary(&structure, CarbohydrateLookupKind::Link);
    let terminal_link_lookups =
        carbohydrate_lookup_summary(&structure, CarbohydrateLookupKind::TerminalLink);

    format!(
        "{{\"molstar_reference_commit\":\"{}\",\"molstar_sources\":[\"artifacts/molstar/src/mol-model/structure/structure/carbohydrates/compute.ts\",\"artifacts/molstar/src/mol-model/structure/structure/carbohydrates/data.ts\",\"artifacts/molstar/src/mol-model/structure/structure/carbohydrates/constants.ts\"],\"fixture\":\"tests/fixtures/cif/carbohydrate-branched.cif\",\"summary_name\":\"carbohydrate-branched-reference\",\"elements\":[{}],\"partial_elements\":[{}],\"links\":[{}],\"terminal_links\":[{}],\"lookups\":{{\"element\":{},\"link\":{},\"terminal_link\":{}}}}}",
        MOLSTAR_REFERENCE_COMMIT,
        elements,
        partial_elements,
        links,
        terminal_links,
        element_lookups,
        link_lookups,
        terminal_link_lookups
    )
}

#[derive(Clone, Copy)]
enum CarbohydrateLookupKind {
    Element,
    Link,
    TerminalLink,
}

fn carbohydrate_lookup_summary(
    structure: &AtomicStructure,
    lookup_kind: CarbohydrateLookupKind,
) -> String {
    let mut entries = Vec::new();
    for unit in &structure.units {
        for &source_atom_index in &unit.atom_indices {
            let indices = match lookup_kind {
                CarbohydrateLookupKind::Element => {
                    structure.carbohydrate_element_indices(unit.id, source_atom_index)
                }
                CarbohydrateLookupKind::Link => {
                    structure.carbohydrate_link_indices(unit.id, source_atom_index)
                }
                CarbohydrateLookupKind::TerminalLink => {
                    structure.carbohydrate_terminal_link_indices(unit.id, source_atom_index)
                }
            };
            if indices.is_empty() {
                continue;
            }
            entries.push(format!(
                "{{\"unit_id\":{},\"source_atom_index\":{},\"indices\":[{}]}}",
                unit.id,
                source_atom_index,
                indices
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            ));
        }
    }
    format!("[{}]", entries.join(","))
}

fn vec3(x: f32, y: f32, z: f32) -> Vec3 {
    Vec3 { x, y, z }
}

fn assert_vec3_close(actual: Vec3, expected: Vec3, tolerance: f32) {
    assert!(
        actual.distance(expected) <= tolerance,
        "actual {actual:?} expected {expected:?}"
    );
}

fn secondary_structure_assignment_summary_json(molecule: &Molecule) -> String {
    let structure = molecule.atomic_structure();
    let secondary = &structure.model.secondary_structure;
    let residues = structure
        .model
        .hierarchy
        .residues
        .iter()
        .enumerate()
        .map(|(residue_index, residue)| {
            let chain = &structure.model.hierarchy.chains[residue.chain_index];
            let secondary_type = secondary.residue_type(residue_index);
            let key = secondary.key[residue_index];
            let element = secondary.elements.get(key).copied();
            let type_name = if secondary_type.contains(SecondaryStructureType::HELIX) {
                "helix"
            } else if secondary_type.contains(SecondaryStructureType::BETA) {
                "sheet"
            } else {
                "coil"
            };
            let element_name = match element {
                Some(SecondaryStructureElement::Helix) => "helix",
                Some(SecondaryStructureElement::Sheet) => "sheet",
                _ => "none",
            };
            format!(
                concat!(
                    "    {{\n",
                    "      \"residue_index\": {},\n",
                    "      \"chain\": \"{}\",\n",
                    "      \"seq_id\": \"{}\",\n",
                    "      \"insertion_code\": \"{}\",\n",
                    "      \"type\": \"{}\",\n",
                    "      \"type_bits\": {},\n",
                    "      \"key\": {},\n",
                    "      \"element\": \"{}\"\n",
                    "    }}"
                ),
                residue_index,
                test_json_escape(&chain.id),
                test_json_escape(&residue.label_seq_id),
                test_json_escape(&residue.insertion_code),
                type_name,
                secondary_type.bits,
                key,
                element_name
            )
        })
        .collect::<Vec<_>>()
        .join(",\n");

    format!(
        concat!(
            "{{\n",
            "  \"molstar_reference_commit\": \"{}\",\n",
            "  \"molstar_sources\": [\n",
            "    \"artifacts/molstar/src/mol-model-formats/structure/property/secondary-structure.ts\",\n",
            "    \"artifacts/molstar/src/mol-model/structure/model/types.ts\"\n",
            "  ],\n",
            "  \"summary_name\": \"helix-sheet-coil-assignment\",\n",
            "  \"residues\": [\n",
            "{}\n",
            "  ]\n",
            "}}"
        ),
        MOLSTAR_REFERENCE_COMMIT, residues
    )
}

fn test_json_escape(value: &str) -> String {
    value
        .chars()
        .flat_map(|ch| match ch {
            '"' => "\\\"".chars().collect::<Vec<_>>(),
            '\\' => "\\\\".chars().collect::<Vec<_>>(),
            '\n' => "\\n".chars().collect::<Vec<_>>(),
            '\r' => "\\r".chars().collect::<Vec<_>>(),
            '\t' => "\\t".chars().collect::<Vec<_>>(),
            _ => vec![ch],
        })
        .collect()
}

fn sorted_pair(a: usize, b: usize) -> (usize, usize) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}
