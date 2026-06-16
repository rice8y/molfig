use super::binary::*;
use super::*;

fn table<'a>(tables: &'a [CifTable], name: &str) -> &'a CifTable {
    tables
        .iter()
        .find(|table| table.name == name)
        .unwrap_or_else(|| panic!("missing table {name}"))
}

fn column<'a>(table: &'a CifTable, name: &str) -> &'a ColumnData {
    let index = table
        .header_index(name)
        .unwrap_or_else(|| panic!("missing column {name}"));
    table.columns[index]
        .as_ref()
        .unwrap_or_else(|| panic!("missing decoded column {name}"))
}

fn field_name(header: &str) -> &str {
    header
        .split_once('.')
        .map(|(_, field)| field)
        .unwrap_or(header)
}

fn decoded_column_kind(data: &ColumnData) -> &'static str {
    match data {
        ColumnData::Int(_) => "int",
        ColumnData::Float(_) => "float",
        ColumnData::Str(_) => "string",
        ColumnData::Bytes(_) => "bytes",
        ColumnData::Masked(data, _) => decoded_column_kind(data),
    }
}

fn molstar_value_kind(data: &ColumnData, row: usize) -> i32 {
    match data {
        ColumnData::Masked(_, mask) => mask.get(row).copied().unwrap_or(0),
        _ => 0,
    }
}

fn molstar_field_str(data: &ColumnData, row: usize) -> String {
    match data {
        ColumnData::Masked(data, mask) => match mask.get(row).copied().unwrap_or(0) {
            0 => molstar_field_str(data, row),
            _ => String::new(),
        },
        _ => data.string_at(row),
    }
}

fn molstar_i32(data: &ColumnData, row: usize) -> Option<i32> {
    match data {
        ColumnData::Masked(data, _) => molstar_i32(data, row),
        _ => data.i32_at(row),
    }
}

fn molstar_f32(data: &ColumnData, row: usize) -> Option<f32> {
    match data {
        ColumnData::Masked(data, _) => molstar_f32(data, row),
        _ => data.f32_at(row),
    }
}

fn summarize_molstar_typed_column(
    out: &mut String,
    fixture: &str,
    table: &CifTable,
    column_name: &str,
    rows: &[usize],
) {
    let data = column(table, column_name);
    out.push_str(&format!(
        "typed|{fixture}|{}.{}|kind={}|rows={}|valueKind={}|str={}",
        table.name,
        field_name(column_name),
        decoded_column_kind(data),
        table.row_count(),
        rows.iter()
            .map(|row| molstar_value_kind(data, *row).to_string())
            .collect::<Vec<_>>()
            .join(","),
        rows.iter()
            .map(|row| molstar_field_str(data, *row))
            .collect::<Vec<_>>()
            .join(",")
    ));
    match decoded_column_kind(data) {
        "int" | "bytes" => out.push_str(&format!(
            "|int={}",
            rows.iter()
                .map(|row| {
                    molstar_i32(data, *row)
                        .map(|value| value.to_string())
                        .unwrap_or_default()
                })
                .collect::<Vec<_>>()
                .join(",")
        )),
        "float" => out.push_str(&format!(
            "|float={}",
            rows.iter()
                .map(|row| {
                    molstar_f32(data, *row)
                        .map(|value| format!("{value}"))
                        .unwrap_or_default()
                })
                .collect::<Vec<_>>()
                .join(",")
        )),
        _ => {}
    }
    out.push('\n');
}

fn summarize_molstar_fixture(
    out: &mut String,
    fixture: &str,
    tables: &[CifTable],
    typed_checks: &[(&str, &str, &[usize])],
) {
    out.push_str(&format!("fixture|{fixture}\n"));
    for table in tables {
        out.push_str(&format!(
            "category|{fixture}|{}|rows={}|columns={}\n",
            table.name,
            table.row_count(),
            table
                .headers
                .iter()
                .map(|header| field_name(header))
                .collect::<Vec<_>>()
                .join(",")
        ));
    }
    for (table_name, column_name, rows) in typed_checks {
        summarize_molstar_typed_column(out, fixture, table(tables, table_name), column_name, rows);
    }
}

fn molstar_binary_cif_golden_summary() -> String {
    let assembly_fixture = "assembly-altloc-helix.bcif";
    let assembly_tables = parse_binary_cif_tables(include_bytes!(
        "../../tests/fixtures/bcif/assembly-altloc-helix.bcif"
    ))
    .unwrap();
    let ihm_fixture = "ihm-only.bcif";
    let ihm_tables =
        parse_binary_cif_tables(include_bytes!("../../tests/fixtures/bcif/ihm-only.bcif")).unwrap();

    let mut out = String::new();
    out.push_str("molstar_commit|1b8117d3f10f7c978aabb5a0d3d47370635aefe4\n");
    out.push_str("molstar_source|artifacts/molstar/src/mol-io/reader/cif/binary/parser.ts:21-52\n");
    out.push_str("molstar_source|artifacts/molstar/src/mol-io/reader/cif/binary/field.ts:15-57\n");
    out.push_str("molstar_source|artifacts/molstar/src/mol-data/db/column.ts:100-115\n");
    summarize_molstar_fixture(
        &mut out,
        assembly_fixture,
        &assembly_tables,
        &[
            ("entry", "_entry.id", &[0]),
            (
                "pdbx_entity_branch_link",
                "_pdbx_entity_branch_link.link_id",
                &[0],
            ),
            ("atom_site", "_atom_site.id", &[0, 1, 13]),
            ("atom_site", "_atom_site.Cartn_x", &[0, 1, 13]),
            ("atom_site", "_atom_site.auth_atom_id", &[0, 1, 13]),
            ("atom_site", "_atom_site.pdbx_PDB_ins_code", &[0, 6, 13]),
            ("chem_comp", "_chem_comp.formula_weight", &[0, 1, 2]),
            ("chem_comp_bond", "_chem_comp_bond.pdbx_ordinal", &[0, 1, 3]),
        ],
    );
    summarize_molstar_fixture(
        &mut out,
        ihm_fixture,
        &ihm_tables,
        &[
            ("ihm_model_list", "_ihm_model_list.model_id", &[0]),
            (
                "ihm_sphere_obj_site",
                "_ihm_sphere_obj_site.object_radius",
                &[0],
            ),
            (
                "ihm_cross_link_restraint",
                "_ihm_cross_link_restraint.distance_threshold",
                &[0],
            ),
        ],
    );
    let masked = ColumnData::Int(vec![11, 22, 33]).with_mask(&[0, 1, 2]);
    out.push_str(&format!(
        "mask|synthetic-int|valueKind={}|str={}|int={}\n",
        (0..3)
            .map(|row| molstar_value_kind(&masked, row).to_string())
            .collect::<Vec<_>>()
            .join(","),
        (0..3)
            .map(|row| molstar_field_str(&masked, row))
            .collect::<Vec<_>>()
            .join(","),
        (0..3)
            .map(|row| molstar_i32(&masked, row).unwrap_or_default().to_string())
            .collect::<Vec<_>>()
            .join(",")
    ));
    out
}

#[test]
fn binary_cif_decoded_columns_match_molstar_golden() {
    let actual = molstar_binary_cif_golden_summary();
    let expected = include_str!("../../tests/fixtures/bcif/molstar-decoded-columns.golden");
    assert_eq!(actual, expected, "\n{actual}");
}

#[test]
fn binary_cif_fixture_preserves_typed_numeric_and_string_columns() {
    let tables = parse_binary_cif_tables(include_bytes!(
        "../../tests/fixtures/bcif/assembly-altloc-helix.bcif"
    ))
    .unwrap();
    assert_eq!(
        tables
            .iter()
            .map(|table| table.name.as_str())
            .collect::<Vec<_>>()[..6],
        [
            "entry",
            "exptl",
            "entity",
            "entity_poly",
            "entity_poly_seq",
            "struct_asym"
        ]
    );

    let entry = table(&tables, "entry");
    assert_eq!(entry.row_count(), 1);
    assert_eq!(entry.headers, vec!["_entry.id"]);
    assert!(matches!(column(entry, "_entry.id"), ColumnData::Str(_)));

    let branch_link = table(&tables, "pdbx_entity_branch_link");
    assert!(matches!(
        column(branch_link, "_pdbx_entity_branch_link.link_id"),
        ColumnData::Int(_)
    ));
    assert!(matches!(
        column(
            branch_link,
            "_pdbx_entity_branch_link.entity_branch_list_num_2"
        ),
        ColumnData::Int(_)
    ));
    assert_eq!(
        branch_link.i32_at(
            0,
            branch_link
                .header_index("_pdbx_entity_branch_link.link_id")
                .unwrap()
        ),
        Some(1)
    );
    let branch_scheme = table(&tables, "pdbx_branch_scheme");
    assert!(matches!(
        column(branch_scheme, "_pdbx_branch_scheme.num"),
        ColumnData::Int(_)
    ));
    assert!(matches!(
        column(branch_scheme, "_pdbx_branch_scheme.auth_seq_num"),
        ColumnData::Str(_)
    ));
    let poly_scheme = table(&tables, "pdbx_poly_seq_scheme");
    assert!(matches!(
        column(poly_scheme, "_pdbx_poly_seq_scheme.seq_id"),
        ColumnData::Int(_)
    ));

    let atom_site = table(&tables, "atom_site");
    assert_eq!(atom_site.row_count(), 14);
    assert_eq!(atom_site.headers.len(), 22);
    assert_eq!(atom_site.headers[0], "_atom_site.group_PDB");
    assert_eq!(atom_site.headers[1], "_atom_site.id");
    assert!(matches!(
        column(atom_site, "_atom_site.id"),
        ColumnData::Int(_)
    ));
    assert!(matches!(
        column(atom_site, "_atom_site.Cartn_x"),
        ColumnData::Float(_)
    ));
    assert!(matches!(
        column(atom_site, "_atom_site.occupancy"),
        ColumnData::Float(_)
    ));
    assert!(matches!(
        column(atom_site, "_atom_site.pdbx_formal_charge"),
        ColumnData::Int(_)
    ));
    assert!(matches!(
        column(atom_site, "_atom_site.auth_atom_id"),
        ColumnData::Str(_)
    ));
    assert!(matches!(
        column(atom_site, "_atom_site.auth_comp_id"),
        ColumnData::Str(_)
    ));
    assert!(matches!(
        column(atom_site, "_atom_site.auth_asym_id"),
        ColumnData::Str(_)
    ));
    assert!(matches!(
        column(atom_site, "_atom_site.label_entity_id"),
        ColumnData::Str(_)
    ));
    assert!(matches!(
        column(atom_site, "_atom_site.auth_seq_id"),
        ColumnData::Str(_)
    ));
    assert!(matches!(
        column(atom_site, "_atom_site.pdbx_PDB_ins_code"),
        ColumnData::Str(_)
    ));
    assert!(matches!(
        column(atom_site, "_atom_site.ihm_model_id"),
        ColumnData::Int(_)
    ));
    assert_eq!(
        atom_site.usize_at(1, atom_site.header_index("_atom_site.id").unwrap()),
        Some(2)
    );
    assert_eq!(
        atom_site.float_at(1, atom_site.header_index("_atom_site.Cartn_x").unwrap()),
        Some(1.45)
    );
    assert_eq!(
        atom_site.i32_at(
            1,
            atom_site
                .header_index("_atom_site.pdbx_formal_charge")
                .unwrap()
        ),
        Some(1)
    );
    assert_eq!(
        atom_site.clean_at(
            1,
            atom_site.header_index("_atom_site.auth_atom_id").unwrap()
        ),
        "CAX"
    );
    assert_eq!(
        atom_site.clean_at(
            1,
            atom_site.header_index("_atom_site.auth_comp_id").unwrap()
        ),
        "ALAX"
    );
    assert_eq!(
        atom_site.clean_at(
            1,
            atom_site.header_index("_atom_site.auth_asym_id").unwrap()
        ),
        "X"
    );
    assert_eq!(
        atom_site.clean_at(
            1,
            atom_site
                .header_index("_atom_site.label_entity_id")
                .unwrap()
        ),
        "1"
    );
    assert_eq!(
        atom_site.clean_at(1, atom_site.header_index("_atom_site.auth_seq_id").unwrap()),
        "101"
    );
    assert_eq!(
        atom_site.clean_at(
            1,
            atom_site
                .header_index("_atom_site.pdbx_PDB_ins_code")
                .unwrap()
        ),
        "A"
    );

    let chem_comp = table(&tables, "chem_comp");
    assert!(matches!(
        column(chem_comp, "_chem_comp.formula_weight"),
        ColumnData::Float(_)
    ));
    assert_eq!(
        chem_comp.float_at(
            0,
            chem_comp.header_index("_chem_comp.formula_weight").unwrap()
        ),
        Some(89.09)
    );

    let chem_comp_atom = table(&tables, "chem_comp_atom");
    assert!(matches!(
        column(chem_comp_atom, "_chem_comp_atom.charge"),
        ColumnData::Int(_)
    ));
    assert!(matches!(
        column(chem_comp_atom, "_chem_comp_atom.model_Cartn_x"),
        ColumnData::Float(_)
    ));

    let chem_comp_bond = table(&tables, "chem_comp_bond");
    assert!(matches!(
        column(chem_comp_bond, "_chem_comp_bond.pdbx_ordinal"),
        ColumnData::Int(_)
    ));
    assert_eq!(
        chem_comp_bond.i32_at(
            3,
            chem_comp_bond
                .header_index("_chem_comp_bond.pdbx_ordinal")
                .unwrap()
        ),
        Some(104)
    );

    let chem_comp_angle = table(&tables, "chem_comp_angle");
    assert!(matches!(
        column(chem_comp_angle, "_chem_comp_angle.value_angle"),
        ColumnData::Float(_)
    ));
    assert_eq!(
        chem_comp_angle.float_at(
            0,
            chem_comp_angle
                .header_index("_chem_comp_angle.value_angle")
                .unwrap()
        ),
        Some(111.0)
    );

    let anisotrop = table(&tables, "atom_site_anisotrop");
    assert!(matches!(
        column(anisotrop, "_atom_site_anisotrop.id"),
        ColumnData::Int(_)
    ));
    assert!(matches!(
        column(anisotrop, "_atom_site_anisotrop.U[1][1]"),
        ColumnData::Float(_)
    ));

    let struct_conn = table(&tables, "struct_conn");
    assert_eq!(
        struct_conn.clean_at(
            0,
            struct_conn
                .header_index("_struct_conn.pdbx_ptnr1_PDB_ins_code")
                .unwrap()
        ),
        "A"
    );
    assert_eq!(
        struct_conn.clean_at(
            8,
            struct_conn
                .header_index("_struct_conn.pdbx_ptnr2_PDB_ins_code")
                .unwrap()
        ),
        "B"
    );

    let struct_conf = table(&tables, "struct_conf");
    assert_eq!(
        struct_conf.clean_at(
            0,
            struct_conf
                .header_index("_struct_conf.pdbx_beg_PDB_ins_code")
                .unwrap()
        ),
        "A"
    );
    assert_eq!(
        struct_conf.clean_at(
            0,
            struct_conf
                .header_index("_struct_conf.pdbx_end_PDB_ins_code")
                .unwrap()
        ),
        "B"
    );

    let struct_sheet_range = table(&tables, "struct_sheet_range");
    assert_eq!(
        struct_sheet_range.clean_at(
            0,
            struct_sheet_range
                .header_index("_struct_sheet_range.pdbx_beg_PDB_ins_code")
                .unwrap()
        ),
        ""
    );
    assert_eq!(
        struct_sheet_range.clean_at(
            0,
            struct_sheet_range
                .header_index("_struct_sheet_range.pdbx_end_PDB_ins_code")
                .unwrap()
        ),
        "B"
    );
}

#[test]
fn binary_cif_ihm_metadata_preserves_typed_columns() {
    let tables =
        parse_binary_cif_tables(include_bytes!("../../tests/fixtures/bcif/ihm-only.bcif")).unwrap();
    assert_eq!(
        tables
            .iter()
            .map(|table| table.name.as_str())
            .collect::<Vec<_>>()[..4],
        [
            "ihm_model_list",
            "ihm_model_group",
            "ihm_model_group_link",
            "ihm_sphere_obj_site"
        ]
    );

    let model_list = table(&tables, "ihm_model_list");
    assert_eq!(
        column(model_list, "_ihm_model_list.model_id"),
        &ColumnData::Int(vec![1])
    );
    assert!(matches!(
        column(model_list, "_ihm_model_list.model_name"),
        ColumnData::Str(_)
    ));
    assert_eq!(
        model_list.i32_at(
            0,
            model_list.header_index("_ihm_model_list.model_id").unwrap()
        ),
        Some(1)
    );

    let group = table(&tables, "ihm_model_group");
    assert!(matches!(
        column(group, "_ihm_model_group.id"),
        ColumnData::Int(_)
    ));

    let link = table(&tables, "ihm_model_group_link");
    assert!(matches!(
        column(link, "_ihm_model_group_link.model_id"),
        ColumnData::Int(_)
    ));
    assert!(matches!(
        column(link, "_ihm_model_group_link.group_id"),
        ColumnData::Int(_)
    ));

    let restraint = table(&tables, "ihm_cross_link_restraint");
    assert!(matches!(
        column(restraint, "_ihm_cross_link_restraint.id"),
        ColumnData::Int(_)
    ));
    assert_eq!(
        column(restraint, "_ihm_cross_link_restraint.distance_threshold"),
        &ColumnData::Float(vec![25.0])
    );
    assert_eq!(
        restraint.float_at(
            0,
            restraint
                .header_index("_ihm_cross_link_restraint.distance_threshold")
                .unwrap()
        ),
        Some(25.0)
    );
}

#[test]
fn binary_cif_masked_values_preserve_present_unknown_and_not_present_semantics() {
    let values = ColumnData::Int(vec![11, 22, 33]).with_mask(&[0, 1, 2]);

    assert_eq!(values.string_at(0), "11");
    assert_eq!(values.string_at(1), ".");
    assert_eq!(values.string_at(2), "?");
    assert_eq!(values.i32_at(0), Some(11));
    assert_eq!(values.i32_at(1), None);
    assert_eq!(values.i32_at(2), None);
    assert_eq!(values.usize_at(0), Some(11));
    assert_eq!(values.usize_at(1), None);
    assert_eq!(values.f32_at(2), None);
}

#[test]
fn cif_text_tokens_preserve_quotes_multiline_and_missing_values_like_molstar() {
    let tokens = cif_tokens(
            "data_demo\r\n_entity.id 'A B # C'\r\n_entity.type .\r\n_entity.pdbx_description\r\n;\r\nfirst line\r\n  second line\r\n;\r\n_entity.formula ?\r\n",
        );
    let tables = cif_tables(&tokens).unwrap();
    let entity = table(&tables, "entity");

    assert_eq!(
        entity.raw_at(0, entity.header_index("_entity.id").unwrap()),
        "A B # C"
    );
    assert_eq!(
        entity.clean_at(0, entity.header_index("_entity.type").unwrap()),
        ""
    );
    assert_eq!(
        entity.raw_at(0, entity.header_index("_entity.pdbx_description").unwrap()),
        "\r\nfirst line\r\n  second line"
    );
    assert_eq!(
        entity.clean_at(0, entity.header_index("_entity.formula").unwrap()),
        ""
    );
}

#[test]
fn binary_cif_decodes_all_molstar_encoding_kinds() {
    let byte_cases = [
        (vec![0xff], 1, ColumnData::Int(vec![-1])),
        (
            i16::to_le_bytes(-258).to_vec(),
            2,
            ColumnData::Int(vec![-258]),
        ),
        (
            i32::to_le_bytes(-65_536).to_vec(),
            3,
            ColumnData::Int(vec![-65_536]),
        ),
        (vec![255], 4, ColumnData::Int(vec![255])),
        (
            u16::to_le_bytes(65_530).to_vec(),
            5,
            ColumnData::Int(vec![65_530]),
        ),
        (
            u32::to_le_bytes(65_530).to_vec(),
            6,
            ColumnData::Int(vec![65_530]),
        ),
        (
            f32::to_le_bytes(1.25).to_vec(),
            32,
            ColumnData::Float(vec![1.25]),
        ),
        (
            f64::to_le_bytes(2.5).to_vec(),
            33,
            ColumnData::Float(vec![2.5]),
        ),
    ];
    for (bytes, ty, expected) in byte_cases {
        assert_eq!(decode_byte_array(&bytes, ty).unwrap(), expected);
    }

    assert_eq!(
        decode_bcif_step(
            ColumnData::Int(vec![15, -25]),
            &encoding("FixedPoint", &[("factor", MpValue::Int(10))])
        )
        .unwrap(),
        ColumnData::Float(vec![1.5, -2.5])
    );
    assert_eq!(
        decode_bcif_step(
            ColumnData::Int(vec![0, 2, 4]),
            &encoding(
                "IntervalQuantization",
                &[
                    ("min", MpValue::Float(-1.0)),
                    ("max", MpValue::Float(1.0)),
                    ("numSteps", MpValue::Int(5)),
                ]
            )
        )
        .unwrap(),
        ColumnData::Float(vec![-1.0, 0.0, 1.0])
    );
    assert_eq!(
        decode_bcif_step(
            ColumnData::Int(vec![5, 3, 2, 2]),
            &encoding("RunLength", &[("srcSize", MpValue::Int(5))])
        )
        .unwrap(),
        ColumnData::Int(vec![5, 5, 5, 2, 2])
    );
    assert_eq!(
        decode_bcif_step(
            ColumnData::Int(vec![1, 2, -3]),
            &encoding("Delta", &[("origin", MpValue::Int(10))])
        )
        .unwrap(),
        ColumnData::Int(vec![11, 13, 10])
    );
    assert_eq!(
        decode_bcif_step(
            ColumnData::Int(vec![127, 5, -128, -2, 3]),
            &encoding(
                "IntegerPacking",
                &[
                    ("byteCount", MpValue::Int(1)),
                    ("isUnsigned", MpValue::Bool(false)),
                    ("srcSize", MpValue::Int(3)),
                ]
            )
        )
        .unwrap(),
        ColumnData::Int(vec![132, -130, 3])
    );

    let string_encoding = encoding(
        "StringArray",
        &[
            ("stringData", MpValue::Str("ABB".to_string())),
            (
                "offsetEncoding",
                MpValue::Array(vec![encoding("ByteArray", &[("type", MpValue::Int(4))])]),
            ),
            ("offsets", MpValue::Bin(vec![0, 1, 3])),
            (
                "dataEncoding",
                MpValue::Array(vec![encoding("ByteArray", &[("type", MpValue::Int(4))])]),
            ),
        ],
    );
    assert_eq!(
        decode_bcif_step(ColumnData::Bytes(vec![0, 1, 0]), &string_encoding).unwrap(),
        ColumnData::Str(vec!["A".to_string(), "BB".to_string(), "A".to_string()])
    );
}

#[test]
fn binary_cif_chained_encoding_applies_in_molstar_order() {
    let encoded = MpValue::Map(vec![
        (
            "encoding".to_string(),
            MpValue::Array(vec![
                encoding("FixedPoint", &[("factor", MpValue::Int(2))]),
                encoding("Delta", &[("origin", MpValue::Int(10))]),
                encoding("ByteArray", &[("type", MpValue::Int(4))]),
            ]),
        ),
        ("data".to_string(), MpValue::Bin(vec![1, 2, 3])),
    ]);

    assert_eq!(
        decode_bcif_data(&encoded).unwrap(),
        ColumnData::Float(vec![5.5, 6.5, 8.0])
    );
}

#[test]
fn binary_cif_block_layer_preserves_all_blocks_and_default_first_block() {
    let value = MpValue::Map(vec![(
        "dataBlocks".to_string(),
        MpValue::Array(vec![
            binary_cif_test_block("first", "A"),
            binary_cif_test_block("second", "B"),
        ]),
    )]);

    let blocks = binary_cif_blocks_from_value(&value).unwrap();
    assert_eq!(blocks.len(), 2);
    assert_eq!(blocks[0].header, "first");
    assert_eq!(blocks[1].header, "second");
    assert_eq!(
        table(&blocks[0].tables, "entry").clean_at(0, 0),
        "A".to_string()
    );
    assert_eq!(
        table(&blocks[1].tables, "entry").clean_at(0, 0),
        "B".to_string()
    );
    assert_eq!(
        table(
            &select_cif_block_tables(blocks.clone(), CifBlockSelection::First)
                .unwrap()
                .tables,
            "entry"
        )
        .clean_at(0, 0),
        "A".to_string()
    );
    assert_eq!(
        table(
            &select_cif_block_tables(blocks.clone(), CifBlockSelection::Index(1))
                .unwrap()
                .tables,
            "entry"
        )
        .clean_at(0, 0),
        "B".to_string()
    );
    assert_eq!(
        table(
            &select_cif_block_tables(blocks.clone(), CifBlockSelection::Header("second"))
                .unwrap()
                .tables,
            "entry"
        )
        .clean_at(0, 0),
        "B".to_string()
    );
    assert_eq!(
        binary_cif_block_selection(Some("second"), Some(0)),
        CifBlockSelection::Header("second")
    );
    assert_eq!(
        select_cif_block_tables(blocks.clone(), CifBlockSelection::Index(2)).unwrap_err(),
        "BinaryCIF has no data block at index 2"
    );
    assert_eq!(
        select_cif_block_tables(blocks, CifBlockSelection::Header("missing")).unwrap_err(),
        "BinaryCIF has no data block named missing"
    );
}

#[test]
fn binary_cif_nil_mask_is_treated_as_absent() {
    let value = MpValue::Map(vec![(
        "dataBlocks".to_string(),
        MpValue::Array(vec![binary_cif_test_block_with_mask(
            "nil-mask",
            "1CRN",
            Some(MpValue::Nil),
        )]),
    )]);

    let blocks = binary_cif_blocks_from_value(&value).unwrap();
    let entry = table(&blocks[0].tables, "entry");
    assert_eq!(entry.clean_at(0, 0), "1CRN");
    assert!(matches!(column(entry, "_entry.id"), ColumnData::Str(_)));
}

#[test]
fn binary_cif_real_1crn_fixture_accepts_nil_masks() {
    let tables =
        parse_binary_cif_tables(include_bytes!("../../../package/examples/data/1crn.bcif"))
            .unwrap();
    let entry = table(&tables, "entry");
    assert_eq!(entry.clean_at(0, 0), "1CRN");
    assert!(table(&tables, "atom_site").rows.len() > 300);
}

#[test]
fn binary_cif_real_1crn_fixture_parses_to_molecule() {
    let molecule = parse_molecule(
        include_bytes!("../../../package/examples/data/1crn.bcif"),
        crate::options::InputFormat::BinaryCif,
    )
    .unwrap();
    assert_eq!(molecule.source_data.name, "1CRN");
    assert_eq!(molecule.atoms.len(), 327);
    assert!(!molecule.bonds.is_empty());
}

#[test]
fn binary_cif_reports_unsupported_byte_array_type() {
    let err = decode_byte_array(&[0, 1, 2, 3], 127).unwrap_err();
    assert_eq!(err, "unsupported BinaryCIF byte array type: 127");
}

#[test]
fn binary_cif_reports_misaligned_byte_array_length() {
    let err = decode_byte_array(&[0, 1, 2], 3).unwrap_err();
    assert_eq!(
        err,
        "BinaryCIF byte array type 3 has 3 bytes, not divisible by 4"
    );
}

#[test]
fn binary_cif_reports_invalid_mask_lengths_and_values() {
    let short = validate_bcif_mask("atom_site", "Cartn_x", 3, &[0, 1]).unwrap_err();
    assert_eq!(
        short,
        "BinaryCIF column atom_site.Cartn_x mask has 2 rows, expected 3"
    );

    let invalid = validate_bcif_mask("atom_site", "Cartn_x", 3, &[0, 3, 2]).unwrap_err();
    assert_eq!(
        invalid,
        "BinaryCIF column atom_site.Cartn_x mask has unsupported value 3"
    );
}

#[test]
fn binary_cif_reports_unsupported_encoding_kind() {
    let encoding = MpValue::Map(vec![(
        "kind".to_string(),
        MpValue::Str("ImaginaryEncoding".to_string()),
    )]);

    let err = decode_bcif_step(ColumnData::Bytes(vec![1, 2, 3]), &encoding).unwrap_err();
    assert_eq!(
        err,
        "unsupported BinaryCIF encoding kind: ImaginaryEncoding"
    );
}

#[test]
fn binary_cif_reports_malformed_messagepack() {
    let err = parse_binary_cif_tables(&[0xde]).unwrap_err();
    assert_eq!(err, "MessagePack ended unexpectedly");
}

#[test]
fn binary_cif_reports_trailing_messagepack_bytes() {
    let err = MsgPack::parse(&[0xc0, 0xc0]).unwrap_err();
    assert_eq!(err, "MessagePack has trailing bytes");
}

fn encoding(kind: &str, entries: &[(&str, MpValue)]) -> MpValue {
    let mut map = vec![("kind".to_string(), MpValue::Str(kind.to_string()))];
    map.extend(
        entries
            .iter()
            .map(|(key, value)| ((*key).to_string(), value.clone())),
    );
    MpValue::Map(map)
}

fn binary_cif_test_block(header: &str, entry_id: &str) -> MpValue {
    binary_cif_test_block_with_mask(header, entry_id, None)
}

fn binary_cif_test_block_with_mask(header: &str, entry_id: &str, mask: Option<MpValue>) -> MpValue {
    let mut column = vec![
        ("name".to_string(), MpValue::Str("id".to_string())),
        (
            "data".to_string(),
            encoding(
                "",
                &[
                    (
                        "encoding",
                        MpValue::Array(vec![encoding(
                            "StringArray",
                            &[
                                ("stringData", MpValue::Str(entry_id.to_string())),
                                (
                                    "offsetEncoding",
                                    MpValue::Array(vec![encoding(
                                        "ByteArray",
                                        &[("type", MpValue::Int(4))],
                                    )]),
                                ),
                                ("offsets", MpValue::Bin(vec![0, entry_id.len() as u8])),
                                (
                                    "dataEncoding",
                                    MpValue::Array(vec![encoding(
                                        "ByteArray",
                                        &[("type", MpValue::Int(4))],
                                    )]),
                                ),
                            ],
                        )]),
                    ),
                    ("data", MpValue::Bin(vec![0])),
                ],
            ),
        ),
    ];
    if let Some(mask) = mask {
        column.push(("mask".to_string(), mask));
    }

    MpValue::Map(vec![
        ("header".to_string(), MpValue::Str(header.to_string())),
        (
            "categories".to_string(),
            MpValue::Array(vec![MpValue::Map(vec![
                ("name".to_string(), MpValue::Str("entry".to_string())),
                ("rowCount".to_string(), MpValue::Int(1)),
                (
                    "columns".to_string(),
                    MpValue::Array(vec![MpValue::Map(column)]),
                ),
            ])]),
        ),
    ])
}
