use super::ColumnData;

#[derive(Clone, Debug)]
pub(super) struct CifTable {
    pub(super) name: String,
    pub(super) headers: Vec<String>,
    pub(super) rows: Vec<Vec<String>>,
    pub(super) columns: Vec<Option<ColumnData>>,
}

impl CifTable {
    pub(super) fn header_index(&self, name: &str) -> Option<usize> {
        self.headers
            .iter()
            .position(|header| header.eq_ignore_ascii_case(name))
    }

    pub(super) fn row_count(&self) -> usize {
        self.rows.len()
    }

    pub(super) fn row_indices(&self) -> std::ops::Range<usize> {
        0..self.row_count()
    }

    pub(super) fn raw_at(&self, row: usize, column: usize) -> String {
        self.columns
            .get(column)
            .and_then(Option::as_ref)
            .map(|data| data.string_at(row))
            .or_else(|| {
                self.rows
                    .get(row)
                    .and_then(|fields| fields.get(column))
                    .cloned()
            })
            .unwrap_or_default()
    }

    pub(super) fn clean_at(&self, row: usize, column: usize) -> String {
        clean_cif_value(&self.raw_at(row, column))
    }

    pub(super) fn float_at(&self, row: usize, column: usize) -> Option<f32> {
        self.columns
            .get(column)
            .and_then(Option::as_ref)
            .and_then(|data| data.f32_at(row))
            .or_else(|| {
                self.rows
                    .get(row)
                    .and_then(|fields| fields.get(column))
                    .and_then(|value| parse_js_number_f32(&clean_cif_value(value)))
            })
    }

    pub(super) fn usize_at(&self, row: usize, column: usize) -> Option<usize> {
        self.columns
            .get(column)
            .and_then(Option::as_ref)
            .and_then(|data| data.usize_at(row))
            .or_else(|| {
                self.rows
                    .get(row)
                    .and_then(|fields| fields.get(column))
                    .and_then(|value| clean_cif_value(value).parse::<usize>().ok())
            })
    }

    pub(super) fn i32_at(&self, row: usize, column: usize) -> Option<i32> {
        self.columns
            .get(column)
            .and_then(Option::as_ref)
            .and_then(|data| data.i32_at(row))
            .or_else(|| {
                self.rows
                    .get(row)
                    .and_then(|fields| fields.get(column))
                    .and_then(|value| clean_cif_value(value).parse::<i32>().ok())
            })
    }
}

pub(super) fn clean_nonempty_at(table: &CifTable, row: usize, column: usize) -> Option<String> {
    let value = table.clean_at(row, column);
    (!value.is_empty()).then_some(value)
}

fn column_has_present_value(table: &CifTable, column: usize) -> bool {
    table
        .row_indices()
        .any(|row| is_present_cif_value(&table.raw_at(row, column)))
}

pub(super) fn normalized_cif_column_pair(
    table: &CifTable,
    a: Option<usize>,
    b: Option<usize>,
) -> (Option<usize>, Option<usize>) {
    let a_present = a.is_some_and(|column| column_has_present_value(table, column));
    let b_present = b.is_some_and(|column| column_has_present_value(table, column));
    let a = if a_present { a } else { b };
    let b = if b_present { b } else { a };
    (a, b)
}

pub(super) fn cif_tables(tokens: &[String]) -> Result<Vec<CifTable>, String> {
    let mut tables = Vec::new();
    let mut i = 0;
    while i < tokens.len() {
        if tokens[i].starts_with("data_") {
            if !tables.is_empty() {
                break;
            }
            i += 1;
            continue;
        }

        if tokens[i] == "loop_" {
            i += 1;
            let header_start = i;
            while i < tokens.len() && tokens[i].starts_with("_") {
                i += 1;
            }
            let headers = tokens[header_start..i].to_vec();
            if headers.is_empty() {
                continue;
            }
            let width = headers.len();
            let mut rows = Vec::new();
            while i + width <= tokens.len() && !is_loop_boundary(&tokens[i]) {
                rows.push(tokens[i..i + width].to_vec());
                i += width;
            }
            let name = cif_category_name(&headers[0]);
            tables.push(CifTable {
                name,
                headers,
                rows,
                columns: vec![None; width],
            });
            continue;
        }

        if tokens[i].starts_with('_') {
            while i < tokens.len() && tokens[i].starts_with('_') {
                let header = tokens[i].clone();
                i += 1;
                if i >= tokens.len()
                    || tokens[i].starts_with('_')
                    || is_single_row_value_boundary(&tokens[i])
                {
                    return Err("Expected value.".to_string());
                }
                let value = tokens[i].clone();
                i += 1;
                merge_single_row_cif_field(&mut tables, header, value);
            }
            continue;
        }

        i += 1;
    }
    Ok(tables)
}

fn merge_single_row_cif_field(tables: &mut Vec<CifTable>, header: String, value: String) {
    let name = cif_category_name(&header);
    if name.is_empty() {
        return;
    }
    if let Some(table) = tables
        .iter_mut()
        .find(|table| table.name == name && table.row_count() == 1)
    {
        if let Some(index) = table
            .headers
            .iter()
            .position(|existing| existing.eq_ignore_ascii_case(&header))
        {
            table.rows[0][index] = value;
        } else {
            table.headers.push(header);
            table.rows[0].push(value);
            table.columns.push(None);
        }
    } else {
        tables.push(CifTable {
            name,
            headers: vec![header],
            rows: vec![vec![value]],
            columns: vec![None],
        });
    }
}

fn cif_category_name(header: &str) -> String {
    header
        .strip_prefix('_')
        .and_then(|h| h.split_once('.').map(|(cat, _)| cat.to_string()))
        .unwrap_or_default()
}

fn is_single_row_value_boundary(token: &str) -> bool {
    token == "loop_" || token == "stop_" || token.starts_with("data_") || token.starts_with("save_")
}

fn is_loop_boundary(token: &str) -> bool {
    token == "loop_"
        || token == "stop_"
        || token.starts_with("data_")
        || token.starts_with("save_")
        || token.starts_with("_")
}

pub(super) fn cif_tokens(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut chars = text.chars().peekable();
    let mut at_line_start = true;

    while let Some(ch) = chars.next() {
        if ch == '\n' || ch == '\r' {
            at_line_start = true;
            continue;
        }
        if ch.is_whitespace() {
            continue;
        }
        if ch == '#' {
            for c in chars.by_ref() {
                if c == '\n' {
                    at_line_start = true;
                    break;
                }
            }
            continue;
        }
        if at_line_start && ch == ';' {
            let mut value = String::new();
            let mut line_start = true;
            while let Some(c) = chars.next() {
                if line_start && c == ';' {
                    for c2 in chars.by_ref() {
                        if c2 == '\n' {
                            break;
                        }
                    }
                    break;
                }
                line_start = c == '\n' || c == '\r';
                value.push(c);
            }
            while value.ends_with('\n') || value.ends_with('\r') {
                value.pop();
            }
            out.push(value);
            at_line_start = true;
            continue;
        }
        at_line_start = false;
        if ch == '\'' || ch == '"' {
            let quote = ch;
            let mut value = String::new();
            for c in chars.by_ref() {
                if c == quote {
                    break;
                }
                value.push(c);
            }
            out.push(value);
        } else {
            let mut value = String::from(ch);
            while let Some(&c) = chars.peek() {
                if c.is_whitespace() || c == '#' {
                    break;
                }
                value.push(c);
                chars.next();
            }
            out.push(value);
        }
    }
    out
}

pub(super) fn clean_cif_value(value: &str) -> String {
    if matches!(value, "." | "?") {
        String::new()
    } else {
        value.to_string()
    }
}

pub(super) fn parse_js_number_f32(value: &str) -> Option<f32> {
    value.parse::<f64>().ok().map(|value| value as f32)
}

pub(super) fn is_present_cif_value(value: &str) -> bool {
    !matches!(value.trim(), "" | "." | "?")
}
