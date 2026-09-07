use super::cif_table::{clean_cif_value, parse_js_number_f32, CifTable};
use super::{parse_cif_tables, source_categories};
use crate::model::{Molecule, SourceData};
use crate::options::MeshOptions;

pub(super) fn parse_binary_cif(data: &[u8], options: &MeshOptions) -> Result<Molecule, String> {
    let block = parse_binary_cif_block_with_selection(
        data,
        binary_cif_block_selection(options.block_header.as_deref(), options.block_index),
    )?;
    let source_data = SourceData::mmcif(block.header.clone(), source_categories(&block.tables));
    parse_cif_tables(&block.tables, source_data)
}

#[allow(dead_code)]
pub(super) fn parse_binary_cif_tables(data: &[u8]) -> Result<Vec<CifTable>, String> {
    parse_binary_cif_tables_with_selection(data, CifBlockSelection::First)
}

#[allow(dead_code)]
pub(super) fn parse_binary_cif_tables_by_index(
    data: &[u8],
    index: usize,
) -> Result<Vec<CifTable>, String> {
    parse_binary_cif_tables_with_selection(data, CifBlockSelection::Index(index))
}

#[allow(dead_code)]
pub(super) fn parse_binary_cif_tables_by_header(
    data: &[u8],
    header: &str,
) -> Result<Vec<CifTable>, String> {
    parse_binary_cif_tables_with_selection(data, CifBlockSelection::Header(header))
}

pub(super) fn parse_binary_cif_tables_with_selection(
    data: &[u8],
    selection: CifBlockSelection<'_>,
) -> Result<Vec<CifTable>, String> {
    Ok(parse_binary_cif_block_with_selection(data, selection)?.tables)
}

fn parse_binary_cif_block_with_selection(
    data: &[u8],
    selection: CifBlockSelection<'_>,
) -> Result<CifBlockTables, String> {
    let value = MsgPack::parse(data)?;
    let blocks = binary_cif_blocks_from_value(&value)?;
    select_cif_block_tables(blocks, selection)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum CifBlockSelection<'a> {
    First,
    Index(usize),
    Header(&'a str),
}

pub(super) fn binary_cif_block_selection(
    header: Option<&str>,
    index: Option<usize>,
) -> CifBlockSelection<'_> {
    if let Some(header) = header.filter(|value| !value.is_empty()) {
        CifBlockSelection::Header(header)
    } else if let Some(index) = index {
        CifBlockSelection::Index(index)
    } else {
        CifBlockSelection::First
    }
}

pub(super) fn select_cif_block_tables(
    blocks: Vec<CifBlockTables>,
    selection: CifBlockSelection<'_>,
) -> Result<CifBlockTables, String> {
    match selection {
        CifBlockSelection::First => blocks
            .into_iter()
            .next()
            .ok_or_else(|| "BinaryCIF has no data blocks".to_string()),
        CifBlockSelection::Index(index) => blocks
            .into_iter()
            .nth(index)
            .ok_or_else(|| format!("BinaryCIF has no data block at index {index}")),
        CifBlockSelection::Header(header) => blocks
            .into_iter()
            .find(|block| block.header == header)
            .ok_or_else(|| format!("BinaryCIF has no data block named {header}")),
    }
}

#[derive(Clone, Debug)]
pub(super) struct CifBlockTables {
    pub(super) header: String,
    pub(super) tables: Vec<CifTable>,
}

pub(super) fn binary_cif_blocks_from_value(value: &MpValue) -> Result<Vec<CifBlockTables>, String> {
    let blocks = value
        .get("dataBlocks")
        .and_then(MpValue::as_array)
        .ok_or("BinaryCIF is missing dataBlocks")?;
    if blocks.is_empty() {
        return Err("BinaryCIF has no data blocks".to_string());
    }

    blocks.iter().map(binary_cif_block_from_value).collect()
}

pub(super) fn binary_cif_block_from_value(block: &MpValue) -> Result<CifBlockTables, String> {
    let header = block
        .get("header")
        .and_then(MpValue::as_str)
        .unwrap_or("")
        .to_string();
    let categories = block
        .get("categories")
        .and_then(MpValue::as_array)
        .ok_or("BinaryCIF block is missing categories")?;

    let mut tables = Vec::new();
    for category in categories {
        let name = category
            .get("name")
            .and_then(MpValue::as_str)
            .unwrap_or("")
            .trim_start_matches('_')
            .to_string();
        let row_count = category
            .get("rowCount")
            .and_then(MpValue::as_i64)
            .unwrap_or(0)
            .max(0) as usize;
        let columns = category
            .get("columns")
            .and_then(MpValue::as_array)
            .ok_or("BinaryCIF category is missing columns")?;
        let mut headers = Vec::new();
        let mut decoded = Vec::<(String, ColumnData)>::new();
        for column in columns {
            let col_name = column
                .get("name")
                .and_then(MpValue::as_str)
                .ok_or("BinaryCIF column is missing name")?;
            let data_value = column
                .get("data")
                .ok_or("BinaryCIF column is missing data")?;
            let mut values = decode_bcif_data(data_value)
                .map_err(|err| format!("BinaryCIF column {name}.{col_name} data: {err}"))?;
            if let Some(mask_value) = column
                .get("mask")
                .filter(|value| !matches!(value, MpValue::Nil))
            {
                let mask = decode_bcif_data(mask_value)
                    .map_err(|err| format!("BinaryCIF column {name}.{col_name} mask: {err}"))?
                    .to_i32_vec();
                validate_bcif_mask(&name, col_name, row_count, &mask)?;
                values = values.with_mask(&mask);
            }
            if values.len() != row_count {
                return Err(format!(
                    "BinaryCIF column {}.{} has {} rows, expected {}",
                    name,
                    col_name,
                    values.len(),
                    row_count
                ));
            }
            headers.push(format!("_{}.{}", name, col_name));
            decoded.push((col_name.to_string(), values));
        }
        let mut rows = vec![Vec::<String>::new(); row_count];
        for (_, values) in &decoded {
            for (row, row_values) in rows.iter_mut().enumerate().take(row_count) {
                row_values.push(values.string_at(row));
            }
        }
        tables.push(CifTable {
            name,
            headers,
            rows,
            columns: decoded
                .into_iter()
                .map(|(_, values)| Some(values))
                .collect(),
        });
    }

    Ok(CifBlockTables { header, tables })
}

pub(super) fn validate_bcif_mask(
    category_name: &str,
    column_name: &str,
    row_count: usize,
    mask: &[i32],
) -> Result<(), String> {
    if mask.len() != row_count {
        return Err(format!(
            "BinaryCIF column {category_name}.{column_name} mask has {} rows, expected {row_count}",
            mask.len()
        ));
    }
    if let Some(value) = mask.iter().find(|value| !matches!(value, 0..=2)) {
        return Err(format!(
            "BinaryCIF column {category_name}.{column_name} mask has unsupported value {value}"
        ));
    }
    Ok(())
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum ColumnData {
    Int(Vec<i32>),
    Float(Vec<f32>),
    Str(Vec<String>),
    Bytes(Vec<u8>),
    Masked(Box<ColumnData>, Vec<i32>),
}

impl ColumnData {
    fn len(&self) -> usize {
        match self {
            Self::Int(v) => v.len(),
            Self::Float(v) => v.len(),
            Self::Str(v) => v.len(),
            Self::Bytes(v) => v.len(),
            Self::Masked(data, _) => data.len(),
        }
    }

    pub(crate) fn string_at(&self, row: usize) -> String {
        match self {
            Self::Int(v) => v.get(row).copied().unwrap_or_default().to_string(),
            Self::Float(v) => format!("{}", v.get(row).copied().unwrap_or_default()),
            Self::Str(v) => v.get(row).cloned().unwrap_or_default(),
            Self::Bytes(v) => v.get(row).copied().unwrap_or_default().to_string(),
            Self::Masked(data, mask) => match mask.get(row).copied().unwrap_or(0) {
                1 => ".".to_string(),
                2 => "?".to_string(),
                _ => data.string_at(row),
            },
        }
    }

    pub(crate) fn f32_at(&self, row: usize) -> Option<f32> {
        match self {
            Self::Int(v) => v.get(row).copied().map(|value| value as f32),
            Self::Float(v) => v.get(row).copied(),
            Self::Str(v) => v
                .get(row)
                .map(|value| clean_cif_value(value))
                .and_then(|value| parse_js_number_f32(&value)),
            Self::Bytes(v) => v.get(row).copied().map(|value| value as f32),
            Self::Masked(data, mask) => match mask.get(row).copied().unwrap_or(0) {
                0 => data.f32_at(row),
                _ => None,
            },
        }
    }

    pub(crate) fn i32_at(&self, row: usize) -> Option<i32> {
        match self {
            Self::Int(v) => v.get(row).copied(),
            Self::Float(v) => v.get(row).copied().and_then(float_to_i32),
            Self::Str(v) => v
                .get(row)
                .map(|value| clean_cif_value(value))
                .and_then(|value| value.parse::<i32>().ok()),
            Self::Bytes(v) => v.get(row).copied().map(i32::from),
            Self::Masked(data, mask) => match mask.get(row).copied().unwrap_or(0) {
                0 => data.i32_at(row),
                _ => None,
            },
        }
    }

    pub(crate) fn usize_at(&self, row: usize) -> Option<usize> {
        match self {
            Self::Int(v) => v
                .get(row)
                .copied()
                .and_then(|value| usize::try_from(value).ok()),
            Self::Float(v) => v.get(row).copied().and_then(float_to_usize),
            Self::Str(v) => v
                .get(row)
                .map(|value| clean_cif_value(value))
                .and_then(|value| value.parse::<usize>().ok()),
            Self::Bytes(v) => v.get(row).copied().map(usize::from),
            Self::Masked(data, mask) => match mask.get(row).copied().unwrap_or(0) {
                0 => data.usize_at(row),
                _ => None,
            },
        }
    }

    fn to_i32_vec(&self) -> Vec<i32> {
        match self {
            Self::Int(v) => v.clone(),
            Self::Float(v) => v.iter().map(|x| *x as i32).collect(),
            Self::Str(v) => v.iter().map(|x| x.parse::<i32>().unwrap_or(0)).collect(),
            Self::Bytes(v) => v.iter().map(|x| *x as i32).collect(),
            Self::Masked(data, mask) => (0..data.len())
                .map(|row| match mask.get(row).copied().unwrap_or(0) {
                    0 => data.string_at(row).parse::<i32>().unwrap_or(0),
                    value => value,
                })
                .collect(),
        }
    }

    pub(crate) fn with_mask(self, mask: &[i32]) -> ColumnData {
        ColumnData::Masked(Box::new(self), mask.to_vec())
    }
}

fn float_to_i32(value: f32) -> Option<i32> {
    if value.is_finite()
        && value >= i32::MIN as f32
        && value <= i32::MAX as f32
        && value.fract().abs() <= f32::EPSILON
    {
        Some(value as i32)
    } else {
        None
    }
}

fn float_to_usize(value: f32) -> Option<usize> {
    if value.is_finite() && value >= 0.0 && value.fract().abs() <= f32::EPSILON {
        Some(value as usize)
    } else {
        None
    }
}

pub(super) fn decode_bcif_data(value: &MpValue) -> Result<ColumnData, String> {
    let encodings = value
        .get("encoding")
        .and_then(MpValue::as_array)
        .ok_or("BinaryCIF encoded data is missing encoding")?;
    let bytes = value
        .get("data")
        .and_then(MpValue::as_bin)
        .ok_or("BinaryCIF encoded data is missing bytes")?;
    let mut current = ColumnData::Bytes(bytes.to_vec());
    for encoding in encodings.iter().rev() {
        current = decode_bcif_step(current, encoding)?;
    }
    Ok(current)
}

pub(super) fn decode_bcif_step(data: ColumnData, encoding: &MpValue) -> Result<ColumnData, String> {
    let kind = encoding
        .get("kind")
        .and_then(MpValue::as_str)
        .ok_or("BinaryCIF encoding is missing kind")?;
    match kind {
        "ByteArray" => {
            let ty = encoding.get("type").and_then(MpValue::as_i64).unwrap_or(4);
            let ColumnData::Bytes(bytes) = data else {
                return Ok(data);
            };
            decode_byte_array(&bytes, ty as i32)
        }
        "FixedPoint" => {
            let factor = encoding
                .get("factor")
                .and_then(MpValue::as_f64)
                .unwrap_or(1.0) as f32;
            Ok(ColumnData::Float(
                data.to_i32_vec()
                    .into_iter()
                    .map(|v| v as f32 / factor)
                    .collect(),
            ))
        }
        "IntervalQuantization" => {
            let min = encoding.get("min").and_then(MpValue::as_f64).unwrap_or(0.0) as f32;
            let max = encoding.get("max").and_then(MpValue::as_f64).unwrap_or(0.0) as f32;
            let steps = encoding
                .get("numSteps")
                .and_then(MpValue::as_f64)
                .unwrap_or(1.0) as f32;
            let delta = if steps > 1.0 {
                (max - min) / (steps - 1.0)
            } else {
                0.0
            };
            Ok(ColumnData::Float(
                data.to_i32_vec()
                    .into_iter()
                    .map(|v| min + delta * v as f32)
                    .collect(),
            ))
        }
        "RunLength" => {
            let src_size = encoding
                .get("srcSize")
                .and_then(MpValue::as_i64)
                .unwrap_or(0) as usize;
            let values = data.to_i32_vec();
            let mut out = Vec::with_capacity(src_size);
            for pair in values.chunks(2) {
                if pair.len() == 2 {
                    out.extend(std::iter::repeat_n(pair[0], pair[1].max(0) as usize));
                }
            }
            Ok(ColumnData::Int(out))
        }
        "Delta" => {
            let origin = encoding
                .get("origin")
                .and_then(MpValue::as_i64)
                .unwrap_or(0) as i32;
            let values = data.to_i32_vec();
            let mut out = Vec::with_capacity(values.len());
            let mut acc = origin;
            for value in values {
                acc += value;
                out.push(acc);
            }
            Ok(ColumnData::Int(out))
        }
        "IntegerPacking" => {
            let byte_count = encoding
                .get("byteCount")
                .and_then(MpValue::as_i64)
                .unwrap_or(1);
            let is_unsigned = encoding
                .get("isUnsigned")
                .and_then(MpValue::as_bool)
                .unwrap_or(false);
            let src_size = encoding
                .get("srcSize")
                .and_then(MpValue::as_i64)
                .unwrap_or(0) as usize;
            Ok(ColumnData::Int(integer_unpack(
                &data.to_i32_vec(),
                byte_count,
                is_unsigned,
                src_size,
            )))
        }
        "StringArray" => {
            let string_data = encoding
                .get("stringData")
                .and_then(MpValue::as_str)
                .unwrap_or("");
            let offsets_value = MpValue::Map(vec![
                (
                    "encoding".to_string(),
                    encoding
                        .get("offsetEncoding")
                        .cloned()
                        .unwrap_or_else(|| MpValue::Array(Vec::new())),
                ),
                (
                    "data".to_string(),
                    encoding
                        .get("offsets")
                        .cloned()
                        .unwrap_or_else(|| MpValue::Bin(Vec::new())),
                ),
            ]);
            let offsets = decode_bcif_data(&offsets_value)?.to_i32_vec();
            let data_encoding = MpValue::Map(vec![
                (
                    "encoding".to_string(),
                    encoding
                        .get("dataEncoding")
                        .cloned()
                        .unwrap_or_else(|| MpValue::Array(Vec::new())),
                ),
                (
                    "data".to_string(),
                    match data {
                        ColumnData::Bytes(bytes) => MpValue::Bin(bytes),
                        other => {
                            MpValue::Bin(other.to_i32_vec().into_iter().map(|v| v as u8).collect())
                        }
                    },
                ),
            ]);
            let indices = decode_bcif_data(&data_encoding)?.to_i32_vec();
            let mut strings = vec![String::new()];
            let offsets = offsets
                .into_iter()
                .map(|v| v.max(0) as usize)
                .collect::<Vec<_>>();
            if offsets.first().copied() == Some(0) {
                for pair in offsets.windows(2) {
                    let (start, end) = (pair[0], pair[1]);
                    if start <= end && end <= string_data.len() {
                        strings.push(string_data[start..end].to_string());
                    }
                }
            } else {
                let mut start = 0usize;
                for end in offsets {
                    if end <= string_data.len() && start <= end {
                        strings.push(string_data[start..end].to_string());
                        start = end;
                    }
                }
            }
            Ok(ColumnData::Str(
                indices
                    .into_iter()
                    .map(|idx| {
                        strings
                            .get((idx + 1).max(0) as usize)
                            .cloned()
                            .unwrap_or_default()
                    })
                    .collect(),
            ))
        }
        _ => Err(format!("unsupported BinaryCIF encoding kind: {kind}")),
    }
}

pub(super) fn decode_byte_array(bytes: &[u8], ty: i32) -> Result<ColumnData, String> {
    match ty {
        1 => Ok(ColumnData::Int(
            bytes.iter().map(|b| *b as i8 as i32).collect(),
        )),
        2 => Ok(ColumnData::Int(
            byte_chunks(bytes, 2, ty)?
                .map(|b| i16::from_le_bytes([b[0], b[1]]) as i32)
                .collect(),
        )),
        3 => Ok(ColumnData::Int(
            byte_chunks(bytes, 4, ty)?
                .map(|b| i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
                .collect(),
        )),
        4 => Ok(ColumnData::Int(bytes.iter().map(|b| *b as i32).collect())),
        5 => Ok(ColumnData::Int(
            byte_chunks(bytes, 2, ty)?
                .map(|b| u16::from_le_bytes([b[0], b[1]]) as i32)
                .collect(),
        )),
        6 => Ok(ColumnData::Int(
            byte_chunks(bytes, 4, ty)?
                .map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as i32)
                .collect(),
        )),
        32 => Ok(ColumnData::Float(
            byte_chunks(bytes, 4, ty)?
                .map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
                .collect(),
        )),
        33 => Ok(ColumnData::Float(
            byte_chunks(bytes, 8, ty)?
                .map(|b| {
                    f64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]) as f32
                })
                .collect(),
        )),
        _ => Err(format!("unsupported BinaryCIF byte array type: {ty}")),
    }
}

fn byte_chunks<'a>(
    bytes: &'a [u8],
    width: usize,
    ty: i32,
) -> Result<std::slice::ChunksExact<'a, u8>, String> {
    let chunks = bytes.chunks_exact(width);
    if !chunks.remainder().is_empty() {
        return Err(format!(
            "BinaryCIF byte array type {ty} has {} bytes, not divisible by {width}",
            bytes.len()
        ));
    }
    Ok(chunks)
}

fn integer_unpack(values: &[i32], byte_count: i64, is_unsigned: bool, src_size: usize) -> Vec<i32> {
    if values.len() == src_size {
        return values.to_vec();
    }
    let upper = match (byte_count, is_unsigned) {
        (1, true) => 0xff,
        (2, true) => 0xffff,
        (1, false) => 0x7f,
        _ => 0x7fff,
    };
    let lower = if is_unsigned { 0 } else { -upper - 1 };
    let mut out = Vec::with_capacity(src_size);
    let mut i = 0;
    while i < values.len() {
        let mut value = 0;
        while i < values.len() && (values[i] == upper || (!is_unsigned && values[i] == lower)) {
            value += values[i];
            i += 1;
        }
        if i < values.len() {
            value += values[i];
            i += 1;
        }
        out.push(value);
    }
    out
}

#[derive(Clone, Debug)]
pub(super) enum MpValue {
    Nil,
    Bool(bool),
    Int(i64),
    Float(f64),
    Str(String),
    Bin(Vec<u8>),
    Array(Vec<MpValue>),
    Map(Vec<(String, MpValue)>),
}

impl MpValue {
    fn get(&self, key: &str) -> Option<&MpValue> {
        match self {
            MpValue::Map(items) => items.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    fn as_array(&self) -> Option<&[MpValue]> {
        match self {
            MpValue::Array(v) => Some(v),
            _ => None,
        }
    }

    fn as_str(&self) -> Option<&str> {
        match self {
            MpValue::Str(v) => Some(v),
            _ => None,
        }
    }

    fn as_bin(&self) -> Option<&[u8]> {
        match self {
            MpValue::Bin(v) => Some(v),
            _ => None,
        }
    }

    fn as_i64(&self) -> Option<i64> {
        match self {
            MpValue::Int(v) => Some(*v),
            MpValue::Float(v) => Some(*v as i64),
            _ => None,
        }
    }

    fn as_f64(&self) -> Option<f64> {
        match self {
            MpValue::Int(v) => Some(*v as f64),
            MpValue::Float(v) => Some(*v),
            _ => None,
        }
    }

    fn as_bool(&self) -> Option<bool> {
        match self {
            MpValue::Bool(v) => Some(*v),
            _ => None,
        }
    }
}

pub(super) struct MsgPack<'a> {
    data: &'a [u8],
    offset: usize,
}

impl<'a> MsgPack<'a> {
    pub(super) fn parse(data: &'a [u8]) -> Result<MpValue, String> {
        let mut parser = Self { data, offset: 0 };
        let value = parser.value()?;
        if parser.offset != parser.data.len() {
            return Err("MessagePack has trailing bytes".to_string());
        }
        Ok(value)
    }

    fn value(&mut self) -> Result<MpValue, String> {
        let ty = self.u8()?;
        if (ty & 0x80) == 0x00 {
            return Ok(MpValue::Int(ty as i64));
        }
        if (ty & 0xf0) == 0x80 {
            return self.map((ty & 0x0f) as usize);
        }
        if (ty & 0xf0) == 0x90 {
            return self.array((ty & 0x0f) as usize);
        }
        if (ty & 0xe0) == 0xa0 {
            return self.string((ty & 0x1f) as usize);
        }
        if (ty & 0xe0) == 0xe0 {
            return Ok(MpValue::Int((ty as i8) as i64));
        }
        match ty {
            0xc0 => Ok(MpValue::Nil),
            0xc2 => Ok(MpValue::Bool(false)),
            0xc3 => Ok(MpValue::Bool(true)),
            0xc4 => {
                let len = self.u8()? as usize;
                self.bin(len)
            }
            0xc5 => {
                let len = self.u16()? as usize;
                self.bin(len)
            }
            0xc6 => {
                let len = self.u32()? as usize;
                self.bin(len)
            }
            0xca => Ok(MpValue::Float(f32::from_bits(self.u32()?) as f64)),
            0xcb => Ok(MpValue::Float(f64::from_bits(self.u64()?))),
            0xcc => Ok(MpValue::Int(self.u8()? as i64)),
            0xcd => Ok(MpValue::Int(self.u16()? as i64)),
            0xce => Ok(MpValue::Int(self.u32()? as i64)),
            0xcf => Ok(MpValue::Int(self.u64()? as i64)),
            0xd0 => Ok(MpValue::Int(self.u8()? as i8 as i64)),
            0xd1 => Ok(MpValue::Int(self.u16()? as i16 as i64)),
            0xd2 => Ok(MpValue::Int(self.u32()? as i32 as i64)),
            0xd3 => Ok(MpValue::Int(self.u64()? as i64)),
            0xd9 => {
                let len = self.u8()? as usize;
                self.string(len)
            }
            0xda => {
                let len = self.u16()? as usize;
                self.string(len)
            }
            0xdb => {
                let len = self.u32()? as usize;
                self.string(len)
            }
            0xdc => {
                let len = self.u16()? as usize;
                self.array(len)
            }
            0xdd => {
                let len = self.u32()? as usize;
                self.array(len)
            }
            0xde => {
                let len = self.u16()? as usize;
                self.map(len)
            }
            0xdf => {
                let len = self.u32()? as usize;
                self.map(len)
            }
            _ => Err(format!("unsupported MessagePack type 0x{ty:02x}")),
        }
    }

    fn array(&mut self, len: usize) -> Result<MpValue, String> {
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            out.push(self.value()?);
        }
        Ok(MpValue::Array(out))
    }

    fn map(&mut self, len: usize) -> Result<MpValue, String> {
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            let key = match self.value()? {
                MpValue::Str(s) => s,
                other => format!("{other:?}"),
            };
            let value = self.value()?;
            out.push((key, value));
        }
        Ok(MpValue::Map(out))
    }

    fn string(&mut self, len: usize) -> Result<MpValue, String> {
        let bytes = self.take(len)?;
        Ok(MpValue::Str(
            std::str::from_utf8(bytes)
                .map_err(|_| "MessagePack string is not UTF-8".to_string())?
                .to_string(),
        ))
    }

    fn bin(&mut self, len: usize) -> Result<MpValue, String> {
        Ok(MpValue::Bin(self.take(len)?.to_vec()))
    }

    fn take(&mut self, len: usize) -> Result<&'a [u8], String> {
        if self.offset + len > self.data.len() {
            return Err("MessagePack ended unexpectedly".to_string());
        }
        let out = &self.data[self.offset..self.offset + len];
        self.offset += len;
        Ok(out)
    }

    fn u8(&mut self) -> Result<u8, String> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16, String> {
        let bytes = self.take(2)?;
        Ok(u16::from_be_bytes([bytes[0], bytes[1]]))
    }

    fn u32(&mut self) -> Result<u32, String> {
        let bytes = self.take(4)?;
        Ok(u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
    }

    fn u64(&mut self) -> Result<u64, String> {
        let bytes = self.take(8)?;
        Ok(u64::from_be_bytes([
            bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
        ]))
    }
}
