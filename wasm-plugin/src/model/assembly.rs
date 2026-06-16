use super::Transform;

#[derive(Clone, Debug)]
pub struct SecondaryRange {
    pub chain: String,
    pub start: i32,
    pub start_insertion_code: String,
    pub end: i32,
    pub end_insertion_code: String,
}

#[derive(Clone, Debug)]
pub struct Assembly {
    pub id: String,
    pub details: String,
    pub oligomeric_details: String,
    pub oligomeric_count: Option<i32>,
    pub asym_ids: Vec<String>,
    pub transforms: Vec<Transform>,
    pub generators: Vec<AssemblyGenerator>,
}

#[derive(Clone, Debug)]
pub struct AssemblyOperator {
    pub name: String,
    pub instance_id: String,
    pub assembly_id: String,
    pub oper_id: usize,
    pub oper_list_ids: Vec<String>,
    pub transform: Transform,
}

impl AssemblyOperator {
    pub fn new(
        assembly_id: impl Into<String>,
        oper_id: usize,
        oper_list_ids: Vec<String>,
        transform: Transform,
    ) -> Self {
        let name = format!("ASM_{oper_id}");
        let instance_id = if oper_list_ids.is_empty() {
            name.clone()
        } else {
            format!("ASM-{}", oper_list_ids.join("-"))
        };
        AssemblyOperator {
            name,
            instance_id,
            assembly_id: assembly_id.into(),
            oper_id,
            oper_list_ids,
            transform,
        }
    }
}

#[derive(Clone, Debug)]
pub struct AssemblyGenerator {
    pub asym_ids: Vec<String>,
    pub transforms: Vec<Transform>,
    pub oper_list_ids: Vec<Vec<String>>,
    pub operators: Vec<AssemblyOperator>,
}

impl AssemblyGenerator {
    pub fn from_transforms(
        assembly_id: &str,
        asym_ids: Vec<String>,
        start_oper_id: usize,
        transforms: Vec<Transform>,
        oper_list_ids: Vec<Vec<String>>,
    ) -> Self {
        let operators = transforms
            .iter()
            .enumerate()
            .map(|(index, transform)| {
                AssemblyOperator::new(
                    assembly_id,
                    start_oper_id + index + 1,
                    oper_list_ids.get(index).cloned().unwrap_or_default(),
                    *transform,
                )
            })
            .collect();
        AssemblyGenerator {
            asym_ids,
            transforms,
            oper_list_ids,
            operators,
        }
    }

    pub fn operators_for_assembly(
        &self,
        assembly_id: &str,
        start_oper_id: usize,
    ) -> Vec<AssemblyOperator> {
        if !self.operators.is_empty() {
            return self.operators.clone();
        }
        self.transforms
            .iter()
            .enumerate()
            .map(|(index, transform)| {
                AssemblyOperator::new(
                    assembly_id,
                    start_oper_id + index + 1,
                    self.oper_list_ids.get(index).cloned().unwrap_or_default(),
                    *transform,
                )
            })
            .collect()
    }
}
