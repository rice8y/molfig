#![cfg_attr(target_arch = "wasm32", no_main)]

mod api;
mod chemistry;
mod diff;
mod export;
mod json;
mod mesh;
mod model;
mod options;
mod parser;
#[cfg(target_arch = "wasm32")]
mod typst_plugin;

pub use api::{
    convert_to_mtl, convert_to_obj, convert_to_obj_bundle, convert_to_ply,
    convert_to_render_object_bundle, convert_to_stl, maquette_material_map, molecule_info,
    stl_export_facet_context, stl_export_facet_context_timed, stl_facet_semantic_context,
    stl_facet_semantic_context_with_vertex_offset,
    stl_facet_semantic_context_with_vertex_offset_timed,
};
pub use diff::{diff_bytes, diff_bytes_with_options, diff_text, DiffBytesOptions, DiffReport};
pub use model::{
    get_saccharide_name, get_saccharide_shape, saccharide_component, saccharide_component_with_map,
};
pub use model::{
    Assembly, AssemblyGenerator, AssemblyOperator, Atom, AtomDerivedData, AtomSiteAnisotrop,
    AtomicDerivedData, AtomicHierarchy, AtomicIndex, AtomicModel, AtomicRanges, AtomicSegmentation,
    AtomicStructure, AtomicUnitKind, Axes3D, Bond, BondFlags, BondMetadata, BondSource, Boundary,
    BoundingSphere, BranchedEntityLinkMap, BranchedEntityLinkPlacement, BranchedSequenceMap,
    CarbohydrateElement, CarbohydrateLink, CarbohydrateSymbolGeometry, CarbohydrateTerminalLink,
    Carbohydrates, ChemicalComponent, ChemicalComponentAngle, ChemicalComponentAtom,
    ChemicalComponentBond, CoarseConformation, CoarseElement, CoarseElementKey, CoarseElementKind,
    CoarseElementReference, CoarseElements, CoarseGaussian, CoarseGaussianConformation,
    CoarseHierarchy, CoarseIndex, CoarseModel, CoarseRange, CoarseSegmentation, CoarseSphere,
    CoarseSphereConformation, CustomPropertyContainer, DelocalizedTriplets, Entity, EntityIndexMap,
    EntityPoly, EntityPolySeq, Entry, Experiment, GlobalModelTransform, IhmCrossLinkRestraint,
    IhmModelGroup, IhmModelGroupLink, IhmModelList, IndexPairBonds, IntraUnitBondProps,
    IntraUnitBonds, LookupHit, ModelPropertyData, Molecule, MoleculeType, MolstarBondSiteEntry,
    Operator, PartialCarbohydrateElement, PdbxBranchScheme, PdbxEntityBranch, PdbxEntityBranchLink,
    PdbxNonpolyScheme, PdbxPolySeqScheme, PolymerType, PrincipalAxes, ResidueDerivedData,
    Resonance, Ring, SaccharideCompIdMapType, SaccharideComponent, SaccharideShape, SaccharideType,
    SecondaryRange, SequenceEntity, SequenceRange, SequenceResidue, SourceCategory, SourceData,
    StructAsym, StructConnMetadata, StructureLookup3D, StructureLookupHit, StructureProperties,
    StructureSequence, StructureUnit, Transform, UnitKind, UnitLookup3D, UnitOperator, UnitProps,
    UnitSymmetryGroup, UnitTraits, Vec3,
};
pub use options::{InputFormat, MeshOptions};
pub use parser::parse_molecule;

pub const MOLSTAR_REFERENCE_COMMIT: &str = "1b8117d3f10f7c978aabb5a0d3d47370635aefe4";

#[cfg(test)]
mod tests;
