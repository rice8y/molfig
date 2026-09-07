pub(super) mod angle_metal_math;
mod antialiasing;
mod composition;
mod fog;
mod outline;
mod packing;
mod smaa;
mod ssao;
mod ssao_sinc;

pub(super) use antialiasing::apply_smaa;
pub(super) use composition::{read_color, write_composed_color};
pub(super) use fog::{fog_factor, fog_opaque_fragment, fog_range};
pub(super) use outline::{apply_outline, apply_transparent_outline};
#[cfg(test)]
pub(super) use packing::packed_depth_alpha_roundtrip;
pub(in crate::render) use packing::{pack_depth_alpha_rgba8, unpack_depth_alpha_rgba};
pub(super) use ssao::{
    apply_occlusion_factors, compute_occlusion_factors,
    compute_occlusion_factors_including_transparency, ssao_target_dimensions, OcclusionFactors,
};
