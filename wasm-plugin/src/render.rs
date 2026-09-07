//! Deterministic Mol*-oriented CPU renderer.
//!
//! This facade is intentionally small. Renderer configuration, camera/depth,
//! rasterization, shading, postprocessing, reporting, and image encoding live
//! in separate modules so each subsystem can be checked against the matching
//! pinned Mol* source independently.

mod camera;
mod color;
mod framebuffer;
mod multisample;
mod options;
mod output;
mod pipeline;
mod postprocessing;
mod raster;
mod report;
mod shading;
mod style;
mod transparency;

pub(crate) use options::RendererOptions;
pub(crate) use output::{encode_render_output, RenderOutputFormat};
pub(crate) use pipeline::render_scene;
pub(crate) use report::native_render_objects_json;
