//! Mol* post-process antialiasing dispatch.

use super::super::framebuffer::Framebuffer;
use super::super::options::SmaaParams;
use super::smaa;

pub(in crate::render) fn apply_smaa(framebuffer: &mut Framebuffer, params: &SmaaParams) {
    smaa::apply(framebuffer, params);
}
