use crate::model::Vec3;

pub(super) struct Framebuffer {
    pub(super) width: usize,
    pub(super) height: usize,
    pub(super) color: Vec<[u8; 4]>,
    pub(super) depth: Vec<f32>,
    pub(super) depth01: Vec<f32>,
    pub(super) normal: Vec<Vec3>,
}

impl Framebuffer {
    pub(super) fn new(width: usize, height: usize, background: [u8; 4]) -> Self {
        let len = width.saturating_mul(height);
        Self {
            width,
            height,
            color: vec![background; len],
            depth: vec![f32::INFINITY; len],
            depth01: vec![1.0; len],
            normal: vec![Vec3::default(); len],
        }
    }
}
