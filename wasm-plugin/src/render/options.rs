use serde::Deserialize;

use super::color::parse_color;

#[derive(Clone, Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(crate) struct RendererOptions {
    pub(super) viewport: ViewportOptions,
    pub(super) camera: CameraOptions,
    pub(super) background: BackgroundOptions,
    pub(super) shading: ShadingOptions,
    pub(super) lighting: LightingOptions,
    pub(super) transparency: TransparencyOptions,
    pub(super) multi_sample: MultiSampleOptions,
    pub(super) postprocessing: PostprocessingOptions,
}

impl RendererOptions {
    pub(crate) fn from_json(bytes: &[u8]) -> Result<Self, String> {
        let mut value: Self = serde_json::from_slice(bytes)
            .map_err(|error| format!("invalid renderer options: {error}"))?;
        value.validate()?;
        value.background.color_value = parse_color(&value.background.color)?;
        Ok(value)
    }

    fn validate(&self) -> Result<(), String> {
        if !(1..=4096).contains(&self.viewport.width) || !(1..=4096).contains(&self.viewport.height)
        {
            return Err("renderer.viewport width and height must be between 1 and 4096".into());
        }
        finite_range64(
            self.viewport.pixel_ratio,
            0.25,
            4.0,
            "renderer.viewport.pixel-ratio",
        )?;
        finite_range64(self.camera.fov, 1.0, 179.0, "renderer.camera.fov")?;
        if !matches!(self.camera.mode.as_str(), "perspective" | "orthographic") {
            return Err(format!(
                "renderer.camera.mode must be \"perspective\" or \"orthographic\"; got {}",
                self.camera.mode
            ));
        }
        self.camera.view.validate()?;
        self.camera.fog.validate("renderer.camera.fog")?;
        finite_range64(
            self.camera.clipping.min_near,
            0.1,
            1000.0,
            "renderer.camera.clipping.min-near",
        )?;
        finite_range64(
            self.camera.clipping.min_far,
            0.0,
            f32::MAX as f64,
            "renderer.camera.clipping.min-far",
        )?;
        parse_color(&self.background.color)?;
        finite_range(
            self.shading.material.metalness,
            0.0,
            1.0,
            "renderer.shading.material.metalness",
        )?;
        finite_range(
            self.shading.material.roughness,
            0.0,
            1.0,
            "renderer.shading.material.roughness",
        )?;
        finite_range(
            self.shading.material.bumpiness,
            0.0,
            1.0,
            "renderer.shading.material.bumpiness",
        )?;
        if self.shading.material.bumpiness != 0.0 {
            return Err(
                "renderer.shading.material.bumpiness is not implemented; use 0 until representation-specific bump frequency is supported"
                    .into(),
            );
        }
        finite_range(
            self.lighting.exposure,
            0.0,
            10.0,
            "renderer.lighting.exposure",
        )?;
        parse_color(&self.lighting.ambient.color)?;
        finite_range64(
            self.lighting.ambient.intensity,
            0.0,
            10.0,
            "renderer.lighting.ambient.intensity",
        )?;
        for (index, light) in self.lighting.directional.iter().enumerate() {
            parse_color(&light.color)?;
            finite_range64(
                light.intensity,
                0.0,
                10.0,
                &format!("renderer.lighting.directional[{index}].intensity"),
            )?;
            finite_range64(
                light.inclination,
                0.0,
                180.0,
                &format!("renderer.lighting.directional[{index}].inclination"),
            )?;
            finite_range64(
                light.azimuth,
                -3600.0,
                3600.0,
                &format!("renderer.lighting.directional[{index}].azimuth"),
            )?;
        }
        if self.transparency.mode != "wboit" {
            return Err(format!(
                "renderer.transparency.mode currently supports only \"wboit\"; got {}",
                self.transparency.mode
            ));
        }
        if !matches!(self.multi_sample.mode.as_str(), "off" | "on" | "temporal") {
            return Err(format!(
                "renderer.multi-sample.mode must be \"off\", \"on\", or \"temporal\"; got {}",
                self.multi_sample.mode
            ));
        }
        if self.multi_sample.sample_level > 5 {
            return Err("renderer.multi-sample.sample-level must be from 0 through 5".into());
        }
        if let Some(pass) = &self.postprocessing.occlusion {
            pass.validate()?;
        }
        if let Some(pass) = &self.postprocessing.outline {
            pass.validate()?;
        }
        if let Some(pass) = &self.postprocessing.shadow {
            pass.validate()?;
        }
        if let Some(antialiasing) = &self.postprocessing.antialiasing {
            antialiasing.validate()?;
        }
        Ok(())
    }
}

fn finite_range(value: f32, min: f32, max: f32, name: &str) -> Result<(), String> {
    if !value.is_finite() || value < min || value > max {
        Err(format!(
            "{name} must be a finite number from {min} through {max}"
        ))
    } else {
        Ok(())
    }
}

fn finite_range64(value: f64, min: f64, max: f64, name: &str) -> Result<(), String> {
    if !value.is_finite() || value < min || value > max {
        Err(format!(
            "{name} must be a finite number from {min} through {max}"
        ))
    } else {
        Ok(())
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct ViewportOptions {
    pub(super) width: u32,
    pub(super) height: u32,
    /// WebGL pixel ratio is a JavaScript Number and participates in drawing
    /// buffer sizing plus outline CPU state before any uniform upload.
    pub(super) pixel_ratio: f64,
}

impl Default for ViewportOptions {
    fn default() -> Self {
        Self {
            width: 800,
            height: 800,
            pixel_ratio: 1.0,
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct CameraOptions {
    pub(super) mode: String,
    pub(super) fov: f64,
    pub(super) view: ViewOptions,
    pub(super) fog: NamedFog,
    pub(super) clipping: ClippingOptions,
}

impl Default for CameraOptions {
    fn default() -> Self {
        Self {
            mode: "perspective".into(),
            fov: 45.0,
            view: ViewOptions::default(),
            fog: NamedFog::default(),
            clipping: ClippingOptions::default(),
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct ViewOptions {
    pub(super) name: String,
    pub(super) params: ViewParams,
}

impl Default for ViewOptions {
    fn default() -> Self {
        Self {
            name: "auto".into(),
            params: ViewParams::default(),
        }
    }
}

impl ViewOptions {
    fn validate(&self) -> Result<(), String> {
        match self.name.as_str() {
            "auto" => {
                if !self.params.is_empty() {
                    return Err(
                        "renderer.camera.view.params must be empty when renderer.camera.view.name is \"auto\""
                            .into(),
                    );
                }
            }
            "orbit" => {
                if self.params.has_snapshot_values() {
                    return Err(
                        "renderer.camera.view.params position, target, up, radius, and radius-max are valid only for the \"snapshot\" view"
                            .into(),
                    );
                }
                finite_range64(
                    self.params.azimuth.unwrap_or(0.0),
                    -3600.0,
                    3600.0,
                    "renderer.camera.view.params.azimuth",
                )?;
                finite_range64(
                    self.params.elevation.unwrap_or(0.0),
                    -89.9,
                    89.9,
                    "renderer.camera.view.params.elevation",
                )?;
                finite_range64(
                    self.params.roll.unwrap_or(0.0),
                    -3600.0,
                    3600.0,
                    "renderer.camera.view.params.roll",
                )?;
            }
            "snapshot" => {
                if self.params.has_orbit_values() {
                    return Err(
                        "renderer.camera.view.params azimuth, elevation, and roll are valid only for the \"orbit\" view"
                            .into(),
                    );
                }
                let position = self.params.position.ok_or_else(|| {
                    "renderer.camera.view.params.position is required for the \"snapshot\" view"
                        .to_string()
                })?;
                let target = self.params.target.ok_or_else(|| {
                    "renderer.camera.view.params.target is required for the \"snapshot\" view"
                        .to_string()
                })?;
                let up = self.params.up.ok_or_else(|| {
                    "renderer.camera.view.params.up is required for the \"snapshot\" view"
                        .to_string()
                })?;
                for (name, value) in [
                    ("position", position),
                    ("target", target),
                    ("up", up),
                ] {
                    if value.iter().any(|v| !v.is_finite()) {
                        return Err(format!("renderer.camera.view.params.{name} must contain finite numbers"));
                    }
                }
                finite_range64(
                    self.params.radius.ok_or_else(|| {
                        "renderer.camera.view.params.radius is required for the \"snapshot\" view"
                            .to_string()
                    })?,
                    0.01,
                    f32::MAX as f64,
                    "renderer.camera.view.params.radius",
                )?;
                finite_range64(
                    self.params.radius_max.ok_or_else(|| {
                        "renderer.camera.view.params.radius-max is required for the \"snapshot\" view"
                            .to_string()
                    })?,
                    0.01,
                    f32::MAX as f64,
                    "renderer.camera.view.params.radius-max",
                )?;
            }
            other => {
                return Err(format!(
                    "renderer.camera.view.name must be \"auto\", \"orbit\", or \"snapshot\"; got {other}"
                ))
            }
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct ViewParams {
    pub(super) azimuth: Option<f64>,
    pub(super) elevation: Option<f64>,
    pub(super) roll: Option<f64>,
    pub(super) position: Option<[f64; 3]>,
    pub(super) target: Option<[f64; 3]>,
    pub(super) up: Option<[f64; 3]>,
    pub(super) radius: Option<f64>,
    pub(super) radius_max: Option<f64>,
}

impl ViewParams {
    fn is_empty(&self) -> bool {
        !self.has_orbit_values() && !self.has_snapshot_values()
    }

    fn has_orbit_values(&self) -> bool {
        self.azimuth.is_some() || self.elevation.is_some() || self.roll.is_some()
    }

    fn has_snapshot_values(&self) -> bool {
        self.position.is_some()
            || self.target.is_some()
            || self.up.is_some()
            || self.radius.is_some()
            || self.radius_max.is_some()
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct NamedFog {
    pub(super) name: String,
    pub(super) params: FogParams,
}

impl Default for NamedFog {
    fn default() -> Self {
        Self {
            name: "on".into(),
            params: FogParams::default(),
        }
    }
}

impl NamedFog {
    fn validate(&self, field: &str) -> Result<(), String> {
        if !matches!(self.name.as_str(), "on" | "off") {
            return Err(format!("{field}.name must be \"on\" or \"off\""));
        }
        finite_range64(
            self.params.intensity,
            0.0,
            100.0,
            &format!("{field}.params.intensity"),
        )
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct FogParams {
    pub(super) intensity: f64,
}

impl Default for FogParams {
    fn default() -> Self {
        Self { intensity: 15.0 }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct ClippingOptions {
    pub(super) far: bool,
    pub(super) min_near: f64,
    pub(super) min_far: f64,
}

impl Default for ClippingOptions {
    fn default() -> Self {
        Self {
            far: true,
            min_near: 1.0,
            min_far: 0.0,
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct BackgroundOptions {
    pub(super) color: String,
    pub(super) transparent: bool,
    /// Parsed once at the renderer API boundary. Opaque fragment fogging runs
    /// before the RGBA8 attachment conversion and must not repeatedly parse
    /// the public string for every covered fragment.
    #[serde(skip)]
    pub(super) color_value: u32,
}

impl Default for BackgroundOptions {
    fn default() -> Self {
        Self {
            color: "#fcfbfa".into(),
            transparent: false,
            color_value: 0xfcfbfa,
        }
    }
}

#[derive(Clone, Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct ShadingOptions {
    pub(super) ignore_light: Option<bool>,
    pub(super) material: MaterialOptions,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct MaterialOptions {
    pub(super) metalness: f32,
    pub(super) roughness: f32,
    pub(super) bumpiness: f32,
}

impl Default for MaterialOptions {
    fn default() -> Self {
        Self {
            metalness: 0.0,
            roughness: 1.0,
            bumpiness: 0.0,
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct LightingOptions {
    pub(super) exposure: f32,
    pub(super) ambient: AmbientLight,
    pub(super) directional: Vec<DirectionalLight>,
}

impl Default for LightingOptions {
    fn default() -> Self {
        Self {
            exposure: 1.0,
            ambient: AmbientLight::default(),
            directional: vec![DirectionalLight::default()],
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct AmbientLight {
    pub(super) color: String,
    pub(super) intensity: f64,
}

impl Default for AmbientLight {
    fn default() -> Self {
        Self {
            color: "#ffffff".into(),
            intensity: 0.4,
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct DirectionalLight {
    pub(super) inclination: f64,
    pub(super) azimuth: f64,
    pub(super) color: String,
    pub(super) intensity: f64,
}

impl Default for DirectionalLight {
    fn default() -> Self {
        Self {
            inclination: 150.0,
            azimuth: 320.0,
            color: "#ffffff".into(),
            intensity: 0.6,
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct TransparencyOptions {
    pub(super) mode: String,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct MultiSampleOptions {
    pub(super) mode: String,
    pub(super) sample_level: usize,
    pub(super) reuse_occlusion: bool,
}

impl Default for MultiSampleOptions {
    fn default() -> Self {
        Self {
            mode: "temporal".into(),
            sample_level: 2,
            reuse_occlusion: true,
        }
    }
}

impl Default for TransparencyOptions {
    fn default() -> Self {
        Self {
            mode: "wboit".into(),
        }
    }
}

#[derive(Clone, Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct PostprocessingOptions {
    pub(super) occlusion: Option<OcclusionPass>,
    pub(super) outline: Option<OutlinePass>,
    pub(super) shadow: Option<ShadowPass>,
    pub(super) antialiasing: Option<AntialiasingOptions>,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct AntialiasingOptions {
    pub(super) name: String,
    pub(super) params: SmaaParams,
}

impl Default for AntialiasingOptions {
    fn default() -> Self {
        Self {
            name: "smaa".into(),
            params: SmaaParams::default(),
        }
    }
}

impl AntialiasingOptions {
    fn validate(&self) -> Result<(), String> {
        if !matches!(self.name.as_str(), "smaa" | "off") {
            return Err(format!(
                "renderer.postprocessing.antialiasing.name must be \"smaa\" or \"off\"; got {}",
                self.name
            ));
        }
        finite_range(
            self.params.edge_threshold,
            0.05,
            0.15,
            "renderer.postprocessing.antialiasing.params.edge-threshold",
        )?;
        if self.params.max_search_steps > 32 {
            return Err(
                "renderer.postprocessing.antialiasing.params.max-search-steps must be from 0 through 32"
                    .into(),
            );
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct SmaaParams {
    pub(super) edge_threshold: f32,
    pub(super) max_search_steps: usize,
}

impl Default for SmaaParams {
    fn default() -> Self {
        Self {
            edge_threshold: 0.1,
            max_search_steps: 16,
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct OcclusionPass {
    pub(super) name: String,
    pub(super) params: OcclusionParams,
}

impl Default for OcclusionPass {
    fn default() -> Self {
        Self {
            name: "off".into(),
            params: OcclusionParams::default(),
        }
    }
}

impl OcclusionPass {
    fn validate(&self) -> Result<(), String> {
        if !matches!(self.name.as_str(), "on" | "off") {
            return Err("renderer.postprocessing.occlusion.name must be \"on\" or \"off\"".into());
        }
        parse_color(&self.params.color)?;
        if !(1..=256).contains(&self.params.samples) {
            return Err(
                "renderer.postprocessing.occlusion.params.samples must be from 1 through 256"
                    .into(),
            );
        }
        self.params.multi_scale.validate()?;
        finite_range64(
            self.params.radius,
            0.0,
            20.0,
            "renderer.postprocessing.occlusion.params.radius",
        )?;
        finite_range64(
            self.params.bias,
            0.0,
            3.0,
            "renderer.postprocessing.occlusion.params.bias",
        )?;
        if !(1..=25).contains(&self.params.blur_kernel_size)
            || self.params.blur_kernel_size.is_multiple_of(2)
        {
            return Err(
                "renderer.postprocessing.occlusion.params.blur-kernel-size must be an odd integer from 1 through 25"
                    .into(),
            );
        }
        finite_range64(
            self.params.blur_depth_bias,
            0.0,
            1.0,
            "renderer.postprocessing.occlusion.params.blur-depth-bias",
        )?;
        finite_range64(
            self.params.resolution_scale,
            0.1,
            1.0,
            "renderer.postprocessing.occlusion.params.resolution-scale",
        )?;
        finite_range64(
            self.params.transparent_threshold,
            0.0,
            1.0,
            "renderer.postprocessing.occlusion.params.transparent-threshold",
        )
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct OcclusionParams {
    pub(super) samples: usize,
    pub(super) multi_scale: OcclusionMultiScale,
    /// Mol* keeps these values as JavaScript Numbers until uploading the
    /// corresponding WebGL uniforms. Retaining f64 here prevents fractional
    /// authoring values from being rounded before `Math.pow(2, radius)` and
    /// the camera-scale multiplication.
    pub(super) radius: f64,
    pub(super) bias: f64,
    pub(super) blur_kernel_size: usize,
    pub(super) blur_depth_bias: f64,
    pub(super) resolution_scale: f64,
    pub(super) color: String,
    /// Keep this as a JavaScript-number-width value because Mol* compares it
    /// against `1 - alpha` before uploading anything to WebGL. Converting the
    /// default 0.4 to f32 first would incorrectly make alpha 0.6 pass Mol*'s
    /// strict `< transparentThreshold` test.
    pub(super) transparent_threshold: f64,
}

impl Default for OcclusionParams {
    fn default() -> Self {
        Self {
            samples: 32,
            multi_scale: OcclusionMultiScale::default(),
            radius: 5.0,
            bias: 0.8,
            blur_kernel_size: 15,
            blur_depth_bias: 0.5,
            resolution_scale: 1.0,
            color: "#000000".into(),
            transparent_threshold: 0.4,
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct OcclusionMultiScale {
    pub(super) name: String,
    pub(super) params: OcclusionMultiScaleParams,
}

impl Default for OcclusionMultiScale {
    fn default() -> Self {
        Self {
            name: "off".into(),
            params: OcclusionMultiScaleParams::default(),
        }
    }
}

impl OcclusionMultiScale {
    fn validate(&self) -> Result<(), String> {
        if !matches!(self.name.as_str(), "on" | "off") {
            return Err(
                "renderer.postprocessing.occlusion.params.multi-scale.name must be \"on\" or \"off\""
                    .into(),
            );
        }
        for (index, level) in self.params.levels.iter().enumerate() {
            finite_range64(
                level.radius,
                0.0,
                20.0,
                &format!(
                    "renderer.postprocessing.occlusion.params.multi-scale.params.levels[{index}].radius"
                ),
            )?;
            finite_range64(
                level.bias,
                0.0,
                3.0,
                &format!(
                    "renderer.postprocessing.occlusion.params.multi-scale.params.levels[{index}].bias"
                ),
            )?;
        }
        finite_range64(
            self.params.near_threshold,
            0.0,
            50.0,
            "renderer.postprocessing.occlusion.params.multi-scale.params.near-threshold",
        )?;
        finite_range64(
            self.params.far_threshold,
            0.0,
            10_000.0,
            "renderer.postprocessing.occlusion.params.multi-scale.params.far-threshold",
        )
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct OcclusionMultiScaleParams {
    pub(super) levels: Vec<OcclusionLevel>,
    pub(super) near_threshold: f64,
    pub(super) far_threshold: f64,
}

impl Default for OcclusionMultiScaleParams {
    fn default() -> Self {
        Self {
            levels: [2.0, 5.0, 8.0, 11.0]
                .map(|radius| OcclusionLevel { radius, bias: 1.0 })
                .to_vec(),
            near_threshold: 10.0,
            far_threshold: 1_500.0,
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct OcclusionLevel {
    pub(super) radius: f64,
    pub(super) bias: f64,
}

impl Default for OcclusionLevel {
    fn default() -> Self {
        Self {
            radius: 5.0,
            bias: 1.0,
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct OutlinePass {
    pub(super) name: String,
    pub(super) params: OutlineParams,
}

impl Default for OutlinePass {
    fn default() -> Self {
        Self {
            name: "off".into(),
            params: OutlineParams::default(),
        }
    }
}

impl OutlinePass {
    fn validate(&self) -> Result<(), String> {
        if !matches!(self.name.as_str(), "on" | "off") {
            return Err("renderer.postprocessing.outline.name must be \"on\" or \"off\"".into());
        }
        parse_color(&self.params.color)?;
        finite_range64(
            self.params.scale,
            1.0,
            5.0,
            "renderer.postprocessing.outline.params.scale",
        )?;
        finite_range64(
            self.params.threshold,
            0.01,
            1.0,
            "renderer.postprocessing.outline.params.threshold",
        )
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct OutlineParams {
    pub(super) color: String,
    /// Both values remain JavaScript Numbers through Mol*'s pixel-ratio
    /// scaling. The threshold converts to float32 only at `uniform1f`.
    pub(super) scale: f64,
    pub(super) threshold: f64,
    pub(super) include_transparent: bool,
}

impl Default for OutlineParams {
    fn default() -> Self {
        Self {
            color: "#000000".into(),
            scale: 1.0,
            threshold: 0.33,
            include_transparent: true,
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub(super) struct ShadowPass {
    pub(super) name: String,
}

impl Default for ShadowPass {
    fn default() -> Self {
        Self { name: "off".into() }
    }
}

impl ShadowPass {
    fn validate(&self) -> Result<(), String> {
        match self.name.as_str() {
            "off" => Ok(()),
            "on" => {
                Err("renderer.postprocessing.shadow is not implemented; use name: \"off\"".into())
            }
            _ => Err("renderer.postprocessing.shadow.name must be \"on\" or \"off\"".into()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renderer_options_reject_unknown_fields() {
        let error = RendererOptions::from_json(br#"{"unknown":true}"#).unwrap_err();
        assert!(error.contains("unknown field"));
    }

    #[test]
    fn renderer_options_validate_colors_and_modes() {
        RendererOptions::from_json(br##"{"background":{"color":"#ffffff","transparent":true}}"##)
            .unwrap();
        assert!(RendererOptions::from_json(br#"{"camera":{"mode":"fish-eye"}}"#).is_err());
        assert!(RendererOptions::from_json(br#"{"background":{"color":"white"}}"#).is_err());
        assert!(RendererOptions::from_json(br#"{"multi-sample":{"sample-level":6}}"#).is_err());
        assert!(
            RendererOptions::from_json(br#"{"multi-sample":{"reduce-flicker":false}}"#)
                .unwrap_err()
                .contains("unknown field")
        );
        RendererOptions::from_json(
            br#"{"postprocessing":{"antialiasing":{"name":"smaa","params":{"edge-threshold":0.1,"max-search-steps":16}}}}"#,
        )
        .unwrap();
        assert!(RendererOptions::from_json(
            br#"{"postprocessing":{"antialiasing":{"name":"smaa","params":{"max-search-steps":33}}}}"#,
        )
        .is_err());
    }

    #[test]
    fn occlusion_multi_scale_uses_the_molstar_mapped_parameters() {
        let options = RendererOptions::from_json(
            br#"{"postprocessing":{"occlusion":{"name":"on","params":{"multi-scale":{"name":"on","params":{"levels":[{"radius":2,"bias":1},{"radius":8.5,"bias":0.75}],"near-threshold":10,"far-threshold":1500}}}}}}"#,
        )
        .unwrap();
        let multi_scale = &options
            .postprocessing
            .occlusion
            .as_ref()
            .unwrap()
            .params
            .multi_scale;
        assert_eq!(multi_scale.name, "on");
        assert_eq!(multi_scale.params.levels.len(), 2);
        assert_eq!(multi_scale.params.levels[1].radius, 8.5);
        assert_eq!(multi_scale.params.levels[1].bias, 0.75);

        assert!(RendererOptions::from_json(
            br#"{"postprocessing":{"occlusion":{"name":"on","params":{"multi-scale":{"name":"on","params":{"levels":[{"radius":21,"bias":1}]}}}}}}"#,
        )
        .unwrap_err()
        .contains("levels[0].radius"));
        assert!(RendererOptions::from_json(
            br#"{"postprocessing":{"occlusion":{"name":"on","params":{"multi-scale":{"name":"adaptive"}}}}}"#,
        )
        .unwrap_err()
        .contains("multi-scale.name"));
    }

    #[test]
    fn camera_view_schema_rejects_params_from_other_mapped_variants() {
        RendererOptions::from_json(br#"{"camera":{"view":{"name":"auto","params":{}}}}"#).unwrap();
        RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"orbit","params":{"elevation":24}}}}"#,
        )
        .unwrap();
        RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"snapshot","params":{"position":[0,0,10],"target":[0,0,0],"up":[0,1,0],"radius":2,"radius-max":3}}}}"#,
        )
        .unwrap();

        assert!(RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"auto","params":{"azimuth":1}}}}"#,
        )
        .unwrap_err()
        .contains("must be empty"));
        assert!(RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"orbit","params":{"position":[0,0,10]}}}}"#,
        )
        .unwrap_err()
        .contains("valid only for the \"snapshot\" view"));
        assert!(RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"snapshot","params":{"azimuth":1,"position":[0,0,10],"target":[0,0,0],"up":[0,1,0],"radius":2,"radius-max":3}}}}"#,
        )
        .unwrap_err()
        .contains("valid only for the \"orbit\" view"));
        assert!(RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"snapshot","params":{"position":[0,0,10]}}}}"#,
        )
        .unwrap_err()
        .contains("target is required"));
    }

    #[test]
    fn camera_clipping_preserves_snapshot_fields_and_canvas3d_ranges() {
        let options = RendererOptions::from_json(
            br#"{"camera":{"clipping":{"far":false,"min-near":0.3,"min-far":12.75}}}"#,
        )
        .unwrap();
        assert!(!options.camera.clipping.far);
        assert_eq!(options.camera.clipping.min_near, 0.3);
        assert_eq!(options.camera.clipping.min_far, 12.75);

        assert!(
            RendererOptions::from_json(br#"{"camera":{"clipping":{"min-near":0.099}}}"#,)
                .unwrap_err()
                .contains("0.1 through 1000")
        );
        assert!(
            RendererOptions::from_json(br#"{"camera":{"clipping":{"min-far":-0.01}}}"#,)
                .unwrap_err()
                .contains("0 through")
        );
    }

    #[test]
    fn pass_schemas_reject_unrelated_and_unimplemented_values() {
        RendererOptions::from_json(
            br##"{"postprocessing":{"occlusion":{"name":"on","params":{"samples":16,"transparent-threshold":0.4}}}}"##,
        )
        .unwrap();
        assert!(RendererOptions::from_json(
            br##"{"postprocessing":{"occlusion":{"name":"on","params":{"scale":2}}}}"##,
        )
        .unwrap_err()
        .contains("unknown field"));
        assert!(RendererOptions::from_json(
            br##"{"postprocessing":{"outline":{"name":"on","params":{"samples":16}}}}"##,
        )
        .unwrap_err()
        .contains("unknown field"));
        assert!(RendererOptions::from_json(
            br##"{"postprocessing":{"shadow":{"name":"off","params":{}}}}"##,
        )
        .unwrap_err()
        .contains("unknown field"));
        let scaled = RendererOptions::from_json(
            br##"{"postprocessing":{"occlusion":{"name":"on","params":{"resolution-scale":0.5}}}}"##,
        )
        .unwrap();
        assert_eq!(
            scaled
                .postprocessing
                .occlusion
                .unwrap()
                .params
                .resolution_scale,
            0.5
        );
        assert!(
            RendererOptions::from_json(br##"{"shading":{"material":{"bumpiness":0.25}}}"##,)
                .unwrap_err()
                .contains("not implemented")
        );
    }
}
