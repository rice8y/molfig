use crate::options::{MeshOptions, RenderStyle};

use super::options::RendererOptions;

#[derive(Clone, Copy)]
pub(super) struct ResolvedStyle {
    pub(super) ignore_light: bool,
    pub(super) occlusion: bool,
    pub(super) outline: bool,
}

pub(super) fn resolve_style(options: &MeshOptions, renderer: &RendererOptions) -> ResolvedStyle {
    let illustrative = options.style == RenderStyle::Illustrative;
    ResolvedStyle {
        ignore_light: renderer.shading.ignore_light.unwrap_or(if illustrative {
            options.illustrative_ignore_light
        } else {
            false
        }),
        occlusion: renderer.postprocessing.occlusion.as_ref().map_or(
            if illustrative {
                options.illustrative_occlusion
            } else {
                true
            },
            |pass| pass.name == "on",
        ),
        outline: renderer.postprocessing.outline.as_ref().map_or(
            if illustrative {
                options.illustrative_outline
            } else {
                false
            },
            |pass| pass.name == "on",
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const QUICK_STYLES_SOURCE: &str =
        include_str!("../../artifacts/molstar/src/mol-plugin-ui/structure/quick-styles.tsx");

    #[test]
    fn style_defaults_and_overrides_match_pinned_quick_style_state() {
        for expected in [
            "ignoreLight: false",
            "ignoreLight: true",
            "scale: 1",
            "threshold: 0.33",
            "includeTransparent: true",
            "radius: 5",
            "bias: 0.8",
            "blurKernelSize: 15",
            "blurDepthBias: 0.5",
            "samples: 32",
            "resolutionScale: 1",
            "transparentThreshold: 0.4",
            "shadow: { name: 'off', params: {} }",
        ] {
            assert!(
                QUICK_STYLES_SOURCE.contains(expected),
                "pinned Quick Style source no longer contains {expected:?}"
            );
        }
        assert!(
            !QUICK_STYLES_SOURCE.contains("antialiasing:"),
            "Quick Style unexpectedly started mutating antialiasing state"
        );

        let default_geometry = MeshOptions::from_json(br#"{"style":"default"}"#).unwrap();
        let illustrative_geometry = MeshOptions::from_json(br#"{"style":"illustrative"}"#).unwrap();
        let renderer = RendererOptions::default();
        let default = resolve_style(&default_geometry, &renderer);
        assert!(!default.ignore_light);
        assert!(default.occlusion);
        assert!(!default.outline);
        let illustrative = resolve_style(&illustrative_geometry, &renderer);
        assert!(illustrative.ignore_light);
        assert!(illustrative.occlusion);
        assert!(illustrative.outline);

        let overrides = RendererOptions::from_json(
            br#"{"shading":{"ignore-light":false},"postprocessing":{"occlusion":{"name":"off"},"outline":{"name":"off"}}}"#,
        )
        .unwrap();
        let illustrative = resolve_style(&illustrative_geometry, &overrides);
        assert!(!illustrative.ignore_light);
        assert!(!illustrative.occlusion);
        assert!(!illustrative.outline);

        let antialiasing_off =
            RendererOptions::from_json(br#"{"postprocessing":{"antialiasing":{"name":"off"}}}"#)
                .unwrap();
        for geometry in [&default_geometry, &illustrative_geometry] {
            let _ = resolve_style(geometry, &antialiasing_off);
            assert_eq!(
                antialiasing_off
                    .postprocessing
                    .antialiasing
                    .as_ref()
                    .expect("explicit antialiasing")
                    .name,
                "off"
            );
        }
    }

    #[test]
    fn illustrative_resolution_is_representation_independent() {
        let renderer = RendererOptions::default();
        for representation in [
            "default",
            "cartoon",
            "polymer-cartoon",
            "spacefill",
            "ball-and-stick",
            "surface",
            "ribbon",
            "backbone",
        ] {
            let json = format!(r#"{{"representation":"{representation}","style":"illustrative"}}"#);
            let options = MeshOptions::from_json(json.as_bytes()).unwrap();
            let style = resolve_style(&options, &renderer);
            assert!(style.ignore_light, "{representation}");
            assert!(style.occlusion, "{representation}");
            assert!(style.outline, "{representation}");
        }
    }
}
