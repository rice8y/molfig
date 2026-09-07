/// Render PDB, mmCIF, BinaryCIF, and XYZ molecular structures in Typst with a
/// deterministic Mol*-oriented native renderer, and export interchange meshes.
///
/// Public functions accept structure bytes, inline PDB, mmCIF, or XYZ text, and
/// Typst 0.15+ path values. File inputs on older Typst versions must be read
/// with `read(..., encoding: none)` before being passed to Molfig.
#let _plugin = plugin("molfig.wasm")

#let _v15-or-later() = {
  sys.version >= version(0, 15, 0)
}

#let _is-path(value) = {
  _v15-or-later() and str(type(value)) == "path"
}

#let _normalize-data(data) = {
  if _is-path(data) {
    read(data, encoding: none)
  } else if type(data) == bytes {
    data
  } else if type(data) == str {
    bytes(data)
  } else {
    panic("molfig expects bytes, inline string data, or a Typst 0.15+ path. Use read(\"molecule.pdb\", encoding: none) for Typst 0.14-compatible files, or path(\"molecule.pdb\") on Typst 0.15 or later.")
  }
}

#let _is-number(value) = type(value) == int or type(value) == float

#let _u32-le(data, offset) = {
  data.at(offset) + data.at(offset + 1) * 256 + data.at(offset + 2) * 65536 + data.at(offset + 3) * 16777216
}

#let _native-render-bundle(bundle) = {
  if bundle.len() < 32 or str(bundle.slice(0, 4)) != "MFRG" {
    panic("invalid Molfig native render bundle")
  }
  let version = _u32-le(bundle, 4)
  if version != 2 {
    panic("unsupported Molfig native render bundle version: " + str(version))
  }
  let pixel-width = _u32-le(bundle, 8)
  let pixel-height = _u32-le(bundle, 12)
  let rgba-len = _u32-le(bundle, 16)
  let image-len = _u32-le(bundle, 20)
  let info-len = _u32-le(bundle, 24)
  let render-info-len = _u32-le(bundle, 28)
  let pixels-start = 32
  let image-start = pixels-start + rgba-len
  let info-start = image-start + image-len
  let render-info-start = info-start + info-len
  let end = render-info-start + render-info-len
  if rgba-len != pixel-width * pixel-height * 4 or end != bundle.len() {
    panic("invalid Molfig native render bundle lengths")
  }
  (
    pixels: bundle.slice(pixels-start, image-start),
    image: bundle.slice(image-start, info-start),
    pixel-width: pixel-width,
    pixel-height: pixel-height,
    info: json(bundle.slice(info-start, render-info-start)),
    render-info: json(bundle.slice(render-info-start, end)),
  )
}

#let _renderer-options(renderer) = {
  if type(renderer) != dictionary {
    panic("renderer must be a dictionary")
  }
  json.encode(renderer)
}

#let _output-format(value) = {
  if value != "svg" and value != "png" {
    panic("output-format must be one of \"svg\" or \"png\"")
  }
  value
}

#let _mesh-options(format, representation, color-theme, style, style-params, theme, sphere-detail, radius-scale, atom-radius, bond-radius, infer-bonds, center, assembly, alt-loc, block-index, block-header, ribbon-radius, ribbon-width, helix-profile, round-cap, sheet-arrow-factor, tubular-helices, linear-segments, radial-segments, quality, decimate, mesh-format: none) = {
  let options = (
    format: format,
    representation: representation,
    color-theme: color-theme,
    style: style,
    style-params: style-params,
    theme: theme,
    sphere-detail: sphere-detail,
    radius-scale: radius-scale,
    atom-radius: atom-radius,
    bond-radius: bond-radius,
    infer-bonds: infer-bonds,
    center: center,
    assembly: assembly,
    alt-loc: alt-loc,
    block-index: block-index,
    block-header: block-header,
    ribbon-radius: ribbon-radius,
    ribbon-width: ribbon-width,
    helix-profile: helix-profile,
    round-cap: round-cap,
    sheet-arrow-factor: sheet-arrow-factor,
    tubular-helices: tubular-helices,
    linear-segments: linear-segments,
    radial-segments: radial-segments,
    quality: quality,
    decimate: decimate,
  )
  if mesh-format != none {
    options += (mesh-format: mesh-format)
  }
  json.encode(options)
}

/// Export a molecular structure as a Wavefront OBJ mesh.
///
/// OBJ preserves Molfig face groups, operator metadata when requested, and
/// material identifiers for color themes.
///
/// - data (any): Structure bytes, inline PDB, mmCIF, or XYZ text, or a Typst 0.15+ path value.
/// - format (str): Input format: `"auto"`, `"pdb"`, `"cif"`, `"mmcif"`, `"bcif"`, or `"xyz"`.
/// - representation (str): Molecular representation, such as `"default"`, `"cartoon"`, `"spacefill"`, `"ball-and-stick"`, `"surface"`, `"ribbon"`, or `"backbone"`.
/// - color-theme (str): Mol\* color theme used to assign OBJ materials.
/// - style (str): Rendering style: `"default"` or `"illustrative"`; independent of representation and color theme.
/// - style-params (dictionary): Illustrative parameters: `ignore-light`, `outline`, and `occlusion`.
/// - theme (dictionary): Mol\* Viewer theme overrides, including `globalName`, `carbonColor`, and `symmetryColor`.
/// - sphere-detail (int): Icosphere subdivision detail used by sphere-based visuals.
/// - radius-scale (int, float): Global multiplier for molecular radii.
/// - atom-radius (int, float): Base atom radius used by ball-and-stick-style visuals.
/// - bond-radius (int, float): Base bond-cylinder radius.
/// - infer-bonds (bool): Whether missing covalent bonds are inferred from molecular geometry.
/// - center (bool): Whether the exported mesh is centered using the visible Mol\* bounding sphere.
/// - assembly (str): Biological assembly identifier, or `"asymmetric-unit"` to use the source asymmetric unit.
/// - alt-loc (str): Alternate-location selector; an empty string follows Mol* and includes all locations.
/// - block-index (none, int): Zero-based CIF or BinaryCIF data-block index.
/// - block-header (str): CIF or BinaryCIF data-block header to select instead of `block-index`.
/// - ribbon-radius (int, float): Polymer tube radius for ribbon-derived geometry.
/// - ribbon-width (int, float): Polymer ribbon width.
/// - helix-profile (str): Helix cross-section profile: `"elliptical"`, `"rounded"`, or `"square"`.
/// - round-cap (bool): Whether polymer segment ends use rounded caps.
/// - sheet-arrow-factor (int, float): Width multiplier for beta-sheet arrow tips.
/// - tubular-helices (bool): Whether helices are emitted as tubes instead of ribbons.
/// - linear-segments (int): Longitudinal curve subdivisions.
/// - radial-segments (int): Radial profile subdivisions.
/// - quality (str): Geometry quality preset from `"lowest"` through `"highest"`, or `"auto"`/`"custom"`.
/// - decimate (int, float): Molfig semantic decimation strength in the range `0` to `1`.
/// -> bytes
#let to-obj(
  data,
  format: "auto",
  representation: "default",
  color-theme: "chain-id",
  style: "default",
  style-params: (:),
  theme: (:),
  sphere-detail: 2,
  radius-scale: 1.0,
  atom-radius: 0.28,
  bond-radius: 0.12,
  infer-bonds: true,
  center: true,
  assembly: "1",
  alt-loc: "",
  block-index: none,
  block-header: "",
  ribbon-radius: 0.2,
  ribbon-width: 0.55,
  helix-profile: "elliptical",
  round-cap: false,
  sheet-arrow-factor: 1.5,
  tubular-helices: false,
  linear-segments: 8,
  radial-segments: 16,
  quality: "custom",
  decimate: 0,
) = _plugin.to_obj(_normalize-data(data), bytes(_mesh-options(format, representation, color-theme, style, style-params, theme, sphere-detail, radius-scale, atom-radius, bond-radius, infer-bonds, center, assembly, alt-loc, block-index, block-header, ribbon-radius, ribbon-width, helix-profile, round-cap, sheet-arrow-factor, tubular-helices, linear-segments, radial-segments, quality, decimate)))

/// Export the OBJ materials for a molecular structure as a Wavefront MTL file.
///
/// Material order and identifiers match the corresponding `to-obj` output.
///
/// - data (any): Structure bytes, inline PDB, mmCIF, or XYZ text, or a Typst 0.15+ path value.
/// - format (str): Input format: `"auto"`, `"pdb"`, `"cif"`, `"mmcif"`, `"bcif"`, or `"xyz"`.
/// - representation (str): Molecular representation used to generate material assignments.
/// - color-theme (str): Mol\* color theme used to assign materials.
/// - style (str): Rendering style: `"default"` or `"illustrative"`; independent of representation and color theme.
/// - style-params (dictionary): Illustrative parameters: `ignore-light`, `outline`, and `occlusion`.
/// - theme (dictionary): Mol\* Viewer theme overrides, including `globalName`, `carbonColor`, and `symmetryColor`.
/// - sphere-detail (int): Icosphere subdivision detail used by sphere-based visuals.
/// - radius-scale (int, float): Global multiplier for molecular radii.
/// - atom-radius (int, float): Base atom radius used by ball-and-stick-style visuals.
/// - bond-radius (int, float): Base bond-cylinder radius.
/// - infer-bonds (bool): Whether missing covalent bonds are inferred from molecular geometry.
/// - center (bool): Whether the associated mesh is centered.
/// - assembly (str): Biological assembly identifier, or `"asymmetric-unit"` to use the source asymmetric unit.
/// - alt-loc (str): Alternate-location selector; an empty string follows Mol* and includes all locations.
/// - block-index (none, int): Zero-based CIF or BinaryCIF data-block index.
/// - block-header (str): CIF or BinaryCIF data-block header to select instead of `block-index`.
/// - ribbon-radius (int, float): Polymer tube radius for ribbon-derived geometry.
/// - ribbon-width (int, float): Polymer ribbon width.
/// - helix-profile (str): Helix cross-section profile: `"elliptical"`, `"rounded"`, or `"square"`.
/// - round-cap (bool): Whether polymer segment ends use rounded caps.
/// - sheet-arrow-factor (int, float): Width multiplier for beta-sheet arrow tips.
/// - tubular-helices (bool): Whether helices are emitted as tubes instead of ribbons.
/// - linear-segments (int): Longitudinal curve subdivisions.
/// - radial-segments (int): Radial profile subdivisions.
/// - quality (str): Geometry quality preset from `"lowest"` through `"highest"`, or `"auto"`/`"custom"`.
/// - decimate (int, float): Molfig semantic decimation strength in the range `0` to `1`.
/// -> bytes
#let to-mtl(
  data,
  format: "auto",
  representation: "default",
  color-theme: "chain-id",
  style: "default",
  style-params: (:),
  theme: (:),
  sphere-detail: 2,
  radius-scale: 1.0,
  atom-radius: 0.28,
  bond-radius: 0.12,
  infer-bonds: true,
  center: true,
  assembly: "1",
  alt-loc: "",
  block-index: none,
  block-header: "",
  ribbon-radius: 0.2,
  ribbon-width: 0.55,
  helix-profile: "elliptical",
  round-cap: false,
  sheet-arrow-factor: 1.5,
  tubular-helices: false,
  linear-segments: 8,
  radial-segments: 16,
  quality: "custom",
  decimate: 0,
) = _plugin.to_mtl(_normalize-data(data), bytes(_mesh-options(format, representation, color-theme, style, style-params, theme, sphere-detail, radius-scale, atom-radius, bond-radius, infer-bonds, center, assembly, alt-loc, block-index, block-header, ribbon-radius, ribbon-width, helix-profile, round-cap, sheet-arrow-factor, tubular-helices, linear-segments, radial-segments, quality, decimate)))

/// Export a molecular structure as a binary STL triangle mesh.
///
/// STL contains geometry only; color themes and face groups cannot be
/// represented by the format.
///
/// - data (any): Structure bytes, inline PDB, mmCIF, or XYZ text, or a Typst 0.15+ path value.
/// - format (str): Input format: `"auto"`, `"pdb"`, `"cif"`, `"mmcif"`, `"bcif"`, or `"xyz"`.
/// - representation (str): Molecular representation used to construct the mesh.
/// - color-theme (str): Color theme used during semantic construction; STL does not store its colors.
/// - style (str): Rendering style used during semantic construction; STL does not store its colors.
/// - style-params (dictionary): Illustrative parameters: `ignore-light`, `outline`, and `occlusion`.
/// - theme (dictionary): Mol\* Viewer theme overrides used during semantic construction.
/// - sphere-detail (int): Icosphere subdivision detail used by sphere-based visuals.
/// - radius-scale (int, float): Global multiplier for molecular radii.
/// - atom-radius (int, float): Base atom radius used by ball-and-stick-style visuals.
/// - bond-radius (int, float): Base bond-cylinder radius.
/// - infer-bonds (bool): Whether missing covalent bonds are inferred from molecular geometry.
/// - center (bool): Whether the exported mesh is centered using the visible Mol\* bounding sphere.
/// - assembly (str): Biological assembly identifier, or `"asymmetric-unit"` to use the source asymmetric unit.
/// - alt-loc (str): Alternate-location selector; an empty string follows Mol* and includes all locations.
/// - block-index (none, int): Zero-based CIF or BinaryCIF data-block index.
/// - block-header (str): CIF or BinaryCIF data-block header to select instead of `block-index`.
/// - ribbon-radius (int, float): Polymer tube radius for ribbon-derived geometry.
/// - ribbon-width (int, float): Polymer ribbon width.
/// - helix-profile (str): Helix cross-section profile: `"elliptical"`, `"rounded"`, or `"square"`.
/// - round-cap (bool): Whether polymer segment ends use rounded caps.
/// - sheet-arrow-factor (int, float): Width multiplier for beta-sheet arrow tips.
/// - tubular-helices (bool): Whether helices are emitted as tubes instead of ribbons.
/// - linear-segments (int): Longitudinal curve subdivisions.
/// - radial-segments (int): Radial profile subdivisions.
/// - quality (str): Geometry quality preset from `"lowest"` through `"highest"`, or `"auto"`/`"custom"`.
/// - decimate (int, float): Molfig semantic decimation strength in the range `0` to `1`.
/// -> bytes
#let to-stl(
  data,
  format: "auto",
  representation: "default",
  color-theme: "chain-id",
  style: "default",
  style-params: (:),
  theme: (:),
  sphere-detail: 2,
  radius-scale: 1.0,
  atom-radius: 0.28,
  bond-radius: 0.12,
  infer-bonds: true,
  center: true,
  assembly: "1",
  alt-loc: "",
  block-index: none,
  block-header: "",
  ribbon-radius: 0.2,
  ribbon-width: 0.55,
  helix-profile: "elliptical",
  round-cap: false,
  sheet-arrow-factor: 1.5,
  tubular-helices: false,
  linear-segments: 8,
  radial-segments: 16,
  quality: "custom",
  decimate: 0,
) = _plugin.to_stl(_normalize-data(data), bytes(_mesh-options(format, representation, color-theme, style, style-params, theme, sphere-detail, radius-scale, atom-radius, bond-radius, infer-bonds, center, assembly, alt-loc, block-index, block-header, ribbon-radius, ribbon-width, helix-profile, round-cap, sheet-arrow-factor, tubular-helices, linear-segments, radial-segments, quality, decimate)))

/// Export a molecular structure as an ASCII PLY triangle mesh.
///
/// PLY preserves Molfig face-group values but does not carry OBJ materials.
///
/// - data (any): Structure bytes, inline PDB, mmCIF, or XYZ text, or a Typst 0.15+ path value.
/// - format (str): Input format: `"auto"`, `"pdb"`, `"cif"`, `"mmcif"`, `"bcif"`, or `"xyz"`.
/// - representation (str): Molecular representation used to construct the mesh.
/// - color-theme (str): Color theme used during semantic construction; PLY does not store its colors.
/// - style (str): Rendering style used during semantic construction; PLY does not store its colors.
/// - style-params (dictionary): Illustrative parameters: `ignore-light`, `outline`, and `occlusion`.
/// - theme (dictionary): Mol\* Viewer theme overrides used during semantic construction.
/// - sphere-detail (int): Icosphere subdivision detail used by sphere-based visuals.
/// - radius-scale (int, float): Global multiplier for molecular radii.
/// - atom-radius (int, float): Base atom radius used by ball-and-stick-style visuals.
/// - bond-radius (int, float): Base bond-cylinder radius.
/// - infer-bonds (bool): Whether missing covalent bonds are inferred from molecular geometry.
/// - center (bool): Whether the exported mesh is centered using the visible Mol\* bounding sphere.
/// - assembly (str): Biological assembly identifier, or `"asymmetric-unit"` to use the source asymmetric unit.
/// - alt-loc (str): Alternate-location selector; an empty string follows Mol* and includes all locations.
/// - block-index (none, int): Zero-based CIF or BinaryCIF data-block index.
/// - block-header (str): CIF or BinaryCIF data-block header to select instead of `block-index`.
/// - ribbon-radius (int, float): Polymer tube radius for ribbon-derived geometry.
/// - ribbon-width (int, float): Polymer ribbon width.
/// - helix-profile (str): Helix cross-section profile: `"elliptical"`, `"rounded"`, or `"square"`.
/// - round-cap (bool): Whether polymer segment ends use rounded caps.
/// - sheet-arrow-factor (int, float): Width multiplier for beta-sheet arrow tips.
/// - tubular-helices (bool): Whether helices are emitted as tubes instead of ribbons.
/// - linear-segments (int): Longitudinal curve subdivisions.
/// - radial-segments (int): Radial profile subdivisions.
/// - quality (str): Geometry quality preset from `"lowest"` through `"highest"`, or `"auto"`/`"custom"`.
/// - decimate (int, float): Molfig semantic decimation strength in the range `0` to `1`.
/// -> bytes
#let to-ply(
  data,
  format: "auto",
  representation: "default",
  color-theme: "chain-id",
  style: "default",
  style-params: (:),
  theme: (:),
  sphere-detail: 2,
  radius-scale: 1.0,
  atom-radius: 0.28,
  bond-radius: 0.12,
  infer-bonds: true,
  center: true,
  assembly: "1",
  alt-loc: "",
  block-index: none,
  block-header: "",
  ribbon-radius: 0.2,
  ribbon-width: 0.55,
  helix-profile: "elliptical",
  round-cap: false,
  sheet-arrow-factor: 1.5,
  tubular-helices: false,
  linear-segments: 8,
  radial-segments: 16,
  quality: "custom",
  decimate: 0,
) = _plugin.to_ply(_normalize-data(data), bytes(_mesh-options(format, representation, color-theme, style, style-params, theme, sphere-detail, radius-scale, atom-radius, bond-radius, infer-bonds, center, assembly, alt-loc, block-index, block-header, ribbon-radius, ribbon-width, helix-profile, round-cap, sheet-arrow-factor, tubular-helices, linear-segments, radial-segments, quality, decimate)))

/// Inspect molecular, Model/Structure/Unit, representation, and mesh-planning metadata.
///
/// This performs Molfig's parsing and semantic representation work without
/// running the native rasterizer.
///
/// - data (any): Structure bytes, inline PDB, mmCIF, or XYZ text, or a Typst 0.15+ path value.
/// - format (str): Input format: `"auto"`, `"pdb"`, `"cif"`, `"mmcif"`, `"bcif"`, or `"xyz"`.
/// - representation (str): Molecular representation whose semantic metadata is inspected.
/// - color-theme (str): Mol\* color theme used during semantic construction.
/// - style (str): Rendering style: `"default"` or `"illustrative"`; independent of representation and color theme.
/// - style-params (dictionary): Illustrative parameters: `ignore-light`, `outline`, and `occlusion`.
/// - theme (dictionary): Mol\* Viewer theme overrides, including `globalName`, `carbonColor`, and `symmetryColor`.
/// - sphere-detail (int): Icosphere subdivision detail used by sphere-based visuals.
/// - radius-scale (int, float): Global multiplier for molecular radii.
/// - atom-radius (int, float): Base atom radius used by ball-and-stick-style visuals.
/// - bond-radius (int, float): Base bond-cylinder radius.
/// - infer-bonds (bool): Whether missing covalent bonds are inferred from molecular geometry.
/// - center (bool): Whether geometry bounds and export planning use centered coordinates.
/// - assembly (str): Biological assembly identifier, or `"asymmetric-unit"` to use the source asymmetric unit.
/// - alt-loc (str): Alternate-location selector; an empty string follows Mol* and includes all locations.
/// - block-index (none, int): Zero-based CIF or BinaryCIF data-block index.
/// - block-header (str): CIF or BinaryCIF data-block header to select instead of `block-index`.
/// - ribbon-radius (int, float): Polymer tube radius for ribbon-derived geometry.
/// - ribbon-width (int, float): Polymer ribbon width.
/// - helix-profile (str): Helix cross-section profile: `"elliptical"`, `"rounded"`, or `"square"`.
/// - round-cap (bool): Whether polymer segment ends use rounded caps.
/// - sheet-arrow-factor (int, float): Width multiplier for beta-sheet arrow tips.
/// - tubular-helices (bool): Whether helices are emitted as tubes instead of ribbons.
/// - linear-segments (int): Longitudinal curve subdivisions.
/// - radial-segments (int): Radial profile subdivisions.
/// - quality (str): Geometry quality preset from `"lowest"` through `"highest"`, or `"auto"`/`"custom"`.
/// - decimate (int, float): Molfig semantic decimation strength in the range `0` to `1`.
/// -> dictionary
#let info(
  data,
  format: "auto",
  representation: "default",
  color-theme: "chain-id",
  style: "default",
  style-params: (:),
  theme: (:),
  sphere-detail: 2,
  radius-scale: 1.0,
  atom-radius: 0.28,
  bond-radius: 0.12,
  infer-bonds: true,
  center: true,
  assembly: "1",
  alt-loc: "",
  block-index: none,
  block-header: "",
  ribbon-radius: 0.2,
  ribbon-width: 0.55,
  helix-profile: "elliptical",
  round-cap: false,
  sheet-arrow-factor: 1.5,
  tubular-helices: false,
  linear-segments: 8,
  radial-segments: 16,
  quality: "custom",
  decimate: 0,
) = json(_plugin.info(_normalize-data(data), bytes(_mesh-options(format, representation, color-theme, style, style-params, theme, sphere-detail, radius-scale, atom-radius, bond-radius, infer-bonds, center, assembly, alt-loc, block-index, block-header, ribbon-radius, ribbon-width, helix-profile, round-cap, sheet-arrow-factor, tubular-helices, linear-segments, radial-segments, quality, decimate))))

/// Render a molecular structure and return pixels plus resolved metadata.
///
/// `renderer` is a Mol*-specific dictionary containing `viewport`, `camera`,
/// `background`, `shading`, `lighting`, `transparency`, `multi-sample`, and
/// `postprocessing`.
/// Unknown renderer keys are errors. Style is independent of representation
/// and color theme; explicit renderer values override the selected style.
/// `output-format` is `"svg"` or `"png"` and defaults to `"svg"`.
///
/// -> dictionary
#let render-result(
  data,
  format: "auto",
  representation: "default",
  color-theme: "chain-id",
  style: "default",
  theme: (:),
  sphere-detail: 2,
  radius-scale: 1.0,
  atom-radius: 0.28,
  bond-radius: 0.12,
  infer-bonds: true,
  assembly: "1",
  alt-loc: "",
  block-index: none,
  block-header: "",
  ribbon-radius: 0.2,
  ribbon-width: 0.55,
  helix-profile: "elliptical",
  round-cap: false,
  sheet-arrow-factor: 1.5,
  tubular-helices: false,
  linear-segments: 8,
  radial-segments: 16,
  quality: "custom",
  decimate: 0,
  renderer: (:),
  output-format: "svg",
  width: auto,
  height: auto,
) = {
  if not _is-number(decimate) {
    panic("decimate must be numeric")
  }
  let options = _mesh-options(format, representation, color-theme, style, (:), theme, sphere-detail, radius-scale, atom-radius, bond-radius, infer-bonds, false, assembly, alt-loc, block-index, block-header, ribbon-radius, ribbon-width, helix-profile, round-cap, sheet-arrow-factor, tubular-helices, linear-segments, radial-segments, quality, decimate)
  let source = _normalize-data(data)
  let resolved-output-format = _output-format(output-format)
  let result = _native-render-bundle(_plugin.render_result(source, bytes(options), bytes(_renderer-options(renderer)), bytes(resolved-output-format)))
  (
    kind: "render-result",
    pixels: result.pixels,
    image: result.image,
    output-format: resolved-output-format,
    pixel-width: result.pixel-width,
    pixel-height: result.pixel-height,
    info: result.info,
    render-info: result.render-info,
    content: image(
      result.image,
      width: width,
      height: height,
    ),
  )
}

/// Render a molecular structure as SVG or PNG Typst content.
///
/// The rendering path does not serialize an intermediate OBJ, STL, or PLY.
/// `output-format` defaults to `"svg"`; SVG embeds the exact lossless PNG
/// raster because Mol* screen-space effects do not have a vector equivalent.
/// `width` and `height` affect Typst layout; pixel dimensions are controlled by
/// `renderer.viewport` and default to 800 by 800.
///
/// -> content
#let render(..args) = render-result(..args).content
