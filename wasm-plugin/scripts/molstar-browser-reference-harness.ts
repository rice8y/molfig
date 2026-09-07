import '../artifacts/molstar/src/mol-util/polyfill';
import { ObjExporter } from '../artifacts/molstar/src/extensions/geo-export/obj-exporter';
import { StlExporter } from '../artifacts/molstar/src/extensions/geo-export/stl-exporter';
import { Box3D } from '../artifacts/molstar/src/mol-math/geometry';
import { BoundaryHelper } from '../artifacts/molstar/src/mol-math/geometry/boundary-helper';
import { Structure } from '../artifacts/molstar/src/mol-model/structure';
import { PluginContext } from '../artifacts/molstar/src/mol-plugin/context';
import { PluginConfig } from '../artifacts/molstar/src/mol-plugin/config';
import { DefaultPluginSpec, PluginSpec } from '../artifacts/molstar/src/mol-plugin/spec';
import { OpenFiles } from '../artifacts/molstar/src/mol-plugin-state/actions/file';
import { PresetStructureRepresentations } from '../artifacts/molstar/src/mol-plugin-state/builder/structure/representation-preset';
import { PresetTrajectoryHierarchy } from '../artifacts/molstar/src/mol-plugin-state/builder/structure/hierarchy-preset';
import { ViewerAutoPreset } from '../artifacts/molstar/src/apps/viewer/presets';
import { MAQualityAssessment } from '../artifacts/molstar/src/extensions/model-archive/quality-assessment/behavior';
import { SbNcbrPartialCharges } from '../artifacts/molstar/src/extensions/sb-ncbr';
import { Asset } from '../artifacts/molstar/src/mol-util/assets';
import { Task } from '../artifacts/molstar/src/mol-task';
import type { ValueCell } from '../artifacts/molstar/src/mol-util/value-cell';

type ExportFormat = 'obj' | 'stl' | 'report'

type ViewerTheme = {
    globalName?: string
    globalColorParams?: Record<string, unknown>
    carbonColor?: 'chain-id' | 'operator-name' | 'element-symbol'
    symmetryColor?: string
    symmetryColorParams?: Record<string, unknown>
}

type GaussianSurfaceParams = {
    quality?: string
    resolution?: number
    smoothness?: number
    radiusOffset?: number
    traceOnly?: boolean
    tryUseGpu?: boolean
    visuals?: string[]
}

type MolecularSurfaceParams = {
    quality?: string
    resolution?: number
    probeRadius?: number
    probePositions?: number
    floodfill?: string
    visuals?: string[]
}

type BrowserReferenceRequest = {
    id: string
    fixtureUrl: string
    fileName: string
    formats: ExportFormat[]
    objExportBasename: string
    dataFormat?: string
    structurePreset?: string
    representation?: string
    theme?: ViewerTheme
    sizeThresholds?: Partial<Structure.SizeThresholds>
    gaussianSurfaceParams?: GaussianSurfaceParams
    molecularSurfaceParams?: MolecularSurfaceParams
    renderObjectReport?: boolean
    debugStlFacets?: number[]
    exporterOptions?: {
        primitivesQuality?: 'auto' | 'high' | 'medium' | 'low'
    }
}

type RenderObjectSummary = {
    index: number
    type: string
    geometry: string
    visible: boolean
    drawCount: number
    vertexCount: number
    groupCount: number
    instanceCount: number
    component?: string
    tag?: string
    representation?: string
    colorTheme?: string
    carbonColorTheme?: string
    representationOrder?: number
    visuals?: string[]
    surfaceParams?: GaussianSurfaceParams
    stlTriangleCount?: number
    primitiveCount?: number
    cylinderCapHistogram?: Record<string, number>
    cylinderSamples?: Array<{
        group: number
        cap: number
        scale: number
        start: number[]
        end: number[]
    }>
    meshVertexSamples?: number[][]
    boundingSphere?: {
        center: number[]
        radius: number
        extremaCount: number
        extrema?: number[][]
    }
}

type BrowserReferenceResult = {
    id: string
    renderObjects: RenderObjectSummary[]
    structures: ReturnType<typeof summarizeStructures>
    sceneBoundingSphere: ReturnType<typeof sphereSummary>
    totalDrawCount: number
    visibleDrawCount: number
    hiddenDrawCount: number
    exportableDrawCount: number
    uploads: Record<string, { byteLength: number }>
    webgl: {
        webgl2: boolean
        fragDepth: boolean
        textureFloat: boolean
    }
    debug?: {
        stlFacets: StlFacetDebug[]
    }
}

type StlFacetDebug = {
    stlFacet: number
    found: boolean
    renderObjectIndex?: number
    renderObjectType?: string
    geometry?: string
    visible?: boolean
    localFacet?: number
    localTriangle?: number
    instanceIndex?: number
    objectSlotStart?: number
    objectSlotEnd?: number
    sparseSlot?: number
    vertexIndices?: number[]
    rawVertices?: number[][]
    centeredVertices?: number[][]
    centeredVertexBits?: string[][]
    triangleNormal?: number[]
    triangleNormalBits?: string[]
    centerOffset?: number[]
    boxMin?: number[]
    boxMax?: number[]
    boxExtremaCount?: number
    boxMinIndices?: (number | null)[]
    boxMaxIndices?: (number | null)[]
    boxMinPoints?: (number[] | null)[]
    boxMaxPoints?: (number[] | null)[]
    error?: string
}

declare global {
    interface Window {
        molfigBrowserReferenceExport: (request: BrowserReferenceRequest) => Promise<BrowserReferenceResult>
    }
}

let plugin: PluginContext | undefined;

async function getPlugin() {
    if (plugin) {
        await plugin.clear(true);
        return plugin;
    }

    const target = document.getElementById('app');
    if (!target) throw new Error('missing #app mount target');

    const spec = DefaultPluginSpec();
    spec.behaviors.push(
        PluginSpec.Behavior(MAQualityAssessment),
        PluginSpec.Behavior(SbNcbrPartialCharges),
    );
    spec.config = [
        ...(spec.config ?? []),
        [PluginConfig.Structure.DefaultRepresentationPreset, ViewerAutoPreset.id],
    ];
    plugin = new PluginContext(spec);
    await plugin.init();
    plugin.builders.structure.representation.registerPreset(ViewerAutoPreset);
    const mounted = await plugin.mountAsync(target, { checkeredCanvasBackground: false });
    if (!mounted) throw new Error('failed to mount Mol* plugin');
    await plugin.canvas3dInitialized;
    return plugin;
}

window.molfigBrowserReferenceExport = async (request: BrowserReferenceRequest) => {
    const plugin = await getPlugin();
    const sizeThresholds = { ...Structure.DefaultSizeThresholds, ...(request.sizeThresholds ?? {}) };
    plugin.config.set(PluginConfig.Structure.SizeThresholds, sizeThresholds);
    const bytes = new Uint8Array(await (await fetch(request.fixtureUrl)).arrayBuffer());
    const file = new File([bytes], request.fileName);
    await plugin.runTask(plugin.state.data.applyAction(OpenFiles, {
        files: [Asset.File(file)],
        format: openFilesFormat(request.dataFormat),
        visuals: !request.structurePreset,
    }));
    await applyStructurePreset(plugin, request.structurePreset);
    await applyRequestedRepresentation(plugin, request.representation, request.theme, !request.structurePreset);
    await updateGaussianSurfaceParams(plugin, request.gaussianSurfaceParams);
    await updateMolecularSurfaceParams(plugin, request.molecularSurfaceParams);
    plugin.canvas3d?.commit(true);
    await animationFrames(3);

    const renderObjects = plugin.canvas3d?.getRenderObjects() ?? [];
    if (renderObjects.length === 0) throw new Error('Mol* produced no render objects');
    const summary = summarizeRenderObjects(renderObjects, renderObjectStateMetadata(plugin));
    const sphere = visibleRenderObjectBoundingSphere(renderObjects);
    const box = Box3D.fromSphere3D(Box3D(), sphere);
    const debug = request.debugStlFacets?.length
        ? { stlFacets: request.debugStlFacets.map(facet => debugStlFacet(renderObjects, box, sphere, facet)) }
        : undefined;
    const uploads: BrowserReferenceResult['uploads'] = {};
    for (const format of request.formats) {
        if (format === 'obj') {
            const exporter = new ObjExporter(request.objExportBasename, box);
            configureExporter(exporter, request.exporterOptions);
            const data = await exportWith(plugin, renderObjects, exporter);
            await uploadText(request.id, 'obj', data.obj);
            await uploadText(request.id, 'mtl', data.mtl);
            uploads.obj = { byteLength: new TextEncoder().encode(data.obj).byteLength };
            uploads.mtl = { byteLength: new TextEncoder().encode(data.mtl).byteLength };
        } else if (format === 'stl') {
            const exporter = new StlExporter(box);
            configureExporter(exporter, request.exporterOptions);
            const data = await exportWith(plugin, renderObjects, exporter, (index, triangleCount) => {
                summary.renderObjects[index].stlTriangleCount = triangleCount;
            });
            await uploadBytes(request.id, 'stl', data.stl);
            uploads.stl = { byteLength: data.stl.byteLength };
        }
    }

    return {
        id: request.id,
        ...summary,
        structures: summarizeStructures(plugin, sizeThresholds),
        sceneBoundingSphere: sphereSummary(sphere),
        uploads,
        webgl: webglSummary(plugin),
        debug,
    };
};

function visibleRenderObjectBoundingSphere(renderObjects: any[]) {
    const helper = new BoundaryHelper('98');
    const spheres = renderObjects
        .filter(renderObject => (renderObject.state?.visible ?? true) && Number(value(renderObject.values.drawCount)) > 0)
        .map(renderObject => value(renderObject.values.boundingSphere))
        .filter(sphere => sphere && sphere.radius > 0);
    helper.reset();
    for (const sphere of spheres) helper.includeSphere(sphere);
    helper.finishedIncludeStep();
    for (const sphere of spheres) helper.radiusSphere(sphere);
    return helper.getSphere();
}

function summarizeStructures(plugin: PluginContext, sizeThresholds: Structure.SizeThresholds = Structure.DefaultSizeThresholds) {
    return plugin.managers.structure.hierarchy.selection.structures.map((entry: any, index: number) => {
        const structure = entry?.cell?.obj?.data;
        const operators = new Map<string, any>();
        for (const unit of structure?.units ?? []) {
            const operator = unit?.conformation?.operator;
            if (!operator) continue;
            const key = `${operator.name}|${operator.spgrOp}|${operator.assembly?.id ?? ''}`;
            if (!operators.has(key)) {
                operators.set(key, {
                    name: operator.name,
                    spgrOp: operator.spgrOp,
                    assemblyId: operator.assembly?.id,
                    isAssembly: !!operator.assembly,
                });
            }
        }
        return {
            index,
            unitCount: structure?.units?.length ?? 0,
            elementCount: structure?.elementCount ?? 0,
            polymerResidueCount: structure?.polymerResidueCount ?? 0,
            sizeClass: structure ? Structure.Size[Structure.getSize(structure, sizeThresholds)] : undefined,
            sizeThresholds,
            operators: [...operators.values()],
        };
    });
}

async function updateGaussianSurfaceParams(plugin: PluginContext, requested: GaussianSurfaceParams | undefined) {
    if (!requested) return;
    const update = plugin.state.data.build();
    let count = 0;

    for (const structure of plugin.managers.structure.hierarchy.selection.structures) {
        for (const component of structure.components) {
            for (const representation of component.representations) {
                const old = representation.cell.transform.params as any;
                if (old?.type?.name !== 'gaussian-surface') continue;
                update.to(representation.cell).update({
                    ...old,
                    type: {
                        ...old.type,
                        params: { ...old.type.params, ...requested },
                    },
                });
                count += 1;
            }
        }
    }
    if (count === 0) throw new Error('requested Gaussian surface parameters but ViewerAuto produced no gaussian-surface representation');
    await update.commit({ revertOnError: true });
}

async function updateMolecularSurfaceParams(plugin: PluginContext, requested: MolecularSurfaceParams | undefined) {
    if (!requested) return;
    const update = plugin.state.data.build();
    let count = 0;

    for (const structure of plugin.managers.structure.hierarchy.selection.structures) {
        for (const component of structure.components) {
            for (const representation of component.representations) {
                const old = representation.cell.transform.params as any;
                if (old?.type?.name !== 'molecular-surface') continue;
                update.to(representation.cell).update({
                    ...old,
                    type: {
                        ...old.type,
                        params: { ...old.type.params, ...requested },
                    },
                });
                count += 1;
            }
        }
    }
    if (count === 0) throw new Error('requested molecular surface parameters but the Surface preset produced no molecular-surface representation');
    await update.commit({ revertOnError: true });
}

async function applyStructurePreset(plugin: PluginContext, value: string | undefined) {
    if (!value) return;
    const trajectories = plugin.managers.structure.hierarchy.selection.trajectories;
    const preset = String(value).trim().toLowerCase().replaceAll('_', '-');
    if (preset === 'crystal-contacts') {
        await plugin.managers.structure.hierarchy.applyPreset(trajectories, PresetTrajectoryHierarchy.crystalContacts, {
            representationPreset: 'empty',
        });
        return;
    }
    throw new Error(`unsupported browser structure preset: ${value}`);
}

async function applyRequestedRepresentation(
    plugin: PluginContext,
    value: string | undefined,
    theme: ViewerTheme | undefined,
    canReuseInitial: boolean,
) {
    if (canReuseInitial && representationMatches(plugin, value)) {
        if (theme) await updateViewerThemes(plugin, theme);
        return;
    }
    await applyRepresentationPreset(plugin, value, theme);
    // Built-in presets do not consistently consume the optional theme
    // parameter (notably molecular-surface), so update the realized
    // representations after the preset has created them as well.
    if (theme) await updateViewerThemes(plugin, theme);
}

function representationMatches(plugin: PluginContext, value: string | undefined) {
    const representation = String(value ?? 'default').trim().toLowerCase().replaceAll('_', '-');
    if (representation === 'default' || representation === 'auto') return true;

    const entries = representationStateEntries(plugin);
    if (representation === 'cartoon' || representation === 'polymer-cartoon') {
        return entries.some(entry => entry.type === 'cartoon');
    }
    if (representation === 'spacefill') {
        return entries.some(entry => entry.type === 'spacefill' && entry.tags.includes('all'));
    }
    if (representation === 'surface' || representation === 'molecular-surface') {
        return entries.some(entry => entry.type === 'molecular-surface' && entry.tags.includes('all'));
    }
    return false;
}

function representationStateEntries(plugin: PluginContext) {
    const entries: Array<{ type?: string, tags: string[] }> = [];
    for (const cell of plugin.state.data.cells.values()) {
        const repr = (cell.obj as any)?.data?.repr;
        if (!repr?.renderObjects?.length) continue;
        entries.push({
            type: (cell.transform.params as any)?.type?.name,
            tags: [...(cell.transform.tags ?? [])],
        });
    }
    return entries;
}

async function updateViewerThemes(plugin: PluginContext, theme: ViewerTheme) {
    await plugin.dataTransaction(async () => {
        for (const structure of plugin.managers.structure.hierarchy.selection.structures) {
            await plugin.managers.structure.component.updateRepresentationsTheme(
                structure.components,
                (component: any, representation: any) => viewerThemeForRepresentation(component, representation, theme),
            );
        }
    }, { canUndo: 'Update Viewer Theme' });
}

function viewerThemeForRepresentation(component: any, representation: any, theme: ViewerTheme) {
    const representationParams = representation.cell.transform.params as any;
    const representationType = representationParams?.type?.name;
    const representationTags = [...(representation.cell.transform.tags ?? [])];
    const componentTags = [...(component.cell.transform.tags ?? [])];
    const componentName = componentTags
        .find((tag: string) => tag.startsWith('structure-component-static-'))
        ?.slice('structure-component-static-'.length);
    const globalName = theme.globalName || representationParams?.colorTheme?.name;
    const globalParams = { ...(theme.globalColorParams ?? {}) };

    if (representationTags.includes('polymer') && theme.symmetryColor && structureHasSymmetry(component.structure.cell.obj?.data)) {
        return {
            color: theme.symmetryColor as any,
            colorParams: { ...globalParams, ...(theme.symmetryColorParams ?? {}) },
        };
    }

    if (representationType === 'ball-and-stick') {
        const carbonColor = componentName === 'water' || componentName === 'ion' || componentName === 'lipid'
            ? 'element-symbol'
            : theme.carbonColor;
        return {
            color: globalName as any,
            colorParams: carbonColor
                ? { carbonColor: { name: carbonColor, params: {} }, ...globalParams }
                : globalParams,
        };
    }

    return { color: globalName as any, colorParams: globalParams };
}

function structureHasSymmetry(structure: any) {
    return (structure?.units ?? []).some((unit: any) => {
        const operator = unit?.conformation?.operator;
        return operator && !operator.assembly && operator.spgrOp >= 0;
    });
}

async function applyRepresentationPreset(plugin: PluginContext, value: string | undefined, theme?: ViewerTheme) {
    const structures = plugin.managers.structure.hierarchy.selection.structures;
    const representation = String(value ?? 'default').trim().toLowerCase().replaceAll('_', '-');
    const params = theme ? { theme } : undefined;
    switch (representation) {
        case 'default':
            await plugin.managers.structure.component.applyPreset(structures, ViewerAutoPreset, params);
            break;
        case 'auto':
            await plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations.auto, params);
            break;
        case 'cartoon':
            await plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations['polymer-and-ligand'], params);
            break;
        case 'polymer-cartoon':
            await plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations['polymer-cartoon'], params);
            break;
        case 'spacefill':
            await plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations.illustrative, params);
            break;
        case 'surface':
        case 'molecular-surface':
            await plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations['molecular-surface'], params);
            break;
        default:
            throw new Error(`unsupported browser representation: ${representation}`);
    }
}

function openFilesFormat(provider: string | undefined) {
    if (!provider || provider === 'auto') return { name: 'auto', params: {} };
    return { name: 'specific', params: provider };
}

function configureExporter(exporter: ObjExporter | StlExporter, options: BrowserReferenceRequest['exporterOptions']) {
    const quality = options?.primitivesQuality;
    if (quality) exporter.options.primitivesQuality = quality;
}

async function exportWith<D>(
    plugin: PluginContext,
    renderObjects: any[],
    exporter: { add: Function, getData: Function, triangleCount?: number },
    onAdded?: (index: number, triangleCount: number) => void,
) {
    return plugin.runTask(Task.create('Export Mol* browser reference geometry', async ctx => {
        for (let i = 0; i < renderObjects.length; i++) {
            await ctx.update({ message: `Exporting object ${i + 1}/${renderObjects.length}` });
            const before = Number((exporter as any).triangleCount ?? 0);
            await exporter.add(renderObjects[i], plugin.canvas3d!.webgl, ctx);
            const after = Number((exporter as any).triangleCount ?? before);
            onAdded?.(i, after - before);
        }
        return exporter.getData();
    })) as Promise<D>;
}

type RenderObjectStateMetadata = {
    component?: string
    tag?: string
    representation?: string
    colorTheme?: string
    carbonColorTheme?: string
    representationOrder?: number
    visuals?: string[]
    surfaceParams?: GaussianSurfaceParams
}

function renderObjectStateMetadata(plugin: PluginContext) {
    const metadata = new Map<any, RenderObjectStateMetadata>();
    const orderByTag = new Map([
        ['polymer', 0],
        ['ligand', 1],
        ['non-standard', 2],
        ['branched-ball-and-stick', 3],
        ['branched-snfg-3d', 4],
        ['water', 5],
        ['ion', 6],
        ['lipid', 7],
        ['coarse', 8],
    ]);

    for (const cell of plugin.state.data.cells.values()) {
        const repr = (cell.obj as any)?.data?.repr;
        if (!repr?.renderObjects?.length) continue;
        const tags = cell.transform.tags ?? [];
        const tag = tags.find((value: string) => orderByTag.has(value)) ?? tags[0];
        const parent = plugin.state.data.cells.get(cell.transform.parent);
        const parentTags = parent?.transform.tags ?? [];
        const componentTag = parentTags.find((value: string) => value.startsWith('structure-component-static-'));
        const component = componentTag?.slice('structure-component-static-'.length);
        const params = cell.transform.params as any;
        const typeParams = params?.type?.params;
        const stateMetadata: RenderObjectStateMetadata = {
            component,
            tag,
            representation: params?.type?.name,
            colorTheme: params?.colorTheme?.name,
            carbonColorTheme: params?.colorTheme?.params?.carbonColor?.name,
            representationOrder: tag === undefined ? undefined : orderByTag.get(tag),
            visuals: Array.isArray(typeParams?.visuals) ? [...typeParams.visuals] : undefined,
            surfaceParams: params?.type?.name === 'gaussian-surface' ? gaussianSurfaceParamsSummary(typeParams) : undefined,
        };
        for (const renderObject of repr.renderObjects) metadata.set(renderObject, stateMetadata);
    }
    return metadata;
}

function gaussianSurfaceParamsSummary(params: any): GaussianSurfaceParams {
    return {
        quality: params?.quality,
        resolution: params?.resolution,
        smoothness: params?.smoothness,
        radiusOffset: params?.radiusOffset,
        traceOnly: params?.traceOnly,
        tryUseGpu: params?.tryUseGpu,
        visuals: Array.isArray(params?.visuals) ? [...params.visuals] : undefined,
    };
}

function summarizeRenderObjects(renderObjects: any[], stateMetadata: Map<any, RenderObjectStateMetadata>) {
    const objects: RenderObjectSummary[] = [];
    let totalDrawCount = 0;
    let visibleDrawCount = 0;
    let hiddenDrawCount = 0;
    let exportableDrawCount = 0;

    for (let index = 0; index < renderObjects.length; index++) {
        const renderObject = renderObjects[index];
        const values = renderObject.values ?? {};
        const drawCount = Number(value(values.drawCount) ?? 0) || 0;
        const vertexCount = Number(value(values.uVertexCount) ?? value(values.vertexCount) ?? 0) || 0;
        const groupCount = Number(value(values.uGroupCount) ?? value(values.groupCount) ?? 0) || 0;
        const instanceCount = Number(value(values.instanceCount) ?? value(values.uInstanceCount) ?? 0) || 0;
        const geometry = String(value(values.dGeometryType) ?? renderObject.type ?? '<unknown>');
        const visible = renderObject.state?.visible ?? true;
        const sphere = value(values.boundingSphere);
        const primitiveSummary = summarizePrimitiveValues(geometry, values, vertexCount);

        totalDrawCount += drawCount;
        if (visible) {
            visibleDrawCount += drawCount;
            if (drawCount > 0 && instanceCount !== 0) exportableDrawCount += drawCount;
        } else {
            hiddenDrawCount += drawCount;
        }
        objects.push({
            index,
            type: renderObject.type,
            geometry,
            visible,
            drawCount,
            vertexCount,
            groupCount,
            instanceCount,
            ...primitiveSummary,
            ...stateMetadata.get(renderObject),
            boundingSphere: sphere ? sphereSummary(sphere) : undefined,
        });
    }

    return { renderObjects: objects, totalDrawCount, visibleDrawCount, hiddenDrawCount, exportableDrawCount };
}

function summarizePrimitiveValues(geometry: string, values: any, vertexCount: number) {
    if (geometry === 'spheres') {
        return { primitiveCount: vertexCount / 6 };
    }
    if (geometry === 'mesh') {
        const positions = value(values.aPosition) as ArrayLike<number> | undefined;
        const samples = [];
        if (positions) {
            for (let i = 0; i < Math.min(vertexCount, 24); i++) {
                samples.push([
                    positions[i * 3],
                    positions[i * 3 + 1],
                    positions[i * 3 + 2],
                ]);
            }
        }
        return { meshVertexSamples: samples };
    }
    if (geometry !== 'cylinders') return {};

    const caps = value(values.aCap) as ArrayLike<number> | undefined;
    const starts = value(values.aStart) as ArrayLike<number> | undefined;
    const ends = value(values.aEnd) as ArrayLike<number> | undefined;
    const scales = value(values.aScale) as ArrayLike<number> | undefined;
    const groups = value(values.aGroup) as ArrayLike<number> | undefined;
    const histogram: Record<string, number> = {};
    const samples = [];
    if (caps) {
        for (let i = 0; i < vertexCount; i += 6) {
            const cap = String(caps[i] ?? 0);
            histogram[cap] = (histogram[cap] ?? 0) + 1;
            if (samples.length < 8 && starts && ends && scales && groups) {
                samples.push({
                    group: groups[i] ?? 0,
                    cap: caps[i] ?? 0,
                    scale: scales[i] ?? 0,
                    start: [starts[i * 3], starts[i * 3 + 1], starts[i * 3 + 2]],
                    end: [ends[i * 3], ends[i * 3 + 1], ends[i * 3 + 2]],
                });
            }
        }
    }
    return {
        primitiveCount: vertexCount / 6,
        cylinderCapHistogram: histogram,
        cylinderSamples: samples,
    };
}

function sphereSummary(sphere: any) {
    const extrema = Array.isArray(sphere.extrema) ? sphere.extrema.map((e: any) => [e[0], e[1], e[2]]) : undefined;
    return {
        center: [sphere.center?.[0] ?? 0, sphere.center?.[1] ?? 0, sphere.center?.[2] ?? 0],
        radius: Number(sphere.radius ?? 0),
        extremaCount: extrema?.length ?? 0,
        extrema,
    };
}

function value<T>(cell: ValueCell<T> | undefined): T | undefined {
    return cell?.ref?.value;
}

function webglSummary(plugin: PluginContext) {
    const webgl = plugin.canvas3d!.webgl;
    return {
        webgl2: Boolean(webgl.isWebGL2),
        fragDepth: Boolean(webgl.extensions.fragDepth),
        textureFloat: Boolean(webgl.extensions.textureFloat),
    };
}

function debugStlFacet(renderObjects: any[], box: Box3D, sphere: any, stlFacet: number): StlFacetDebug {
    let slotStart = 0;
    const boxMin = toNumberVec3(box.min);
    const boxMax = toNumberVec3(box.max);
    const boxExtrema = boxExtremaDebug(sphere);
    const centerOffset = [
        -((boxMin[0] as number) + (boxMax[0] as number)) * 0.5,
        -((boxMin[1] as number) + (boxMax[1] as number)) * 0.5,
        -((boxMin[2] as number) + (boxMax[2] as number)) * 0.5,
    ];

    for (let renderObjectIndex = 0; renderObjectIndex < renderObjects.length; renderObjectIndex++) {
        const renderObject = renderObjects[renderObjectIndex];
        if (renderObject.state && renderObject.state.visible === false) continue;
        const values = renderObject.values ?? {};
        const drawCount = Number(value(values.drawCount) ?? 0) || 0;
        const instanceCount = Number(value(values.uInstanceCount) ?? value(values.instanceCount) ?? 1) || 1;
        if (drawCount === 0 || instanceCount === 0) continue;
        const objectSlots = drawCount * instanceCount;
        const objectSlotEnd = slotStart + objectSlots;
        if (stlFacet < objectSlotEnd) {
            const localObjectFacet = stlFacet - slotStart;
            const instanceIndex = drawCount > 0 ? Math.floor(localObjectFacet / drawCount) : 0;
            const localFacet = drawCount > 0 ? localObjectFacet % drawCount : 0;
            const base = baseFacetDebug(renderObject, renderObjectIndex, stlFacet, slotStart, objectSlotEnd, localFacet, instanceIndex, centerOffset, boxMin, boxMax, boxExtrema);
            if (String(value(values.dGeometryType) ?? renderObject.type) !== 'mesh') {
                return { ...base, found: false, error: 'debugStlFacet currently supports mesh render objects only' };
            }
            if (localFacet % 3 !== 0) {
                return { ...base, found: false, sparseSlot: localFacet % 3, error: 'STL sparse slot is empty because Mol* stores one triangle every three draw slots' };
            }
            return debugMeshStlFacet(renderObject, base, localFacet, instanceIndex, centerOffset);
        }
        slotStart = objectSlotEnd;
    }

    return {
        stlFacet,
        found: false,
        centerOffset,
        boxMin,
        boxMax,
        ...boxExtrema,
        error: `facet outside exported render-object slot range (${slotStart})`,
    };
}

function baseFacetDebug(
    renderObject: any,
    renderObjectIndex: number,
    stlFacet: number,
    objectSlotStart: number,
    objectSlotEnd: number,
    localFacet: number,
    instanceIndex: number,
    centerOffset: number[],
    boxMin: number[],
    boxMax: number[],
    boxExtrema: ReturnType<typeof boxExtremaDebug>,
): StlFacetDebug {
    const values = renderObject.values ?? {};
    return {
        stlFacet,
        found: true,
        renderObjectIndex,
        renderObjectType: String(renderObject.type),
        geometry: String(value(values.dGeometryType) ?? renderObject.type ?? '<unknown>'),
        visible: renderObject.state?.visible ?? true,
        localFacet,
        localTriangle: Math.floor(localFacet / 3),
        instanceIndex,
        objectSlotStart,
        objectSlotEnd,
        sparseSlot: localFacet % 3,
        centerOffset,
        boxMin,
        boxMax,
        ...boxExtrema,
    };
}

function debugMeshStlFacet(
    renderObject: any,
    base: StlFacetDebug,
    localFacet: number,
    instanceIndex: number,
    centerOffset: number[],
): StlFacetDebug {
    const values = renderObject.values ?? {};
    const positions = value(values.aPosition) as Float32Array | undefined;
    const elements = value(values.elements) as Uint32Array | undefined;
    const transforms = value(values.aTransform) as Float32Array | undefined;
    if (!positions || !elements) return { ...base, found: false, error: 'mesh render object is missing aPosition/elements buffers' };
    const vertexIndices = [elements[localFacet], elements[localFacet + 1], elements[localFacet + 2]];
    if (vertexIndices.some(index => index === undefined)) {
        return { ...base, found: false, vertexIndices, error: 'local facet is outside mesh element buffer' };
    }

    const transform = centerTransformForInstance(centerOffset, transforms, instanceIndex);
    const rawVertices = vertexIndices.map(index => readVec3(positions, index * 3));
    const centeredVertices = rawVertices.map(vertex => toF32Vec3(transformMat4(vertex, transform)));
    const triangleNormal = toF32Vec3(triangleNormalFromVertices(centeredVertices[0], centeredVertices[1], centeredVertices[2]));
    return {
        ...base,
        vertexIndices,
        rawVertices,
        centeredVertices,
        centeredVertexBits: centeredVertices.map(vec3Bits),
        triangleNormal,
        triangleNormalBits: vec3Bits(triangleNormal),
    };
}

function centerTransformForInstance(centerOffset: number[], transforms: Float32Array | undefined, instanceIndex: number) {
    const center = [
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        centerOffset[0], centerOffset[1], centerOffset[2], 1,
    ];
    const instance = transforms && transforms.length >= (instanceIndex + 1) * 16
        ? Array.from(transforms.slice(instanceIndex * 16, instanceIndex * 16 + 16))
        : [
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1,
        ];
    return mat4Mul(center, instance);
}

function mat4Mul(a: number[], b: number[]) {
    const out = new Array<number>(16);
    const a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3];
    const a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7];
    const a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11];
    const a30 = a[12], a31 = a[13], a32 = a[14], a33 = a[15];
    let b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
    out[0] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
    out[1] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
    out[2] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
    out[3] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;
    b0 = b[4]; b1 = b[5]; b2 = b[6]; b3 = b[7];
    out[4] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
    out[5] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
    out[6] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
    out[7] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;
    b0 = b[8]; b1 = b[9]; b2 = b[10]; b3 = b[11];
    out[8] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
    out[9] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
    out[10] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
    out[11] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;
    b0 = b[12]; b1 = b[13]; b2 = b[14]; b3 = b[15];
    out[12] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
    out[13] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
    out[14] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
    out[15] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;
    return out;
}

function transformMat4(vertex: number[], m: number[]) {
    const x = vertex[0], y = vertex[1], z = vertex[2];
    const w = 1 / ((m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0);
    return [
        (m[0] * x + m[4] * y + m[8] * z + m[12]) * w,
        (m[1] * x + m[5] * y + m[9] * z + m[13]) * w,
        (m[2] * x + m[6] * y + m[10] * z + m[14]) * w,
    ];
}

function triangleNormalFromVertices(a: number[], b: number[], c: number[]) {
    const ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    const ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    const n = [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ];
    const lenSq = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
    if (lenSq <= 0) return [0, 0, 0];
    const scale = 1 / Math.sqrt(lenSq);
    return [n[0] * scale, n[1] * scale, n[2] * scale];
}

function readVec3(array: Float32Array, offset: number) {
    return [array[offset], array[offset + 1], array[offset + 2]];
}

function toNumberVec3(values: ArrayLike<number>) {
    return [Number(values[0]), Number(values[1]), Number(values[2])];
}

function boxExtremaDebug(sphere: any) {
    const extrema = Array.isArray(sphere?.extrema) ? sphere.extrema : [];
    const min = [Infinity, Infinity, Infinity];
    const max = [-Infinity, -Infinity, -Infinity];
    const minIndices: (number | null)[] = [null, null, null];
    const maxIndices: (number | null)[] = [null, null, null];
    const minPoints: (number[] | null)[] = [null, null, null];
    const maxPoints: (number[] | null)[] = [null, null, null];
    for (let index = 0; index < extrema.length; index++) {
        const point = toNumberVec3(extrema[index]);
        for (let axis = 0; axis < 3; axis++) {
            if (point[axis] < min[axis]) {
                min[axis] = point[axis];
                minIndices[axis] = index;
                minPoints[axis] = point;
            }
            if (point[axis] > max[axis]) {
                max[axis] = point[axis];
                maxIndices[axis] = index;
                maxPoints[axis] = point;
            }
        }
    }
    return {
        boxExtremaCount: extrema.length,
        boxMinIndices: minIndices,
        boxMaxIndices: maxIndices,
        boxMinPoints: minPoints,
        boxMaxPoints: maxPoints,
    };
}

function toF32Vec3(values: number[]) {
    const f32 = new Float32Array(3);
    f32[0] = values[0];
    f32[1] = values[1];
    f32[2] = values[2];
    return [f32[0], f32[1], f32[2]];
}

function vec3Bits(values: number[]) {
    return values.map(value => float32Hex(value));
}

function float32Hex(value: number) {
    const bytes = new ArrayBuffer(4);
    new DataView(bytes).setFloat32(0, value, true);
    const bits = new DataView(bytes).getUint32(0, true);
    return `0x${bits.toString(16).padStart(8, '0')}`;
}

async function uploadText(id: string, format: string, text: string) {
    await uploadBytes(id, format, new TextEncoder().encode(text));
}

async function uploadBytes(id: string, format: string, bytes: Uint8Array) {
    const response = await fetch(`/upload/${encodeURIComponent(id)}/${encodeURIComponent(format)}`, {
        method: 'POST',
        headers: { 'content-type': 'application/octet-stream' },
        body: bytes,
    });
    if (!response.ok) {
        throw new Error(`upload ${format} failed: ${response.status} ${await response.text()}`);
    }
}

function animationFrames(count: number) {
    return new Promise<void>(resolve => {
        const tick = () => {
            count -= 1;
            if (count <= 0) resolve();
            else requestAnimationFrame(tick);
        };
        requestAnimationFrame(tick);
    });
}
