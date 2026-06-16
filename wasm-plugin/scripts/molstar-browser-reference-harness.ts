import '../artifacts/molstar/src/mol-util/polyfill';
import { ObjExporter } from '../artifacts/molstar/src/extensions/geo-export/obj-exporter';
import { StlExporter } from '../artifacts/molstar/src/extensions/geo-export/stl-exporter';
import { Box3D } from '../artifacts/molstar/src/mol-math/geometry';
import { PluginContext } from '../artifacts/molstar/src/mol-plugin/context';
import { DefaultPluginSpec } from '../artifacts/molstar/src/mol-plugin/spec';
import { OpenFiles } from '../artifacts/molstar/src/mol-plugin-state/actions/file';
import { Asset } from '../artifacts/molstar/src/mol-util/assets';
import { Task } from '../artifacts/molstar/src/mol-task';
import type { ValueCell } from '../artifacts/molstar/src/mol-util/value-cell';

type ExportFormat = 'obj' | 'stl'

type BrowserReferenceRequest = {
    id: string
    fixtureUrl: string
    fileName: string
    formats: ExportFormat[]
    objExportBasename: string
    dataFormat?: string
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

    plugin = new PluginContext(DefaultPluginSpec());
    await plugin.init();
    const mounted = await plugin.mountAsync(target, { checkeredCanvasBackground: false });
    if (!mounted) throw new Error('failed to mount Mol* plugin');
    await plugin.canvas3dInitialized;
    return plugin;
}

window.molfigBrowserReferenceExport = async (request: BrowserReferenceRequest) => {
    const plugin = await getPlugin();
    const bytes = new Uint8Array(await (await fetch(request.fixtureUrl)).arrayBuffer());
    const file = new File([bytes], request.fileName);
    await plugin.runTask(plugin.state.data.applyAction(OpenFiles, {
        files: [Asset.File(file)],
        format: openFilesFormat(request.dataFormat),
        visuals: true,
    }));
    plugin.canvas3d?.commit(true);
    await animationFrames(3);

    const renderObjects = plugin.canvas3d?.getRenderObjects() ?? [];
    if (renderObjects.length === 0) throw new Error('Mol* produced no render objects');
    const summary = summarizeRenderObjects(renderObjects);
    const sphere = plugin.canvas3d!.boundingSphereVisible;
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
            const data = await exportWith(plugin, renderObjects, exporter);
            await uploadBytes(request.id, 'stl', data.stl);
            uploads.stl = { byteLength: data.stl.byteLength };
        }
    }

    return {
        id: request.id,
        ...summary,
        uploads,
        webgl: webglSummary(plugin),
        debug,
    };
};

function openFilesFormat(provider: string | undefined) {
    if (!provider || provider === 'auto') return { name: 'auto', params: {} };
    return { name: 'specific', params: provider };
}

function configureExporter(exporter: ObjExporter | StlExporter, options: BrowserReferenceRequest['exporterOptions']) {
    const quality = options?.primitivesQuality;
    if (quality) exporter.options.primitivesQuality = quality;
}

async function exportWith<D>(plugin: PluginContext, renderObjects: any[], exporter: { add: Function, getData: Function }) {
    return plugin.runTask(Task.create('Export Mol* browser reference geometry', async ctx => {
        for (let i = 0; i < renderObjects.length; i++) {
            await ctx.update({ message: `Exporting object ${i + 1}/${renderObjects.length}` });
            await exporter.add(renderObjects[i], plugin.canvas3d!.webgl, ctx);
        }
        return exporter.getData();
    })) as Promise<D>;
}

function summarizeRenderObjects(renderObjects: any[]) {
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
            boundingSphere: sphere ? sphereSummary(sphere) : undefined,
        });
    }

    return { renderObjects: objects, totalDrawCount, visibleDrawCount, hiddenDrawCount, exportableDrawCount };
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
