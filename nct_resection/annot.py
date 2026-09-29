"""Minimal FreeSurfer .annot reader, and parcel adjacency built from it.

nibabel would do the reading in one call, but the format is simple enough to
parse directly and it keeps this project dependency-free.

Binary layout, big-endian:
    int32                    vertex count n
    n * (int32, int32)       (vertex index, packed RGB annotation value)
    int32                    colortable present flag
    int32                    version marker (-2 for the current format)
    int32                    entry count
    int32 + chars            original colortable filename
    int32                    entries to read
    per entry:
        int32                structure index
        int32 + chars        structure name
        int32 R, G, B, A

A vertex's parcel is found by matching its packed value against
R + G*2^8 + B*2^16 for each colortable entry.
"""

import base64
import gzip
import struct
import xml.etree.ElementTree as ET
import zlib

import numpy as np

_GIFTI_DTYPES = {
    "NIFTI_TYPE_UINT8": np.uint8,
    "NIFTI_TYPE_INT32": np.int32,
    "NIFTI_TYPE_FLOAT32": np.float32,
    "NIFTI_TYPE_FLOAT64": np.float64,
}


def read_gifti_surface(path):
    """Return (coords, faces) from a GIFTI .surf.gii file.

    Handles the encodings the format actually uses in the wild: raw ASCII,
    plain Base64, and GZip+Base64. nibabel would be the normal way to do this;
    parsing directly keeps the project dependency-free.
    """
    root = ET.parse(path).getroot()
    coords = faces = None

    for da in root.iter("DataArray"):
        intent = da.get("Intent")
        dtype = _GIFTI_DTYPES.get(da.get("DataType"))
        encoding = da.get("Encoding")
        endian = da.get("Endian", "LittleEndian")
        n_dim = int(da.get("Dimensionality", "2"))
        shape = tuple(int(da.get(f"Dim{i}")) for i in range(n_dim))

        raw = da.find("Data").text.strip()
        if encoding == "ASCII":
            arr = np.fromstring(raw, sep=" ", dtype=dtype)
        else:
            blob = base64.b64decode(raw)
            if encoding == "GZipBase64Binary":
                try:
                    blob = gzip.decompress(blob)
                except OSError:
                    blob = zlib.decompress(blob)
            arr = np.frombuffer(blob, dtype=dtype)

        if endian == "BigEndian":
            arr = arr.byteswap()
        arr = arr.reshape(shape)

        if intent == "NIFTI_INTENT_POINTSET":
            coords = np.asarray(arr, dtype=float)
        elif intent == "NIFTI_INTENT_TRIANGLE":
            faces = np.asarray(arr, dtype=np.int64)

    if coords is None or faces is None:
        raise ValueError(f"{path}: missing pointset or triangle array")
    return coords, faces


def read_annot(path):
    """Return (labels, names) where labels[i] is the parcel index of vertex i.

    Vertices with no assigned parcel get -1.
    """
    with open(path, "rb") as f:
        buf = f.read()

    pos = 0

    def i32():
        nonlocal pos
        (val,) = struct.unpack_from(">i", buf, pos)
        pos += 4
        return val

    n_vertices = i32()
    raw = np.frombuffer(buf, dtype=">i4", count=n_vertices * 2, offset=pos)
    pos += n_vertices * 2 * 4
    vertex_ids = raw[0::2].astype(np.int64)
    packed = raw[1::2].astype(np.int64)

    if not i32():
        raise ValueError("annot file has no embedded colortable")

    version = i32()
    if version != -2:
        raise ValueError(f"unsupported colortable version {version}")

    n_entries = i32()
    orig_len = i32()
    pos += orig_len  # original colortable path, unused
    n_read = i32()

    names = [None] * n_entries
    codes = np.full(n_entries, -1, dtype=np.int64)
    for _ in range(n_read):
        idx = i32()
        name_len = i32()
        name = buf[pos:pos + name_len - 1].decode("utf-8", "replace")
        pos += name_len
        r, g, b, a = i32(), i32(), i32(), i32()
        names[idx] = name
        codes[idx] = r + (g << 8) + (b << 16) + (a << 24)

    # Map each vertex's packed value to its colortable entry.
    lookup = {int(c): i for i, c in enumerate(codes) if c >= 0}
    labels = np.full(n_vertices, -1, dtype=np.int64)
    for value, entry in lookup.items():
        labels[packed == value] = entry

    ordered = np.full(n_vertices, -1, dtype=np.int64)
    ordered[vertex_ids] = labels
    return ordered, names


def adjacency_from_faces(labels, faces, n_parcels):
    """Two parcels are adjacent if any mesh triangle touches both.

    This is the real definition of anatomical contiguity: it uses the surface
    topology rather than guessing from connectivity strength.
    """
    adj = np.zeros((n_parcels, n_parcels), dtype=bool)
    tri = labels[faces]
    for a, b in ((0, 1), (1, 2), (0, 2)):
        pairs = tri[:, [a, b]]
        valid = (pairs[:, 0] >= 0) & (pairs[:, 1] >= 0)
        pairs = pairs[valid]
        differing = pairs[pairs[:, 0] != pairs[:, 1]]
        adj[differing[:, 0], differing[:, 1]] = True
        adj[differing[:, 1], differing[:, 0]] = True
    return adj


def build_hcp_mmp1_adjacency(data_dir="data", cache="data/hcp_mmp1_adjacency.npy"):
    """Real 360x360 parcel adjacency in the connectome's row order.

    Verified rather than assumed: `originalParcelIDs` in the Rosen & Halgren
    ordering file matches the left annot colortable entries 1..180 followed by
    the right entries 1..180, exactly, all 360 in order.

    Adjacency is block diagonal by hemisphere. The two surfaces are separate
    meshes, and cortical tissue is not contiguous across the midline, so a
    resection cannot grow from one hemisphere into the other along the surface.
    """
    import os

    if cache and os.path.exists(cache):
        return np.load(cache)

    blocks = []
    for hemi, surf in (("lh", "L"), ("rh", "R")):
        labels, _ = read_annot(f"{data_dir}/{hemi}.HCP-MMP1.annot")
        _, faces = read_gifti_surface(f"{data_dir}/fsaverage_{surf}_white.surf.gii")
        full = adjacency_from_faces(labels, faces, 181)
        blocks.append(full[1:181, 1:181])  # drop the '???' entry

    n = 360
    adj = np.zeros((n, n), dtype=bool)
    adj[:180, :180] = blocks[0]
    adj[180:, 180:] = blocks[1]

    if cache:
        np.save(cache, adj)
    return adj


def adjacency_from_coords(labels, coords, n_parcels, tolerance=1.5):
    """Fallback when face topology is unavailable.

    Two parcels are adjacent if any vertex of one lies within `tolerance` of a
    vertex of the other. Slower and less exact than using faces, but it only
    needs vertex coordinates.
    """
    from scipy.spatial import cKDTree

    adj = np.zeros((n_parcels, n_parcels), dtype=bool)
    trees = {}
    for p in range(n_parcels):
        pts = coords[labels == p]
        if pts.size:
            trees[p] = cKDTree(pts)

    for p, tree_p in trees.items():
        for q, tree_q in trees.items():
            if q <= p:
                continue
            if tree_p.query_ball_tree(tree_q, r=tolerance, p=2):
                if any(tree_p.query_ball_tree(tree_q, r=tolerance)):
                    adj[p, q] = adj[q, p] = True
    return adj
