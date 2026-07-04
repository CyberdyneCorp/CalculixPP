"""Shared toolkit for the mesh-based FE examples.

Turns a surface mesh (STL / OBJ) into a solid C3D8 voxel mesh, builds a linear-static
CalculiX++ deck, and renders / exports the result. Two examples use it: the table
(table_stress/) and the antenna (Antenna_1/).

Dependency-light: numpy + scipy, plus matplotlib for the renders.

    IMPORTANT: any script using this must `import calculixpp` BEFORE importing this
    module (which pulls in numpy). The compiled module links the system libstdc++;
    importing numpy first can load an older libstdc++ and shadow the symbols it needs.
"""

import numpy as np
from scipy import ndimage

# --------------------------------------------------------------------------- IO

def read_stl_triangles(path):
    """Parse an ASCII STL into an (F, 3, 3) array of triangle vertices."""
    verts = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if parts and parts[0] == "vertex":
                verts.append([float(v) for v in parts[1:4]])
    return np.array(verts).reshape(-1, 3, 3)


def read_obj_triangles(path):
    """Parse a Wavefront OBJ into an (F, 3, 3) triangle array (polygons fan-triangulated)."""
    verts, faces = [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "v":
                verts.append([float(x) for x in parts[1:4]])
            elif parts[0] == "f":
                idx = [int(t.split("/")[0]) - 1 for t in parts[1:]]
                for k in range(1, len(idx) - 1):        # fan-triangulate quads/n-gons
                    faces.append([idx[0], idx[k], idx[k + 1]])
    return np.array(verts)[np.array(faces)]


# ------------------------------------------------------------------- voxelising

def voxelize_winding(tri, voxel, chunk=4000):
    """Occupancy grid via the generalized winding number (robust inside test).

    A voxel centre is inside when |w(p)| > 0.5. The winding number sums each triangle's
    signed solid angle at p, so it is robust to *overlapping closed shells* (common in
    game-art assets) and small non-watertightness — it evaluates the true union. Use this
    for genuinely solid bodies. Returns (occ, origin).
    """
    lo = tri.reshape(-1, 3).min(0)
    hi = tri.reshape(-1, 3).max(0)
    nx, ny, nz = np.maximum(1, np.ceil((hi - lo) / voxel).astype(int))
    gx = lo[0] + (np.arange(nx) + 0.5) * voxel
    gy = lo[1] + (np.arange(ny) + 0.5) * voxel
    gz = lo[2] + (np.arange(nz) + 0.5) * voxel
    pts = np.stack(np.meshgrid(gx, gy, gz, indexing="ij"), -1).reshape(-1, 3)
    A, B, C = tri[:, 0], tri[:, 1], tri[:, 2]
    w = np.empty(len(pts))
    for s in range(0, len(pts), chunk):
        p = pts[s:s + chunk][:, None, :]
        a, b, c = A - p, B - p, C - p
        la = np.linalg.norm(a, axis=2)
        lb = np.linalg.norm(b, axis=2)
        lc = np.linalg.norm(c, axis=2)
        num = (a * np.cross(b, c)).sum(2)
        den = (la * lb * lc + (a * b).sum(2) * lc
               + (b * c).sum(2) * la + (c * a).sum(2) * lb)
        w[s:s + chunk] = 2.0 * np.arctan2(num, den).sum(1)
    occ = (np.abs(w) > 0.5 * (4 * np.pi)).reshape(nx, ny, nz)
    return occ, lo


def voxelize_solid_fill(tri, voxel, axis=2):
    """Occupancy of the SOLID bounded by a surface, via surface raster + per-slice fill.

    For a thin-walled / hollow surface (e.g. a tubular mast) an inside test finds almost
    nothing — the wall is thinner than a voxel. Instead rasterise the surface into the
    grid (a shell of occupied cells) and, on each slice perpendicular to `axis`, fill the
    interior enclosed by that shell (`binary_fill_holes`). The result is the solid body
    the outer surface bounds. Returns (occ, origin).
    """
    lo = tri.reshape(-1, 3).min(0)
    hi = tri.reshape(-1, 3).max(0)
    dims = np.maximum(1, np.ceil((hi - lo) / voxel).astype(int))
    occ = np.zeros(tuple(dims), bool)
    for a, b, c in tri:                                 # rasterise the surface
        step = max(np.linalg.norm(b - a), np.linalg.norm(c - a), np.linalg.norm(c - b))
        n = max(1, int(step / (voxel * 0.5)))
        u, v = np.meshgrid(np.linspace(0, 1, n + 1), np.linspace(0, 1, n + 1))
        m = u + v <= 1.0
        pts = a + u[m][:, None] * (b - a) + v[m][:, None] * (c - a)
        ijk = np.clip(np.floor((pts - lo) / voxel).astype(int), 0, dims - 1)
        occ[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = True
    for s in range(dims[axis]):                         # fill each cross-section
        sl = [slice(None)] * 3
        sl[axis] = s
        occ[tuple(sl)] = ndimage.binary_fill_holes(occ[tuple(sl)])
    return occ, lo


def largest_component(occ):
    """Keep the largest 6-connected (face-adjacent) component — a single solid, no specks."""
    labels, count = ndimage.label(occ)
    if count <= 1:
        return occ
    sizes = ndimage.sum(np.ones_like(labels), labels, range(1, count + 1))
    return labels == (int(np.argmax(sizes)) + 1)


def hollow_shell(occ, wall_voxels=1):
    """Turn a solid occupancy into a hollow shell: the solid minus an eroded core.

    Leaves a wall `wall_voxels` cells thick (the thinnest robust wall for solid elements
    is 1 voxel). Note the voxel wall is far thicker than a real thin steel wall, so match
    the *mass* via a calibrated density rather than trusting the modelled volume.
    """
    core = ndimage.binary_erosion(occ, iterations=wall_voxels)
    return occ & ~core


def surface_area(tri):
    """Total area of a triangle mesh (F, 3, 3) — e.g. to size a real thin-wall mass."""
    cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    return 0.5 * np.linalg.norm(cross, axis=1).sum()


# ------------------------------------------------------------------ FE hex mesh

class VoxelMesh:
    """A C3D8 hex mesh built from an occupancy grid, with shared (coord-keyed) nodes."""

    def __init__(self, occ, origin, voxel):
        self.voxel = voxel
        self._origin = np.asarray(origin, float)
        self._nodes = {}
        self.elements = []
        for i, j, k in np.argwhere(occ):
            self.elements.append([
                self._node(i, j, k),         self._node(i + 1, j, k),
                self._node(i + 1, j + 1, k), self._node(i, j + 1, k),
                self._node(i, j, k + 1),     self._node(i + 1, j, k + 1),
                self._node(i + 1, j + 1, k + 1), self._node(i, j + 1, k + 1),
            ])
        self.coords = self._coords()

    def _node(self, i, j, k):
        key = (int(i), int(j), int(k))
        if key not in self._nodes:
            self._nodes[key] = len(self._nodes) + 1
        return self._nodes[key]

    def _coords(self):
        out = np.zeros((len(self._nodes), 3))
        for (i, j, k), nid in self._nodes.items():
            out[nid - 1] = self._origin + np.array([i, j, k]) * self.voxel
        return out

    def nodes_on_plane(self, axis, value, tol=1e-6):
        """1-based ids of nodes with coords[axis] ≈ value."""
        return [i + 1 for i in range(len(self.coords))
                if abs(self.coords[i, axis] - value) < tol]

    def nodes_where(self, mask):
        """1-based ids where mask(coords) is True (mask: (N,3) array -> (N,) bool)."""
        return list(np.where(mask(self.coords))[0] + 1)


def static_deck(coords, elements, fixed_ids, cloads, E=210000.0, nu=0.3,
                density=7.85e-9, material="STEEL", gravity=None):
    """Assemble a linear-static .inp deck.

    cloads: iterable of (node_id, dof, value). gravity: optional (g, (nx, ny, nz)) body
    load — emits `*DLOAD EALL,GRAV,g,nx,ny,nz` so the material density gives self-weight.
    """
    lines = ["*NODE"]
    lines += [f"{i + 1},{coords[i, 0]},{coords[i, 1]},{coords[i, 2]}"
              for i in range(len(coords))]
    lines.append("*ELEMENT,TYPE=C3D8,ELSET=EALL")
    lines += [f"{e},{','.join(map(str, nodes))}" for e, nodes in enumerate(elements, 1)]
    lines.append("*BOUNDARY")
    lines += [f"{n},1,3" for n in fixed_ids]
    lines.append("*CLOAD")
    lines += [f"{n},{dof},{val}" for n, dof, val in cloads]
    lines += [
        f"*MATERIAL,NAME={material}",
        "*ELASTIC", f"{E},{nu}",
        "*DENSITY", f"{density}",
        f"*SOLID SECTION,ELSET=EALL,MATERIAL={material}",
        "*STEP", "*STATIC",
    ]
    if gravity is not None:
        g, (nx, ny, nz) = gravity
        lines += ["*DLOAD", f"EALL,GRAV,{g},{nx},{ny},{nz}"]
    lines.append("*END STEP")
    return "\n".join(lines)


# ------------------------------------------------------------- post-processing

_HEX_FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
              (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


def von_mises(stress):
    """Nodal von Mises stress from an (N, 6) tensor (xx, yy, zz, xy, xz, yz)."""
    sxx, syy, szz, sxy, sxz, syz = stress.T
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
                   + 3.0 * (sxy ** 2 + sxz ** 2 + syz ** 2))


def align_to_ids(result, n_nodes):
    """Row order mapping a solve() result (solver node order) onto 1-based ids 1..N."""
    id_to_row = {int(nid): i for i, nid in enumerate(result["node_ids"])}
    return np.array([id_to_row[i + 1] for i in range(n_nodes)])


def boundary_faces(elements):
    """Outer (once-referenced) quad faces of a hex mesh, as node-id quads."""
    seen, owner = {}, {}
    for elem in elements:
        for face in _HEX_FACES:
            quad = tuple(elem[i] for i in face)
            key = tuple(sorted(quad))
            seen[key] = seen.get(key, 0) + 1
            owner[key] = quad
    return [owner[k] for k, count in seen.items() if count == 1]


def _new_3d_axes(fig, pos, coords):
    ax = fig.add_subplot(pos, projection="3d")
    ax.set_xlim(coords[:, 0].min(), coords[:, 0].max())
    ax.set_ylim(coords[:, 1].min(), coords[:, 1].max())
    ax.set_zlim(coords[:, 2].min(), coords[:, 2].max())
    ax.set_box_aspect(np.ptp(coords, axis=0))
    return ax


def preview_geometry(coords, elements, path, title,
                     views=((20, -60), (15, 30), (20, 120), (-10, 30))):
    """Render the undeformed hex mesh in flat grey from several angles (a quality check)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    polys = [[coords[n - 1] for n in quad] for quad in boundary_faces(elements)]
    fig = plt.figure(figsize=(13, 9))
    for k, (elev, azim) in enumerate(views):
        ax = _new_3d_axes(fig, 221 + k, coords)
        ax.add_collection3d(Poly3DCollection(polys, facecolor="#b0b8c8",
                            edgecolor=(0, 0, 0, 0.25), linewidths=0.2))
        ax.view_init(elev=elev, azim=azim)
        ax.set_title(f"elev={elev}° azim={azim}°", fontsize=9)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def render_stress(coords, elements, vm, disp, path, title, warp_fraction=0.1,
                  max_scale=None, elev=22, azim=-60, figsize=(11, 7)):
    """Render the deformed surface coloured by von Mises stress -> PNG at `path`.

    Deformation is exaggerated by ``warp_fraction × model_size / max|disp|`` so bending
    is visible, but capped at ``max_scale`` — without a cap, a near-rigid response (a
    stiff part barely deflecting) blows up to a misleading multi-million× "bent" shape.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    reach = np.ptp(coords, axis=0).max()
    scale = warp_fraction * reach / max(np.abs(disp).max(), 1e-12)
    if max_scale is not None:
        scale = min(scale, max_scale)
    warped = coords + scale * disp

    faces = boundary_faces(elements)
    polys = [[warped[n - 1] for n in quad] for quad in faces]
    face_vm = np.array([np.mean([vm[n - 1] for n in quad]) for quad in faces])
    norm = plt.Normalize(face_vm.min(), face_vm.max())
    cmap = matplotlib.colormaps["viridis"]

    fig = plt.figure(figsize=figsize)
    ax = _new_3d_axes(fig, 111, warped)
    ax.add_collection3d(Poly3DCollection(polys, facecolors=cmap(norm(face_vm)),
                        edgecolors=(0, 0, 0, 0.15), linewidths=0.3))
    ax.set(xlabel="x [mm]", ylabel="y [mm]", zlabel="z [mm]")
    ax.set_title(f"{title}\n(deformation ×{scale:.0f})")
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                 shrink=0.6, label="von Mises [MPa]")
    ax.view_init(elev=elev, azim=azim)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def write_vtk(path, coords, elements, vm, disp):
    """Legacy VTK UNSTRUCTURED_GRID (VTK_HEXAHEDRON) for ParaView."""
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\nCalculiX++ stress\nASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(coords)} double\n")
        for x, y, z in coords:
            f.write(f"{x} {y} {z}\n")
        f.write(f"\nCELLS {len(elements)} {len(elements) * 9}\n")
        for nodes in elements:
            f.write("8 " + " ".join(str(n - 1) for n in nodes) + "\n")
        f.write(f"\nCELL_TYPES {len(elements)}\n")
        f.write("\n".join(["12"] * len(elements)) + "\n")
        f.write(f"\nPOINT_DATA {len(coords)}\n")
        f.write("SCALARS von_mises double 1\nLOOKUP_TABLE default\n")
        f.write("\n".join(f"{v}" for v in vm) + "\n")
        f.write("VECTORS displacement double\n")
        for ux, uy, uz in disp:
            f.write(f"{ux} {uy} {uz}\n")
