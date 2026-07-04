#!/usr/bin/env python3
"""Stress analysis of a real FBX table model under a 100 kg central object.

Unlike `table_stress.py` (which builds a clean parametric mesh), this reads the actual
`table.fbx` asset and turns it into a finite-element model, then runs the same
linear-static analysis with CalculiX++.

An FBX is a *surface* mesh (a triangle skin), while FEA needs a *solid* (volume) mesh.
The bridge here is **voxelisation** (`femkit.voxelize_winding`): sample a 3-D grid over
the model, keep the cells inside the watertight table shell, and emit one C3D8
hexahedron per inside cell. The winding-number inside test is robust to the overlapping
closed shells (tabletop + apron + legs) this asset is built from. The watertight surface
is extracted from the FBX by Blender (see `extract_surface.py`).

Pipeline:
    table.fbx --(Blender)--> table_solid.stl --(voxelise)--> C3D8 mesh --(CalculiX++)--> stress

Run:
    # 1) one-time: extract the watertight table body from the FBX (needs Blender)
    blender --background --python extract_surface.py
    # 2) analyse + visualise (needs numpy, scipy, matplotlib, calculixpp)
    PYTHONPATH=/path/to/CalculixPP/build/python python3 fbx_table_stress.py [--preview]
"""

# Import calculixpp before numpy (see femkit for the libstdc++ rationale).
import calculixpp as cx

import os
import sys
import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import femkit  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
STL_PATH = os.path.join(HERE, "table_solid.stl")

# The FBX is modelled in metres; scale to the CalculiX consistent set (N, mm, MPa, t).
SCALE_TO_MM = 1000.0
VOXEL_MM = 15.0                 # voxel edge length -> mesh resolution
OBJECT_MASS = 0.1               # tonne (= 100 kg)
GRAVITY = 9810.0                # mm/s^2
OBJECT_HALF_FOOTPRINT = 120.0   # object rests on a ~240 x 240 mm central patch


def main():
    if not os.path.exists(STL_PATH):
        raise SystemExit(
            f"{STL_PATH} not found.\nExtract it from the FBX first:\n"
            f"    blender --background --python {os.path.join(HERE, 'extract_surface.py')}")

    tri = femkit.read_stl_triangles(STL_PATH) * SCALE_TO_MM
    print(f"FBX surface: {len(tri)} triangles, "
          f"bounding box {np.ptp(tri.reshape(-1, 3), axis=0).round(0)} mm")

    print("Voxelising (winding number)...", flush=True)
    occ, origin = femkit.voxelize_winding(tri, VOXEL_MM)
    solid = femkit.largest_component(occ)
    mesh = femkit.VoxelMesh(solid, origin, VOXEL_MM)
    print(f"Voxel mesh: {occ.sum()} raw -> {solid.sum()} solid hexes "
          f"({VOXEL_MM:.0f} mm voxels), {len(mesh.coords)} nodes")

    # Sanity-check before solving: a single connected solid standing on its legs.
    print(f"  legs on the floor    : {ndimage.label(solid[:, :, 0])[1]}")
    if "--preview" in sys.argv:
        femkit.preview_geometry(mesh.coords, mesh.elements,
                                os.path.join(HERE, "fbx_table_mesh.png"),
                                f"FBX table voxel mesh (undeformed) — {len(mesh.elements)} hexes")
        print("Wrote fbx_table_mesh.png")
        return

    coords = mesh.coords
    z_max = coords[:, 2].max()
    cx_, cy_ = coords[:, 0].mean(), coords[:, 1].mean()
    fixed = mesh.nodes_on_plane(axis=2, value=coords[:, 2].min())  # clamp leg bottoms
    load = mesh.nodes_where(lambda c: (np.abs(c[:, 2] - z_max) < 1e-6)
                            & (np.abs(c[:, 0] - cx_) < OBJECT_HALF_FOOTPRINT)
                            & (np.abs(c[:, 1] - cy_) < OBJECT_HALF_FOOTPRINT))
    force = -OBJECT_MASS * GRAVITY / len(load)
    print(f"  clamped bottom nodes : {len(fixed)}")
    print(f"  loaded top nodes     : {len(load)}  (100 kg = {OBJECT_MASS * GRAVITY:.0f} N)")

    deck = femkit.static_deck(coords, mesh.elements, fixed,
                              [(n, 3, force) for n in load])
    result = cx.solve_text(deck)

    order = femkit.align_to_ids(result, len(coords))
    out_coords = result["node_coords"][order]
    disp = result["displacement"][order]
    vm = femkit.von_mises(result["stress"][order])

    print("\nResults")
    print(f"  max vertical deflection : {disp[:, 2].min():+.4f} mm")
    print(f"  max von Mises stress    : {vm.max():.2f} MPa")
    print(f"  vertical reaction sum   : {result['reaction'][:, 2].sum():+.1f} N "
          f"(balances the {OBJECT_MASS * GRAVITY:.0f} N weight)")

    png = os.path.join(HERE, "fbx_table_von_mises.png")
    vtk = os.path.join(HERE, "fbx_table_stress.vtk")
    femkit.render_stress(out_coords, mesh.elements, vm, disp, png,
                         "FBX table — von Mises under 100 kg at centre")
    femkit.write_vtk(vtk, out_coords, mesh.elements, vm, disp)
    print(f"\nWrote {png}\n      {vtk}")


if __name__ == "__main__":
    main()
