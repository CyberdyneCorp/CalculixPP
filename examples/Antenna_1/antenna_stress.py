#!/usr/bin/env python3
"""Linear-static stress of a HOLLOW cell-tower antenna mast under self-weight + payload.

Reads the `Antenna.obj` art asset (a tall, tapered, red/white-banded tubular mast),
builds a HOLLOW finite-element mesh (a wall, not a solid slug), clamps the base, applies
the mast's own weight (gravity) plus a 140 kg antenna payload on top, solves with
CalculiX++, and visualises the von Mises stress.

Why hollow + self-weight: a real mast is a thin-walled tube whose own weight is a
significant load. A *solid* fill would be ~15x too heavy and mask that. We hollow the
voxel model (solid minus an eroded core) so the geometry is a tube.

Wall thickness / mass: the thinnest robust voxel wall is 1 cell (here 150 mm) — still far
thicker than a real steel wall (~15 mm), which voxels can't reach. So we don't trust the
modelled volume for mass; instead we CALIBRATE the density so the model weighs what the
real hollow tube would (outer surface area x 15 mm wall x steel density). Gravity then
produces a realistic self-weight.

Solver: a thin-walled shell is ill-conditioned for CG, so we force the sparse DIRECT
solver (the structure is slender -> narrow bandwidth -> direct is efficient).

Units: CalculiX consistent set (N, mm, MPa, tonne). The model is ~44.7 m tall.

Run:
    PYTHONPATH=/path/to/CalculixPP/build/python python3 antenna_stress.py [--preview]
"""

# Import calculixpp before numpy (see femkit for the libstdc++ rationale).
import calculixpp as cx

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import femkit  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OBJ_PATH = os.path.join(HERE, "Antenna.obj")

SCALE_TO_MM = 1000.0           # the OBJ is modelled in metres
# Voxel edge -> mesh resolution (smaller = sharper). The hollow shell has ~1/voxel^2
# hexes and the sparse DIRECT solve grows steeply: 150 mm ~135k DOF/~18 s,
# 100 mm ~322k DOF/~2 min, 80 mm ~600k DOF/~8 min, 50 mm ~1.4M DOF (~40 GB, impractical
# single-threaded; CG can't be used — a thin wall is ill-conditioned). 100 mm is the
# practical sweet spot for a sharp mesh here.
VOXEL_MM = 100.0
WALL_VOXELS = 1                # hollow-wall thickness in voxels
WALL_REAL_MM = 15.0            # realistic steel wall thickness, for mass calibration
E, NU = 210000.0, 0.3          # structural steel
STEEL_DENSITY = 7.85e-9        # tonne/mm^3
GRAVITY = 9810.0               # mm/s^2 (acts along -Z, the mast's height)
PAYLOAD_MASS = 0.30            # tonne (= 300 kg antenna array on top)

# Wind (applied with --wind): steady drag on the mast blowing along +X.
WIND_SPEED_KMH = 200.0
AIR_DENSITY = 1.225e-12        # tonne/mm^3 (1.225 kg/m^3)
DRAG_CD = 0.7                  # drag coefficient of a circular cylinder (high Re)


def wind_cloads(coords, q_cd, voxel):
    """Horizontal (+X) wind line-load as nodal forces, distributed by height.

    Drag force per height strip = q*Cd * projected_width * strip_height, where the
    projected width (seen by a +X wind) is the Y-extent of that cross-section. Each
    strip's force is shared over its nodes -> the net section drag acts at that height,
    and the tension(windward)/compression(leeward) split emerges from the bending.
    """
    levels = {}
    for i in range(len(coords)):
        levels.setdefault(round(coords[i, 2] / voxel), []).append(i)
    out = []
    for idxs in levels.values():
        width = np.ptp(coords[idxs, 1])                 # Y-extent ⟂ to the +X wind
        if width <= 0:
            continue
        f = q_cd * width * voxel / len(idxs)            # per-node force, +X (dof 1)
        out += [(i + 1, 1, f) for i in idxs]
    return out


def build_mesh(tri):
    """Hollow solid-element mesh of the mast + a density calibrated to a real tube mass."""
    occ, origin = femkit.voxelize_solid_fill(tri, VOXEL_MM, axis=2)   # fill per Z-slice
    hollow = femkit.largest_component(femkit.hollow_shell(occ, WALL_VOXELS))
    mesh = femkit.VoxelMesh(hollow, origin, VOXEL_MM)
    # Calibrate density so the (thick-walled) model weighs like the real thin-walled tube.
    # The OBJ is a closed double-walled tube (outer + inner skins), so its total area is
    # ~2x the outer wall -> use half for the real wall mass (outer_area x wall x density).
    real_mass = 0.5 * femkit.surface_area(tri) * WALL_REAL_MM * STEEL_DENSITY   # tonne
    model_volume = hollow.sum() * VOXEL_MM ** 3                            # mm^3
    density = real_mass / model_volume
    return mesh, hollow.sum(), real_mass, density


def main():
    tri = femkit.read_obj_triangles(OBJ_PATH)[:, :, [0, 2, 1]] * SCALE_TO_MM  # Y-up -> Z-up
    mesh, n_hex, real_mass, density = build_mesh(tri)
    z = mesh.coords[:, 2]
    print(f"Hollow antenna mast: {n_hex} hexes ({VOXEL_MM:.0f} mm voxels, "
          f"{WALL_VOXELS}-voxel wall), {len(mesh.coords)} nodes, height {np.ptp(z) / 1000:.1f} m")
    print(f"  calibrated to a real {WALL_REAL_MM:.0f} mm-wall tube: "
          f"mass {real_mass * 1000:.0f} kg (density {density:.3e} t/mm^3)")

    if "--preview" in sys.argv:
        femkit.preview_geometry(mesh.coords, mesh.elements,
                                os.path.join(HERE, "antenna_mesh.png"),
                                f"Antenna mast hollow mesh — {n_hex} hexes",
                                views=((8, -70), (8, 20), (8, 110), (75, -90)))
        print("Wrote antenna_mesh.png")
        return

    fixed = mesh.nodes_on_plane(axis=2, value=z.min())       # clamp the base ring
    top = mesh.nodes_on_plane(axis=2, value=z.max())         # payload on the top ring
    self_weight = real_mass * GRAVITY
    cloads = [(n, 3, -PAYLOAD_MASS * GRAVITY / len(top)) for n in top]
    print(f"  clamped base nodes : {len(fixed)}")
    print(f"  loads: self-weight {self_weight:.0f} N + payload "
          f"{PAYLOAD_MASS * GRAVITY:.0f} N (down)")

    wind = "--wind" in sys.argv
    if wind:
        v = WIND_SPEED_KMH * 1e6 / 3600.0                    # km/h -> mm/s
        q_cd = 0.5 * AIR_DENSITY * v ** 2 * DRAG_CD           # N/mm^2
        wind_loads = wind_cloads(mesh.coords, q_cd, VOXEL_MM)
        cloads += wind_loads
        wind_force = sum(f for _, _, f in wind_loads)
        print(f"         + {WIND_SPEED_KMH:.0f} km/h wind (q·Cd={q_cd * 1e3:.2f} kPa) "
              f"= {wind_force / 1000:.0f} kN sideways (+X)")

    deck = femkit.static_deck(mesh.coords, mesh.elements, fixed, cloads,
                              density=density, gravity=(GRAVITY, (0, 0, -1)))
    # Thin-walled shell is ill-conditioned for CG; the slender mast has a narrow
    # bandwidth, so the sparse direct solver is both robust and efficient here.
    result = cx.solve_text(deck, solver="direct")

    order = femkit.align_to_ids(result, len(mesh.coords))
    coords = result["node_coords"][order]
    disp = result["displacement"][order]
    vm = femkit.von_mises(result["stress"][order])

    print("\nResults")
    print(f"  max lateral sway (+X)   : {disp[:, 0].max():.2f} mm")
    print(f"  max vertical shortening : {disp[:, 2].min():+.4f} mm")
    print(f"  max von Mises stress    : {vm.max():.3f} MPa")
    print(f"  horizontal reaction (X) : {result['reaction'][:, 0].sum():+.0f} N")

    tag = "wind_" if wind else ""
    subtitle = (f"self-weight + {PAYLOAD_MASS * 1000:.0f} kg"
                + (f" + {WIND_SPEED_KMH:.0f} km/h wind" if wind else " payload on top"))
    title = f"Antenna mast (hollow) — von Mises\n{subtitle}"
    png = os.path.join(HERE, f"antenna_{tag}von_mises.png")
    vtk = os.path.join(HERE, f"antenna_{tag}stress.vtk")
    femkit.render_stress(coords, mesh.elements, vm, disp, png, title,
                         max_scale=(800 if wind else 3000), elev=8, azim=-65,
                         figsize=(7, 11))
    femkit.write_vtk(vtk, coords, mesh.elements, vm, disp)
    print(f"\nWrote {png}\n      {vtk}")


if __name__ == "__main__":
    main()
