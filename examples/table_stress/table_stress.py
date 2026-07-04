#!/usr/bin/env python3
"""Linear-static stress analysis of a 3D table under a heavy central object.

Builds a steel table (tabletop + four legs) as a conforming C3D8 hexahedral mesh,
clamps the leg bottoms, places a 100 kg object at the centre of the tabletop, solves
the linear-static step with CalculiX++, and visualises the von Mises stress on the
deformed body.

Units are the standard CalculiX consistent set: N, mm, MPa, tonne (so density is in
tonne/mm^3 and gravity is 9810 mm/s^2). A 100 kg mass is 0.1 tonne, i.e. a 981 N weight.

Run:
    PYTHONPATH=/path/to/CalculixPP/build/python python3 table_stress.py

Outputs (next to this script):
    table_von_mises.png   matplotlib render of the deformed table coloured by stress
    table_stress.vtk      legacy VTK for ParaView / any VTK viewer
"""

# IMPORTANT: import calculixpp BEFORE numpy/matplotlib. The compiled module links the
# system libstdc++; importing numpy first can load an older libstdc++ (e.g. from a
# conda env) and shadow the newer symbols the module needs.
import calculixpp as cx

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import femkit  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Geometry / material (mm, N, MPa, tonne)
# ---------------------------------------------------------------------------
LENGTH, WIDTH, THICK = 800.0, 500.0, 40.0   # tabletop size
LEG_HEIGHT = 710.0                          # clearance under the tabletop
LEG_SIZE = 50.0                             # square leg cross-section
LEG_INSET = 40.0                            # leg inset from the tabletop edge

E, NU, DENSITY = 210000.0, 0.3, 7.85e-9     # structural steel
GRAVITY = 9810.0                            # mm/s^2
OBJECT_MASS = 0.1                           # tonne (= 100 kg)
OBJECT_HALF_FOOTPRINT = 120.0               # object rests on a ~240 x 240 mm patch


def _sorted_union(*sequences, tol=1e-6):
    """Merge coordinate lists into one sorted list without near-duplicates."""
    merged = []
    for seq in sequences:
        for value in seq:
            if not any(abs(value - m) < tol for m in merged):
                merged.append(float(value))
    return sorted(merged)


class TableMesh:
    """A conforming hex mesh of the table, with the CalculiX deck it produces."""

    def __init__(self):
        self._nodes = {}          # (x, y, z) rounded -> 1-based id
        self.elements = []        # list of 8 node ids (C3D8 ordering)
        self._build()

    def _node(self, x, y, z):
        key = (round(x, 4), round(y, 4), round(z, 4))
        if key not in self._nodes:
            self._nodes[key] = len(self._nodes) + 1
        return self._nodes[key]

    def _hex(self, x0, x1, y0, y1, z0, z1):
        self.elements.append([
            self._node(x0, y0, z0), self._node(x1, y0, z0),
            self._node(x1, y1, z0), self._node(x0, y1, z0),
            self._node(x0, y0, z1), self._node(x1, y0, z1),
            self._node(x1, y1, z1), self._node(x0, y1, z1),
        ])

    def _build(self):
        leg_x = [(LEG_INSET, LEG_INSET + LEG_SIZE),
                 (LENGTH - LEG_INSET - LEG_SIZE, LENGTH - LEG_INSET)]
        leg_y = [(LEG_INSET, LEG_INSET + LEG_SIZE),
                 (WIDTH - LEG_INSET - LEG_SIZE, WIDTH - LEG_INSET)]
        # Grid lines: the tabletop grid must contain the leg footprints so leg and
        # tabletop nodes coincide at the interface (a conforming mesh).
        xs = _sorted_union(np.linspace(0, LENGTH, 9), [v for p in leg_x for v in p])
        ys = _sorted_union(np.linspace(0, WIDTH, 7), [v for p in leg_y for v in p])
        z_top = [LEG_HEIGHT, LEG_HEIGHT + THICK / 2, LEG_HEIGHT + THICK]
        z_leg = list(np.linspace(0, LEG_HEIGHT, 8))

        cells = lambda c: list(zip(c[:-1], c[1:]))
        # Tabletop: full footprint over the top z-layers.
        for x0, x1 in cells(xs):
            for y0, y1 in cells(ys):
                for z0, z1 in cells(z_top):
                    self._hex(x0, x1, y0, y1, z0, z1)
        # Legs: four corner columns sharing the tabletop's bottom nodes.
        span = lambda c, a, b: [v for v in c if a - 1e-6 <= v <= b + 1e-6]
        for xa, xb in leg_x:
            for ya, yb in leg_y:
                for x0, x1 in cells(span(xs, xa, xb)):
                    for y0, y1 in cells(span(ys, ya, yb)):
                        for z0, z1 in cells(z_leg):
                            self._hex(x0, x1, y0, y1, z0, z1)

    @property
    def coords(self):
        """(N, 3) node coordinates indexed by (id - 1)."""
        out = np.zeros((len(self._nodes), 3))
        for key, nid in self._nodes.items():
            out[nid - 1] = key
        return out

    def ids_on_plane(self, axis, value, tol=1e-6):
        return [nid for key, nid in self._nodes.items() if abs(key[axis] - value) < tol]

    def top_center_ids(self):
        z_top = LEG_HEIGHT + THICK
        cx_, cy_ = LENGTH / 2, WIDTH / 2
        return [nid for key, nid in self._nodes.items()
                if abs(key[2] - z_top) < 1e-6
                and abs(key[0] - cx_) < OBJECT_HALF_FOOTPRINT
                and abs(key[1] - cy_) < OBJECT_HALF_FOOTPRINT]

    def to_inp(self):
        coords = self.coords
        lines = ["*NODE"]
        lines += [f"{i + 1},{coords[i, 0]},{coords[i, 1]},{coords[i, 2]}"
                  for i in range(len(coords))]
        lines.append("*ELEMENT,TYPE=C3D8,ELSET=EALL")
        lines += [f"{e},{','.join(map(str, nodes))}"
                  for e, nodes in enumerate(self.elements, 1)]

        fixed = self.ids_on_plane(axis=2, value=0.0)         # clamp leg bottoms
        lines.append("*BOUNDARY")
        lines += [f"{n},1,3" for n in fixed]

        load_nodes = self.top_center_ids()                   # 100 kg at the centre
        force = -OBJECT_MASS * GRAVITY / len(load_nodes)     # split equally, downward
        lines.append("*CLOAD")
        lines += [f"{n},3,{force}" for n in load_nodes]

        lines += [
            "*MATERIAL,NAME=STEEL",
            "*ELASTIC", f"{E},{NU}",
            "*DENSITY", f"{DENSITY}",
            "*SOLID SECTION,ELSET=EALL,MATERIAL=STEEL",
            "*STEP", "*STATIC", "*END STEP",
        ]
        return "\n".join(lines), fixed, load_nodes


def main():
    table = TableMesh()
    deck, fixed, load_nodes = table.to_inp()
    print(f"Table mesh: {len(table.coords)} nodes, {len(table.elements)} C3D8 elements")
    print(f"  clamped leg-bottom nodes : {len(fixed)}")
    print(f"  loaded centre-top nodes  : {len(load_nodes)}  "
          f"(100 kg = {OBJECT_MASS * GRAVITY:.0f} N total, downward)")

    result = cx.solve_text(deck)

    # Align returned fields (solver node order) with our element node ids.
    order = femkit.align_to_ids(result, len(table.coords))
    coords = result["node_coords"][order]
    disp = result["displacement"][order]
    vm = femkit.von_mises(result["stress"][order])

    print("\nResults")
    print(f"  max vertical deflection : {disp[:, 2].min():+.4f} mm")
    print(f"  max von Mises stress    : {vm.max():.2f} MPa "
          f"(steel yield ~250 MPa -> safety factor ~{250 / max(vm.max(), 1e-9):.0f})")
    print(f"  vertical reaction sum   : {result['reaction'][:, 2].sum():+.1f} N "
          f"(balances the {OBJECT_MASS * GRAVITY:.0f} N weight)")

    png = os.path.join(HERE, "table_von_mises.png")
    vtk = os.path.join(HERE, "table_stress.vtk")
    femkit.render_stress(coords, table.elements, vm, disp, png,
                         "Table von Mises stress under 100 kg at centre")
    femkit.write_vtk(vtk, coords, table.elements, vm, disp)
    print(f"\nWrote {png}\n      {vtk}")


if __name__ == "__main__":
    main()
