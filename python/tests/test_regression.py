"""Regression: solve reference decks with CalculiX++ and compare nodal
displacements against stock CalculiX .dat.ref output within tolerance.

Displacements are directly comparable (nodal); CalculiX *EL PRINT stresses are
at integration points, so stress comparison is deferred to a later harness.
"""
import math
import os

import pytest

CCX_TEST = "/home/leonardo/work/CalculiX/test"

# Small, clean linear-static tetrahedral decks with stock-CalculiX .dat.ref output.
# (beam10p is the practical reference case: contact4tet uses S8 shells, circ11p uses
# *USE REFINED MESH, and segmentunsmooth is too large for a fast CI gate.)
REFERENCE_DECKS = ["beam10p"]

# Minimal C3D4 deck: base fixed, pressure 100 on face 3 (nodes 2,3,4). The outward
# face normal is (1,1,1)/sqrt(3) and the area sqrt(3)/2, so the total applied force
# is (-50,-50,-50) and the reactions must sum to (50,50,50).
C3D4_PRESSURE_DECK = """
*NODE
1, 0., 0., 0.
2, 1., 0., 0.
3, 0., 1., 0.
4, 0., 0., 1.
*ELEMENT, TYPE=C3D4, ELSET=EALL
100, 1, 2, 3, 4
*BOUNDARY
1, 1, 3
2, 1, 3
3, 1, 3
*MATERIAL, NAME=EL
*ELASTIC
210000., 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=EL
*STEP
*STATIC
*DLOAD
100, P3, 100.
*END STEP
"""


def parse_ref_displacements(path):
    """Return {node_id: (ux, uy, uz)} from the displacement block of a .dat."""
    disp = {}
    with open(path) as fh:
        lines = fh.readlines()
    i = 0
    while i < len(lines) and "displacements" not in lines[i]:
        i += 1
    for line in lines[i + 1:]:
        parts = line.split()
        if len(parts) == 4:
            try:
                disp[int(parts[0])] = [float(x) for x in parts[1:]]
            except ValueError:
                if disp:
                    break
        elif parts and disp:
            break
    return disp


@pytest.mark.parametrize("deck", REFERENCE_DECKS)
def test_displacements_match_calculix(deck):
    inp = os.path.join(CCX_TEST, deck + ".inp")
    ref_path = os.path.join(CCX_TEST, deck + ".dat.ref")
    if not (os.path.exists(inp) and os.path.exists(ref_path)):
        pytest.skip(f"reference deck {deck} not available")
    import calculixpp

    res = calculixpp.solve(inp)
    ids = res["node_ids"]
    U = res["displacement"]
    ref = parse_ref_displacements(ref_path)
    assert len(ref) > 0

    num = den = 0.0
    max_abs = 0.0
    compared = 0
    for k in range(len(ids)):
        nid = int(ids[k])
        if nid not in ref:
            continue
        compared += 1
        for c in range(3):
            d = float(U[k][c]) - ref[nid][c]
            num += d * d
            den += ref[nid][c] ** 2
            max_abs = max(max_abs, abs(d))

    assert compared == len(ref), f"{deck}: compared {compared} of {len(ref)} nodes"
    rel_l2 = math.sqrt(num / den) if den > 0 else math.sqrt(num)
    assert rel_l2 < 1e-4, f"{deck}: relative L2 displacement error {rel_l2:.3e}"
    assert max_abs < 1e-3, f"{deck}: max abs displacement error {max_abs:.3e}"


def test_c3d4_pressure_equilibrium():
    """C3D4 case: DLOAD pressure, validated analytically via global equilibrium."""
    import calculixpp

    res = calculixpp.solve_text(C3D4_PRESSURE_DECK)
    assert int(res["num_elements"]) == 1
    rf = res["reaction"]
    total = [sum(float(rf[k][c]) for k in range(len(rf))) for c in range(3)]
    for c in range(3):
        assert abs(total[c] - 50.0) < 1e-6, f"reaction sum[{c}]={total[c]}"


def test_exception_propagation():
    import calculixpp

    with pytest.raises(RuntimeError):
        calculixpp.solve_text("*ELEMENT, TYPE=C3D4\n1, 1, 2, 3\n")  # too few nodes
