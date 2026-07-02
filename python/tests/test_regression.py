"""Regression: solve reference decks with CalculiX++ and compare nodal
displacements against stock CalculiX .dat.ref output within tolerance.

Displacements are directly comparable (nodal); CalculiX *EL PRINT stresses are
at integration points, so stress comparison is deferred to a later harness.
"""
import math
import os

import pytest

CCX_TEST = "/home/leonardo/work/CalculiX/test"
BEAM_INP = os.path.join(CCX_TEST, "beam10p.inp")
BEAM_REF = os.path.join(CCX_TEST, "beam10p.dat.ref")

pytestmark = pytest.mark.skipif(
    not (os.path.exists(BEAM_INP) and os.path.exists(BEAM_REF)),
    reason="CalculiX reference deck not available",
)


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


def test_beam10p_displacements():
    import calculixpp

    res = calculixpp.solve(BEAM_INP)
    assert int(res["num_nodes"]) == 90
    assert int(res["num_elements"]) == 31

    ids = res["node_ids"]
    U = res["displacement"]
    ref = parse_ref_displacements(BEAM_REF)
    assert len(ref) == 90, f"parsed {len(ref)} reference nodes"

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

    assert compared == 90
    rel_l2 = math.sqrt(num / den) if den > 0 else math.sqrt(num)
    assert rel_l2 < 1e-4, f"relative L2 displacement error {rel_l2:.3e}"
    assert max_abs < 1e-4, f"max abs displacement error {max_abs:.3e}"


def test_exception_propagation():
    import calculixpp

    with pytest.raises(RuntimeError):
        calculixpp.solve_text("*ELEMENT, TYPE=C3D4\n1, 1, 2, 3\n")  # too few nodes
