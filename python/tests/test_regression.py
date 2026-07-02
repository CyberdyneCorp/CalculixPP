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


def _rel_l2_disp(a, b):
    num = den = 0.0
    for k in range(len(a)):
        for c in range(3):
            d = float(a[k][c]) - float(b[k][c])
            num += d * d
            den += float(b[k][c]) ** 2
    return math.sqrt(num / den) if den else math.sqrt(num)


def test_nonlinear_reproduces_linear_single_tet():
    """GATE (spec: nonlinear-solution-control 1.5): the Newton-Raphson driver
    reproduces the linear solve on the single-tet model to rel-L2 < 1e-10 in one
    increment / two iterations."""
    import calculixpp

    lin = calculixpp.solve_text(C3D4_PRESSURE_DECK)
    nl = calculixpp.solve_nonlinear_text(C3D4_PRESSURE_DECK)
    assert nl["converged"] is True
    assert int(nl["newton_increments"]) == 1
    assert int(nl["newton_iterations"]) == 2
    assert _rel_l2_disp(nl["displacement"], lin["displacement"]) < 1e-10

    # Line search on: identical solution.
    ls = calculixpp.solve_nonlinear_text(C3D4_PRESSURE_DECK, line_search=True)
    assert _rel_l2_disp(ls["displacement"], lin["displacement"]) < 1e-10


def test_nonlinear_reproduces_linear_beam10p():
    """GATE: the nonlinear driver reproduces the linear beam10p solve exactly."""
    import calculixpp

    inp = os.path.join(CCX_TEST, "beam10p.inp")
    if not os.path.exists(inp):
        pytest.skip("beam10p deck not available")
    lin = calculixpp.solve(inp)
    nl = calculixpp.solve_nonlinear(inp)
    assert nl["converged"] is True
    assert _rel_l2_disp(nl["displacement"], lin["displacement"]) < 1e-10


def test_nonlinear_controls_and_increments_parse():
    """*CONTROLS, *STATIC increment data, and DIRECT parse and still reproduce the
    linear solve. DIRECT with increment 0.5 takes exactly two increments."""
    import calculixpp

    deck = C3D4_PRESSURE_DECK.replace(
        "*STATIC",
        "*CONTROLS, PARAMETERS=FIELD\n1e-4,,1e-4\n*STATIC, DIRECT\n0.5, 1.0",
    )
    lin = calculixpp.solve_text(C3D4_PRESSURE_DECK)
    nl = calculixpp.solve_nonlinear_text(deck)
    assert int(nl["newton_increments"]) == 2
    assert int(nl["newton_cutbacks"]) == 0
    assert _rel_l2_disp(nl["displacement"], lin["displacement"]) < 1e-10


def test_c3d4_pressure_equilibrium():
    """C3D4 case: DLOAD pressure, validated analytically via global equilibrium."""
    import calculixpp

    res = calculixpp.solve_text(C3D4_PRESSURE_DECK)
    assert int(res["num_elements"]) == 1
    rf = res["reaction"]
    total = [sum(float(rf[k][c]) for k in range(len(rf))) for c in range(3)]
    for c in range(3):
        assert abs(total[c] - 50.0) < 1e-6, f"reaction sum[{c}]={total[c]}"


def test_large_deck_solvers_agree():
    """Scaling + convergence (SciPP#10, v1.2.0): on a larger FE system the sparse
    direct solve and the IC0-preconditioned CG agree. Skipped if the deck is absent."""
    import calculixpp

    inp = os.path.join(CCX_TEST, "segmentunsmooth.inp")
    if not os.path.exists(inp):
        pytest.skip("segmentunsmooth deck not available")
    d = calculixpp.solve(inp, solver="direct")
    c = calculixpp.solve(inp, solver="cg")
    ud, uc = d["displacement"], c["displacement"]
    num = den = 0.0
    for k in range(len(ud)):
        for j in range(3):
            num += (float(ud[k][j]) - float(uc[k][j])) ** 2
            den += float(ud[k][j]) ** 2
    rel = math.sqrt(num / den) if den else 0.0
    assert rel < 1e-5, f"direct vs IC0-CG disagree: relL2={rel:.3e}"


def test_exception_propagation():
    import calculixpp

    with pytest.raises(RuntimeError):
        calculixpp.solve_text("*ELEMENT, TYPE=C3D4\n1, 1, 2, 3\n")  # too few nodes


def test_unknown_solver_raises():
    """An unrecognized SOLVER= on *STATIC must raise (solver not available)."""
    import calculixpp

    deck = C3D4_PRESSURE_DECK.replace("*STATIC", "*STATIC, SOLVER=BOGUS")
    with pytest.raises(RuntimeError):
        calculixpp.solve_text(deck)


def test_available_and_selected_backend():
    """Backend introspection: cpu is always present and is the resolved default."""
    import calculixpp

    backends = calculixpp.available_backends()
    assert "cpu" in backends
    # Empty request resolves to the default (cpu); an unimplemented backend falls
    # back to cpu rather than erroring; an unknown name raises.
    assert calculixpp.selected_backend() == "cpu"
    assert calculixpp.selected_backend("cpu") == "cpu"
    assert calculixpp.selected_backend("cuda") == "cpu"
    with pytest.raises(ValueError):
        calculixpp.selected_backend("bogus")


def test_solve_records_backend():
    """solve_text echoes the backend that ran; an unimplemented one falls back."""
    import calculixpp

    res = calculixpp.solve_text(C3D4_PRESSURE_DECK)
    assert res["backend"] == "cpu"

    # Explicit unimplemented backend must not change results and reports cpu.
    fallback = calculixpp.solve_text(C3D4_PRESSURE_DECK, backend="metal")
    assert fallback["backend"] == "cpu"
    for c in range(3):
        a = sum(float(res["reaction"][k][c]) for k in range(len(res["reaction"])))
        b = sum(float(fallback["reaction"][k][c]) for k in range(len(fallback["reaction"])))
        assert abs(a - b) < 1e-12

    with pytest.raises(ValueError):
        calculixpp.solve_text(C3D4_PRESSURE_DECK, backend="bogus")


def test_summary_without_solving():
    """summary_text reports parsed counts/materials without running a solve."""
    import calculixpp

    s = calculixpp.summary_text(C3D4_PRESSURE_DECK)
    assert int(s["num_nodes"]) == 4
    assert int(s["num_elements"]) == 1
    assert int(s["num_materials"]) == 1
    assert list(s["materials"]) == ["EL"]
    assert s["requested_solver"] == "auto"

    cg = calculixpp.summary_text(
        C3D4_PRESSURE_DECK.replace("*STATIC", "*STATIC, SOLVER=CG")
    )
    assert cg["requested_solver"] == "cg"


def test_summary_file_matches_solve():
    """summary() on the beam10p deck agrees with solve()'s reported counts."""
    import calculixpp

    inp = os.path.join(CCX_TEST, "beam10p.inp")
    if not os.path.exists(inp):
        pytest.skip("beam10p deck not available")
    s = calculixpp.summary(inp)
    res = calculixpp.solve(inp)
    assert int(s["num_nodes"]) == int(res["num_nodes"])
    assert int(s["num_elements"]) == int(res["num_elements"])
    assert int(s["num_materials"]) >= 1


def test_solver_selection_matches_default():
    """SOLVER=SPOOLES (direct family) yields the same result as the default."""
    import calculixpp

    default = calculixpp.solve_text(C3D4_PRESSURE_DECK)
    spooles = calculixpp.solve_text(
        C3D4_PRESSURE_DECK.replace("*STATIC", "*STATIC, SOLVER=SPOOLES")
    )
    for c in range(3):
        d = sum(float(default["reaction"][k][c]) for k in range(len(default["reaction"])))
        s = sum(float(spooles["reaction"][k][c]) for k in range(len(spooles["reaction"])))
        assert abs(d - s) < 1e-9
