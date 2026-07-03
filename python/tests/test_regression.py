"""Regression: solve reference decks with CalculiX++ and compare nodal
displacements against stock CalculiX .dat.ref output within tolerance.

Displacements are directly comparable (nodal); CalculiX *EL PRINT stresses are
at integration points, so stress comparison is deferred to a later harness.
"""
import math
import os

import pytest

CCX_TEST = "/home/leonardo/work/CalculiX/test"

# Small, clean linear-static decks with stock-CalculiX .dat.ref output.
#   beam10p : C3D10 quadratic-tet cantilever (Phase-1 reference case).
#   beam8p  : C3D8 linear-hex cantilever  (validates the Phase-2 hex family).
#   beam20p : C3D20 quadratic-hex cantilever (validates the Phase-2 hex family).
# (contact4tet uses S8 shells, circ11p uses *USE REFINED MESH, and segmentunsmooth
# is too large for a fast CI gate.)
REFERENCE_DECKS = ["beam10p", "beam8p", "beam20p"]

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


def test_equation_constraint_matches_calculix():
    """GATE (spec: constraints 5.1): the achtel2 deck (C3D20R cube with 2-term
    *EQUATION cards tying node 178 to node 78 in all 3 DOFs) reproduces stock
    CalculiX .dat.ref displacements within tolerance — validates dependent-DOF
    elimination against a real reference deck."""
    inp = os.path.join(CCX_TEST, "achtel2.inp")
    ref_path = os.path.join(CCX_TEST, "achtel2.dat.ref")
    if not (os.path.exists(inp) and os.path.exists(ref_path)):
        pytest.skip("achtel2 reference deck not available")
    import calculixpp

    res = calculixpp.solve(inp)
    ids = res["node_ids"]
    U = res["displacement"]
    id2idx = {int(n): i for i, n in enumerate(ids)}
    ref = parse_ref_displacements(ref_path)
    assert len(ref) > 0

    num = den = 0.0
    for nid, ru in ref.items():
        if nid not in id2idx:
            continue
        u = U[id2idx[nid]]
        for c in range(3):
            d = float(u[c]) - ru[c]
            num += d * d
            den += ru[c] ** 2
    rel_l2 = math.sqrt(num / den) if den > 0 else math.sqrt(num)
    assert rel_l2 < 1e-4, f"achtel2: relative L2 displacement error {rel_l2:.3e}"

    # The eliminated dependent DOFs (node 178) equal their masters (node 78).
    for c in range(3):
        assert abs(float(U[id2idx[178]][c]) - float(U[id2idx[78]][c])) < 1e-10


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


# Single C3D8 unit cube, uniaxial STRESS along z: symmetry planes fix the three back
# faces, a *CLOAD along z on the top face drives it past yield. Analytic hardening:
#   sigma = sy0 + H*eps_p,  eps_p = (sigma - sy0)/H,  eps = sigma/E + eps_p.
# With sy0=800, H=(960-800)/0.02=8000, a total axial force 900 (area 1) gives
# sigma=900 exactly. (No clean *PLASTIC *STATIC C3D4/C3D10 .dat.ref exists in the
# reference corpus — those decks use NLGEOM, multi-step, orthotropic elasticity, or
# Johnson-Cook — so plasticity is validated against the analytic uniaxial solution.)
_PLASTIC_CUBE = """
*NODE
1, 0., 0., 0.
2, 1., 0., 0.
3, 1., 1., 0.
4, 0., 1., 0.
5, 0., 0., 1.
6, 1., 0., 1.
7, 1., 1., 1.
8, 0., 1., 1.
*ELEMENT, TYPE=C3D8, ELSET=EALL
1, 1, 2, 3, 4, 5, 6, 7, 8
*BOUNDARY
1, 1, 1
4, 1, 1
5, 1, 1
8, 1, 1
1, 2, 2
2, 2, 2
5, 2, 2
6, 2, 2
1, 3, 3
2, 3, 3
3, 3, 3
4, 3, 3
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*PLASTIC
800., 0.0
960., 0.02
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*STEP
*STATIC
0.25, 1.0
*CLOAD
5, 3, 225.
6, 3, 225.
7, 3, 225.
8, 3, 225.
*END STEP
"""


def test_plasticity_uniaxial_hardening():
    """Validation (a): a single C3D8 past yield follows the analytic uniaxial
    hardening curve. Also checks that solve() auto-routes a *PLASTIC deck to the
    Newton driver (converged/newton_* keys present)."""
    import calculixpp

    E, sy0, H = 210000.0, 800.0, 8000.0
    sigma = 900.0  # total force 900 over unit area
    eps_p = (sigma - sy0) / H
    eps = sigma / E + eps_p

    # solve() must auto-dispatch to the nonlinear driver for a plastic model.
    res = calculixpp.solve_text(_PLASTIC_CUBE)
    assert res["converged"] is True
    assert int(res["newton_increments"]) >= 1

    S = res["stress"]
    szz = sum(float(S[k][2]) for k in range(len(S))) / len(S)
    assert abs(szz - sigma) < 1e-3 * sigma, f"axial stress {szz} vs {sigma}"

    U = res["displacement"]
    uz = max(float(U[k][2]) for k in range(len(U)))
    assert abs(uz - eps) < 1e-3 * eps, f"axial strain {uz} vs {eps}"


def test_plasticity_matches_nonlinear_call():
    """solve() (auto-dispatch) and solve_nonlinear() give the same plastic result."""
    import calculixpp

    a = calculixpp.solve_text(_PLASTIC_CUBE)
    b = calculixpp.solve_nonlinear_text(_PLASTIC_CUBE)
    assert _rel_l2_disp(a["displacement"], b["displacement"]) < 1e-12


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
