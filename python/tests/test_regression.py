"""Phase-2 regression corpus (tasks 6.3 / 6.4).

Two kinds of validation, both parametrized over a curated corpus:

* REFERENCE decks — clean, in-scope stock-CalculiX decks with a matching
  ``.dat.ref``. We solve them with CalculiX++ and compare nodal displacements
  against the reference within a per-deck relative-L2 tolerance. Displacements
  are directly comparable (nodal); CalculiX ``*EL PRINT`` stresses are at
  integration points, so we compare displacements only.

* ANALYTICAL decks — physics for which the reference corpus has no clean,
  in-scope deck (J2 plasticity, ``*EQUATION`` / MPC, body loads driven purely
  by equilibrium). These carry an embedded deck plus a closed-form expected
  value checked directly.

Per-deck tolerances live in the manifests below (REFERENCE_CORPUS /
ANALYTICAL_CORPUS). Reference files are skipped gracefully when absent so the
suite still runs on a checkout without the CalculiX ``test/`` tree.

Corpus scope (why decks are IN or OUT): only linear ``*STATIC`` (no NLGEOM),
single-step, solid C3D4/6/8/10/15/20 (+ reduced) elements. Excluded families we
verified against the reference tree: shells (S*) / beams (B*) / plane
(CPS/CPE/CAX), heat / dynamic / frequency / contact, ``*REFINE MESH`` /
``*USE REFINED MESH``, multi-step, and anisotropic ``*ELASTIC,TYPE=`` — none are
implemented in Phase 2. Very large decks (segmentunsmooth) are left out to keep
the gate fast.
"""
import math
import os

import pytest

CCX_TEST = "/home/leonardo/work/CalculiX/test"

# --- Reference corpus: (deck, rel-L2 tolerance, what it validates) ------------
# Every entry is a stock-CalculiX deck with a .dat.ref we reproduce. Tolerances
# are set a comfortable margin above the observed rel-L2 (recorded in the note)
# so the gate catches regressions without being brittle to last-bit noise.
#
#   deck          tol       validates
#   beam10p       1e-4      C3D10 quadratic tet   (obs 5.36e-8) — Phase-1 base
#   beam8p        1e-4      C3D8 linear hex       (obs 5.71e-8) — hex family
#   beam20p       1e-4      C3D20 quadratic hex   (obs 5.76e-8) — hex family
#   achtelg       1e-5      C3D20R + GRAV body load     (obs 1.54e-7)
#   achtelc       1e-5      C3D20R + CENTRIF body load  (obs 1.61e-7)
#   achtel2       1e-4      C3D20R + 2-term *EQUATION   (obs 1.27e-7)
#   achtel9       1e-4      C3D20R + 9-term *EQUATION   (obs 1.45e-7)
#   beam8p_mpc    1e-4      C3D8 + *EQUATION on a beam  (obs 1.15e-8, 1 node ref)
REFERENCE_CORPUS = [
    ("beam10p", 1e-4, "C3D10 quadratic-tet cantilever"),
    ("beam8p", 1e-4, "C3D8 linear-hex cantilever"),
    ("beam20p", 1e-4, "C3D20 quadratic-hex cantilever"),
    ("achtelg", 1e-5, "C3D20R cube, GRAV body load"),
    ("achtelc", 1e-5, "C3D20R cube, CENTRIF body load"),
    ("achtel2", 1e-4, "C3D20R cube, 2-term *EQUATION"),
    ("achtel9", 1e-4, "C3D20R cube, 9-term *EQUATION"),
    ("beam8p_mpc", 1e-4, "C3D8 beam with *EQUATION"),
]

# Back-compat alias kept for any external reference to the old name.
REFERENCE_DECKS = [deck for deck, _tol, _desc in REFERENCE_CORPUS]

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


def _reference_rel_l2(deck):
    """Solve a reference deck and return (rel_l2, max_abs, compared, total, U, idx).

    Skips the test gracefully if the deck or its .dat.ref is absent so the suite
    runs on a checkout without the CalculiX test tree.
    """
    inp = os.path.join(CCX_TEST, deck + ".inp")
    ref_path = os.path.join(CCX_TEST, deck + ".dat.ref")
    if not (os.path.exists(inp) and os.path.exists(ref_path)):
        pytest.skip(f"reference deck {deck} not available")
    import calculixpp

    res = calculixpp.solve(inp)
    ids = res["node_ids"]
    U = res["displacement"]
    id2idx = {int(n): i for i, n in enumerate(ids)}
    ref = parse_ref_displacements(ref_path)
    assert len(ref) > 0, f"{deck}: empty reference displacement block"

    # Compare every model node that the reference reports. A model node missing
    # from the reference block would be an unexpected coverage gap (guarded below).
    # The reference may additionally list ids our model does not define (CalculiX
    # can print a padded NSET range) — those are simply not compared.
    num = den = 0.0
    max_abs = 0.0
    compared = 0
    for nid in id2idx:
        if nid not in ref:
            continue
        compared += 1
        u = U[id2idx[nid]]
        ru = ref[nid]
        for c in range(3):
            d = float(u[c]) - ru[c]
            num += d * d
            den += ru[c] ** 2
            max_abs = max(max_abs, abs(d))
    rel_l2 = math.sqrt(num / den) if den > 0 else math.sqrt(num)
    return rel_l2, max_abs, compared, len(id2idx), U, id2idx


@pytest.mark.parametrize(
    "deck, tol, desc",
    REFERENCE_CORPUS,
    ids=[d for d, _t, _s in REFERENCE_CORPUS],
)
def test_reference_corpus_matches_calculix(deck, tol, desc):
    """(6.4) Each reference deck reproduces stock-CalculiX .dat.ref displacements
    within its per-deck rel-L2 tolerance (manifest: REFERENCE_CORPUS)."""
    rel_l2, _max_abs, compared, _num_model, _U, _idx = _reference_rel_l2(deck)
    # Guard against a silent no-op: a deck whose *NODE PRINT selects only a subset
    # (e.g. beam8p_mpc prints a single node) is fine, but we must compare at least
    # one node against the reference.
    assert compared > 0, f"{deck}: no reference nodes matched the model"
    assert rel_l2 < tol, f"{deck} ({desc}): rel-L2 {rel_l2:.3e} >= tol {tol:.0e}"


def test_equation_constraint_matches_calculix():
    """GATE (spec: constraints 5.1): the achtel2 deck (C3D20R cube with 2-term
    *EQUATION cards tying node 178 to node 78 in all 3 DOFs) reproduces stock
    CalculiX .dat.ref displacements within tolerance — validates dependent-DOF
    elimination against a real reference deck."""
    rel_l2, _max_abs, _compared, _total, U, id2idx = _reference_rel_l2("achtel2")
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


# --- Analytical corpus (6.3): physics with no clean in-scope reference deck ----
# Each analytical deck pairs an embedded input with a closed-form expected value.
# These cover the Phase-2 physics the reference tree cannot exercise cleanly:
#   plasticity      : reference *PLASTIC *STATIC decks all use NLGEOM / multi-step
#                     / orthotropic elasticity / Johnson-Cook (out of scope), so
#                     J2 hardening is checked against the analytic uniaxial curve.
#   body load       : gravity global equilibrium (reactions balance the weight).
#   *EQUATION       : a tied DOF must equal its master (analytic identity), on top
#                     of the reference *EQUATION decks (achtel2 / achtel9).
# (Plasticity + C3D4 pressure equilibrium have their own dedicated tests above /
# below; this manifest documents the full analytical set for the corpus report.)

# Single C3D8 unit cube under self-weight along -z, base fixed. Global vertical
# equilibrium: sum of z-reactions = rho * g * V (V = 1). With rho = 7.85e-9 and
# g = 9810 (consistent tonne-mm-s units), weight = 7.6979e-5.
_GRAVITY_CUBE = """
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
1, 1, 3
2, 1, 3
3, 1, 3
4, 1, 3
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*DENSITY
7.85e-9
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*STEP
*STATIC
*DLOAD
EALL, GRAV, 9810., 0., 0., -1.
*END STEP
"""


def test_gravity_body_load_equilibrium():
    """ANALYTICAL (body load): the z-reactions balance the cube's self-weight
    rho*g*V exactly (global vertical equilibrium)."""
    import calculixpp

    rho, g, vol = 7.85e-9, 9810.0, 1.0
    weight = rho * g * vol  # total downward body force

    res = calculixpp.solve_text(_GRAVITY_CUBE)
    rf = res["reaction"]
    rz = sum(float(rf[k][2]) for k in range(len(rf)))
    assert abs(rz - weight) < 1e-6 * weight, f"reaction_z {rz} vs weight {weight}"


# Two C3D8 cubes side by side sharing no nodes; an *EQUATION ties the free block's
# loaded corner to the anchored block so it cannot rigid-body drift. The tied
# (dependent) DOF must equal its master DOF to machine precision — an analytic
# identity independent of the stiffness. Left block (nodes 1-8) fully fixed at
# x=0..1; right block (nodes 11-18) fixed only at its shared face, its far-top
# node 17 (dof 1) tied to left node 7 (dof 1) via *EQUATION.
_EQUATION_TIE = """
*NODE
1, 0., 0., 0.
2, 1., 0., 0.
3, 1., 1., 0.
4, 0., 1., 0.
5, 0., 0., 1.
6, 1., 0., 1.
7, 1., 1., 1.
8, 0., 1., 1.
11, 2., 0., 0.
12, 3., 0., 0.
13, 3., 1., 0.
14, 2., 1., 0.
15, 2., 0., 1.
16, 3., 0., 1.
17, 3., 1., 1.
18, 2., 1., 1.
*ELEMENT, TYPE=C3D8, ELSET=EL
1, 1, 2, 3, 4, 5, 6, 7, 8
2, 11, 12, 13, 14, 15, 16, 17, 18
*BOUNDARY
1, 1, 3
4, 1, 3
5, 1, 3
8, 1, 3
11, 1, 3
14, 1, 3
15, 1, 3
18, 1, 3
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*SOLID SECTION, ELSET=EL, MATERIAL=STEEL
*EQUATION
2
17, 1, 1.0, 7, 1, -1.0
*STEP
*STATIC
*CLOAD
7, 1, 100.
*END STEP
"""


def test_equation_analytical_tie():
    """ANALYTICAL (*EQUATION): the dependent DOF (node 17, dof 1) exactly equals
    its master (node 7, dof 1) — an identity the constraint elimination must hold
    regardless of stiffness. Complements the achtel2 / achtel9 reference checks."""
    import calculixpp

    res = calculixpp.solve_text(_EQUATION_TIE)
    ids = res["node_ids"]
    U = res["displacement"]
    idx = {int(n): i for i, n in enumerate(ids)}
    u_master = float(U[idx[7]][0])
    u_dep = float(U[idx[17]][0])
    assert abs(u_dep - u_master) < 1e-10, f"dep {u_dep} != master {u_master}"
    assert abs(u_master) > 1e-9, "load should produce a nonzero master displacement"


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


# --- Phase-2 workstream 6: parser sweep + bindings introspection + API parity ---

# A hex + connector + constraint + plasticity + amplitude + body-load deck: one
# deck exercising the Phase-2 introspection surface for summary().
_PHASE2_SUMMARY_DECK = """
*AMPLITUDE, NAME=RAMP, DEFINITION=TABULAR
0., 0., 1., 1.
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
1, 1,2,3,4,5,6,7,8
*ELEMENT, TYPE=SPRINGA, ELSET=ESPR
50, 7, 8
*NSET, NSET=TOP
5,6,7,8
*SPRING, ELSET=ESPR
1000.
*EQUATION
2
6, 1, 1.0, 7, 1, -1.0
*MPC
BEAM, 5, 6
*RIGID BODY, NSET=TOP, REF NODE=5
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*DENSITY
7.8e-9
*PLASTIC
800., 0.
960., 0.02
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*STEP
*STATIC
*DLOAD
EALL, GRAV, 9810., 0., 0., -1.
*END STEP
"""


def test_summary_phase2_introspection():
    """(6.2) summary() reflects the Phase-2 capabilities the deck exercises:
    element-type counts, has_plasticity / has_nonlinear_material, constraint
    counts, and connector/amplitude/body-load counts — without solving."""
    import calculixpp

    s = calculixpp.summary_text(_PHASE2_SUMMARY_DECK)
    assert dict(s["element_type_counts"]) == {"C3D8": 1}
    assert s["has_plasticity"] is True
    assert s["has_nonlinear_material"] is True  # solve() will auto-route to Newton
    assert int(s["num_constraints"]) == 3  # 1 equation + 1 mpc + 1 rigid body
    assert s["has_constraints"] is True
    assert int(s["num_equations"]) == 1
    assert int(s["num_mpcs"]) == 1
    assert int(s["num_rigid_bodies"]) == 1
    assert int(s["num_couplings"]) == 0
    assert int(s["num_ties"]) == 0
    assert int(s["num_springs"]) == 1
    assert int(s["num_amplitudes"]) == 1
    assert int(s["num_body_loads"]) == 1


def test_summary_elastic_deck_is_linear():
    """A purely elastic deck reports no nonlinear material and no constraints, so
    solve() keeps the linear path (has_nonlinear_material drives auto-dispatch)."""
    import calculixpp

    s = calculixpp.summary_text(C3D4_PRESSURE_DECK)
    assert s["has_nonlinear_material"] is False
    assert s["has_plasticity"] is False
    assert s["has_constraints"] is False
    assert int(s["num_constraints"]) == 0
    assert dict(s["element_type_counts"]) == {"C3D4": 1}


@pytest.mark.parametrize(
    "card, needle",
    [
        ("*HYPERFOAM, N=1\n1.,2.,0.1", "hyperelastic foam"),
        ("*CREEP, LAW=NORTON\n1e-10, 5., 0.", "creep"),
        ("*VISCO", "viscoelasticity"),
        ("*MOHR COULOMB\n30., 5.", "Mohr-Coulomb"),
        ("*DAMAGE INITIATION, CRITERION=DUCTILE\n0.1, 0., 0.", "damage initiation"),
        ("*DEFORMATION PLASTICITY\n2e5,0.3,4e2,10.,0.2", "deformation"),
    ],
)
def test_deferred_material_cards_reject_clearly(card, needle):
    """(6.1) A deferred Phase-2 material card fails with a clear, actionable error
    that names the capability and says it is not yet implemented — never a silent
    no-op (which would give a wrong solve) or a crash."""
    import calculixpp

    deck = (
        "*NODE\n1,0.,0.,0.\n2,1.,0.,0.\n3,0.,1.,0.\n4,0.,0.,1.\n"
        "*ELEMENT, TYPE=C3D4, ELSET=EALL\n100, 1, 2, 3, 4\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n210000., 0.3\n" + card + "\n"
        "*SOLID SECTION, ELSET=EALL, MATERIAL=M\n"
    )
    with pytest.raises(RuntimeError) as exc:
        calculixpp.summary_text(deck)
    msg = str(exc.value)
    assert needle in msg
    assert "not yet implemented" in msg


@pytest.mark.parametrize(
    "card, needle",
    [
        ("*SHELL SECTION, ELSET=EALL, MATERIAL=M\n0.01", "shell sections"),
        ("*BEAM SECTION, ELSET=EALL, MATERIAL=M, SECTION=RECT\n0.1, 0.1", "beam sections"),
        ("*MEMBRANE SECTION, ELSET=EALL, MATERIAL=M\n0.01", "membrane sections"),
    ],
)
def test_deferred_section_cards_reject_clearly(card, needle):
    """(6.1) Deferred structural-section cards fail with a clear, actionable error."""
    import calculixpp

    deck = (
        "*NODE\n1,0.,0.,0.\n2,1.,0.,0.\n3,0.,1.,0.\n4,0.,0.,1.\n"
        "*ELEMENT, TYPE=C3D4, ELSET=EALL\n100, 1, 2, 3, 4\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n210000., 0.3\n" + card + "\n"
    )
    with pytest.raises(RuntimeError) as exc:
        calculixpp.summary_text(deck)
    msg = str(exc.value)
    assert needle in msg
    assert "not yet implemented" in msg


def test_api_parity_phase2():
    """(6.2) API parity for Phase 2: every key public solve/inspect capability has
    a Python entry point. 'Parity' here means the Phase-2 solver and introspection
    surfaces are each reachable from Python:
      - solve / solve_text      : linear-static solve (auto-routes nonlinear decks)
      - solve_nonlinear[_text]  : explicit Newton-Raphson driver + options
      - summary / summary_text  : solve-free deck inspection (Phase-2 capabilities)
      - backend introspection   : available_backends / selected_backend
    plus the introspection payload each surface must expose. This test fails if a
    capability loses its binding (a regression guard on the public API surface)."""
    import calculixpp

    required_callables = [
        "solve",
        "solve_text",
        "solve_nonlinear",
        "solve_nonlinear_text",
        "summary",
        "summary_text",
        "available_backends",
        "selected_backend",
    ]
    for name in required_callables:
        assert hasattr(calculixpp, name), f"missing Python entry point: {name}"
        assert callable(getattr(calculixpp, name)), f"{name} is not callable"

    # solve_nonlinear must expose the Newton options (line_search) and the solver/
    # backend selectors — parity with the C++ NonlinearOptions surface. Exercise
    # the option via a keyword call (pybind11 builtins aren't inspect-able).
    nl = calculixpp.solve_nonlinear_text(C3D4_PRESSURE_DECK, line_search=True)
    for key in ("displacement", "converged", "newton_increments",
                "newton_iterations", "newton_cutbacks", "backend"):
        assert key in nl, f"solve_nonlinear result missing '{key}'"

    # The summary introspection surface must expose the Phase-2 capability fields.
    s = calculixpp.summary_text(C3D4_PRESSURE_DECK)
    for key in ("element_type_counts", "has_plasticity", "has_nonlinear_material",
                "num_constraints", "has_constraints", "num_springs",
                "num_amplitudes", "num_body_loads"):
        assert key in s, f"summary result missing Phase-2 field '{key}'"

    # A nonlinear (plastic) deck must be reachable through BOTH solve() (auto) and
    # solve_nonlinear() — the auto-dispatch capability has a Python entry point.
    auto = calculixpp.solve_text(_PLASTIC_CUBE)
    explicit = calculixpp.solve_nonlinear_text(_PLASTIC_CUBE)
    assert auto["converged"] is True and explicit["converged"] is True
    assert _rel_l2_disp(auto["displacement"], explicit["displacement"]) < 1e-12
