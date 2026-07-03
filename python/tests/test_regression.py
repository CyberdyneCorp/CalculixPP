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
#   beamt         1e-4      C3D20R + *EXPANSION thermal stress (obs 2.43e-6) — Phase 3
REFERENCE_CORPUS = [
    ("beam10p", 1e-4, "C3D10 quadratic-tet cantilever"),
    ("beam8p", 1e-4, "C3D8 linear-hex cantilever"),
    ("beam20p", 1e-4, "C3D20 quadratic-hex cantilever"),
    ("achtelg", 1e-5, "C3D20R cube, GRAV body load"),
    ("achtelc", 1e-5, "C3D20R cube, CENTRIF body load"),
    ("achtel2", 1e-4, "C3D20R cube, 2-term *EQUATION"),
    ("achtel9", 1e-4, "C3D20R cube, 9-term *EQUATION"),
    ("beam8p_mpc", 1e-4, "C3D8 beam with *EQUATION"),
    # Phase-3 thermal expansion (task 4.3): stock "heated beam between two fixed
    # walls" deck (*EXPANSION,ZERO=293 + *TEMPERATURE 500, C3D20R). One-way thermal
    # strain eps_th = alpha (T - Tref) reproduces the reference displacements.
    ("beamt", 1e-4, "C3D20R beam, *EXPANSION thermal stress"),
]

# Back-compat alias kept for any external reference to the old name.
REFERENCE_DECKS = [deck for deck, _tol, _desc in REFERENCE_CORPUS]

# --- Thermal reference corpus (Phase 3, tasks 3.1 / 3.2 / 3.3 / 6.6) ----------
# Heat-transfer decks with a nodal-temperature .dat.ref. All are single-C3D20R-cube
# decks (conduction is element-agnostic, sharing the mechanical shape functions/
# Gauss rules) that isolate one thermal load path each. Each entry is
# (deck, tol, mode, desc) where mode is "steady" (one .dat block, compared
# directly) or "transient" (backward-Euler march to the end-of-step steady field,
# compared against the FINAL .dat increment):
#   oneel20cf   steady     *CFLUX  concentrated nodal heat flux (obs abs err ~1.8e-15)
#   oneel20df   steady     *DFLUX  distributed surface flux S1 (hex face, ~1.3e-15)
#   oneel20fi   steady     *FILM   convective film q=h(T-T_sink) on F1 (~1e-14)
#   oneel20fi2  transient  *FILM   transient convection (backward Euler) relaxing to
#                                  the fixed-face/film steady field (obs ~2e-6)
# Steady fields are nodal and exact, so their tolerance is machine-precision-tight;
# the transient deck compares the relaxed end-of-step field (looser, ~1e-4, because
# the reference uses adaptive sub-stepping while we march fixed sub-steps).
THERMAL_CORPUS = [
    ("oneel20cf", 1e-9, "steady", "C3D20R cube, *CFLUX concentrated heat flux"),
    ("oneel20df", 1e-9, "steady", "C3D20R cube, *DFLUX surface heat flux"),
    ("oneel20fi", 1e-9, "steady", "C3D20R cube, *FILM convective boundary"),
    ("oneel20fi2", 1e-4, "transient",
     "C3D20R cube, transient *FILM convection relaxing to steady"),
]

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


def parse_ref_step_displacements(path):
    """Return a list (one per 'S T E P' block) of {node_id: (ux, uy, uz)} from a .dat.

    A multi-step deck writes one displacement block per step; this splits them so each
    step's *NODE PRINT can be compared against the driver's per-step total field."""
    steps = []
    cur = None
    indisp = False
    with open(path) as fh:
        for line in fh:
            if "S T E P" in line:
                cur = {}
                steps.append(cur)
                indisp = False
                continue
            if "displacements" in line:
                indisp = True
                continue
            if indisp:
                parts = line.split()
                if len(parts) == 4:
                    try:
                        cur[int(parts[0])] = [float(x) for x in parts[1:]]
                    except ValueError:
                        indisp = False
                elif parts:
                    indisp = False
    return steps


def test_multistep_beampfix_matches_calculix():
    """(multi-step 5.1 / 2.4) The stock two-*STEP deck ``beampfix`` reproduces stock
    CalculiX's per-step *NODE PRINT displacements. Step 1 is a *CLOAD; step 2 resets
    the load (*CLOAD, OP=NEW), applies a new *CLOAD, and holds the LOAD set at its
    step-1 (deformed) value via *BOUNDARY, FIXED. The step-loop driver carries the
    converged displacement forward, so both steps match the reference within rel-L2.
    This is the stock-CalculiX validation of multi-step analysis (the enabler for
    cross-step birth-death 5.1 and OP=MOD/NEW load accumulation 2.4)."""
    inp = os.path.join(CCX_TEST, "beampfix.inp")
    ref_path = os.path.join(CCX_TEST, "beampfix.dat.ref")
    if not (os.path.exists(inp) and os.path.exists(ref_path)):
        pytest.skip("reference deck beampfix not available")
    import calculixpp

    ref_steps = parse_ref_step_displacements(ref_path)
    res = calculixpp.solve_multistep(inp)
    assert len(res) == 2, f"beampfix has 2 *STEP blocks, got {len(res)}"
    assert len(ref_steps) == 2, "beampfix.dat.ref must have two S T E P blocks"

    for si, r in enumerate(res):
        ids = list(r["node_ids"])
        U = r["displacement"]
        rd = ref_steps[si]
        assert rd, f"step {si + 1}: empty reference displacement block"
        num = den = 0.0
        compared = 0
        for idx, nid in enumerate(ids):
            if nid not in rd:
                continue
            compared += 1
            for c in range(3):
                d = float(U[idx][c]) - rd[nid][c]
                num += d * d
                den += rd[nid][c] ** 2
        rel_l2 = math.sqrt(num / den) if den > 0 else math.sqrt(num)
        assert compared > 0, f"step {si + 1}: no reference nodes matched"
        assert rel_l2 < 1e-4, f"beampfix step {si + 1}: rel-L2 {rel_l2:.3e} >= 1e-4"


def test_multistep_single_step_equals_solve():
    """A single-*STEP deck solved via solve_multistep() returns a one-element list
    whose displacement equals solve()'s — the multi-step path never perturbs the
    single-step result (the critical gate, from Python)."""
    import calculixpp

    deck = """
*NODE, NSET=NALL
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
*MATERIAL, NAME=EL
*ELASTIC
210000., 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=EL
*NSET, NSET=FIX
1, 4, 5, 8
*NSET, NSET=LOAD
2, 3, 6, 7
*BOUNDARY
FIX, 1, 3, 0.
*STEP
*STATIC
*CLOAD
LOAD, 1, 100.
*END STEP
"""
    single = calculixpp.solve_text(deck)
    multi = calculixpp.solve_multistep_text(deck)
    assert len(multi) == 1
    Us = single["displacement"]
    Um = multi[0]["displacement"]
    for i in range(len(Us)):
        for c in range(3):
            assert abs(float(Us[i][c]) - float(Um[i][c])) < 1e-12


def test_multistep_two_load_steps_superpose():
    """Two linear load steps of F each (OP=MOD accumulation) sum to the single-step 2F
    result — the core multi-step superposition property, from Python."""
    import calculixpp

    header = """
*NODE, NSET=NALL
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
*MATERIAL, NAME=EL
*ELASTIC
210000., 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=EL
*NSET, NSET=FIX
1, 4, 5, 8
*NSET, NSET=LOAD
2, 3, 6, 7
*BOUNDARY
FIX, 1, 3, 0.
"""
    two = header + (
        "*STEP\n*STATIC\n*CLOAD\nLOAD, 1, 50.\n*END STEP\n"
        "*STEP\n*STATIC\n*CLOAD\nLOAD, 1, 50.\n*END STEP\n"
    )
    one = header + "*STEP\n*STATIC\n*CLOAD\nLOAD, 1, 100.\n*END STEP\n"
    res2 = calculixpp.solve_multistep_text(two)
    res1 = calculixpp.solve_text(one)
    assert len(res2) == 2
    U2 = res2[-1]["displacement"]  # end of step 2 == total 2F state
    U1 = res1["displacement"]
    for i in range(len(U1)):
        for c in range(3):
            assert abs(float(U2[i][c]) - float(U1[i][c])) < 1e-9 * (
                1.0 + abs(float(U1[i][c]))
            )


def test_multistep_element_birth_death_python():
    """(5.1) Cross-step element birth-death: two stacked confined cubes, cube 2
    removed in step 1 and re-added in step 2, is strain-free relative to the deformed
    config — its step-2 stress is built only from the step-2 increment, so it is
    non-zero yet smaller than the never-removed reference's step-2 total, and exactly
    zero if it stays removed. Validates the *MODEL CHANGE active mask carried per step."""
    import calculixpp

    def deck(s1_mc, s2_mc, d1, d2):
        header = """
*NODE, NSET=NALL
1, 0., 0., 0.
2, 1., 0., 0.
3, 1., 1., 0.
4, 0., 1., 0.
5, 0., 0., 1.
6, 1., 0., 1.
7, 1., 1., 1.
8, 0., 1., 1.
9,  0., 2., 0.
10, 1., 2., 0.
11, 0., 2., 1.
12, 1., 2., 1.
*ELEMENT, TYPE=C3D8, ELSET=CUBE1
1, 1, 2, 3, 4, 5, 6, 7, 8
*ELEMENT, TYPE=C3D8, ELSET=CUBE2
2, 4, 3, 10, 9, 8, 7, 12, 11
*ELSET, ELSET=EALL
1, 2
*MATERIAL, NAME=EL
*ELASTIC
210000., 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=EL
*NSET, NSET=TOP
9, 10, 11, 12
*NSET, NSET=MID
3, 4, 7, 8
*NSET, NSET=BOT
1, 2, 5, 6
*NSET, NSET=ALLN
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
*BOUNDARY
ALLN, 1, 1, 0.
ALLN, 3, 3, 0.
BOT, 2, 2, 0.
MID, 2, 2, 0.
"""
        return header + (
            f"*STEP\n*STATIC\n{s1_mc}*BOUNDARY\nTOP, 2, 2, {d1}\n*END STEP\n"
            f"*STEP\n*STATIC\n{s2_mc}*BOUNDARY\nTOP, 2, 2, {d2}\n*END STEP\n"
        )

    bd = calculixpp.solve_multistep_text(
        deck("*MODEL CHANGE, TYPE=ELEMENT, REMOVE\n2\n",
             "*MODEL CHANGE, TYPE=ELEMENT, ADD\n2\n", 0.01, 0.02))
    ref = calculixpp.solve_multistep_text(deck("", "", 0.01, 0.02))
    dead = calculixpp.solve_multistep_text(
        deck("*MODEL CHANGE, TYPE=ELEMENT, REMOVE\n2\n",
             "*MODEL CHANGE, TYPE=ELEMENT, REMOVE\n2\n", 0.01, 0.02))

    ids = list(bd[0]["node_ids"])
    i10 = ids.index(10)  # cube-2-exclusive node
    bd_syy_1 = float(bd[0]["stress"][i10][1])
    bd_syy_2 = float(bd[1]["stress"][i10][1])
    ref_syy_2 = float(ref[1]["stress"][i10][1])
    dead_syy_2 = float(dead[1]["stress"][i10][1])

    assert abs(bd_syy_1) < 1e-6            # dead in step 1 -> no stress
    assert abs(bd_syy_2) > 1e-3            # live in step 2 -> stressed
    assert abs(bd_syy_2) < abs(ref_syy_2)  # but less than the never-removed total
    assert abs(dead_syy_2) < 1e-9         # removed throughout -> exactly zero


def parse_ref_temperatures(path, last=False):
    """Return {node_id: temperature} from a 'temperatures' block of a .dat.

    A steady deck writes one block (last=False, the default). A transient deck
    writes one block per increment; pass last=True to read the final increment
    (the relaxed steady field at the end of the step period)."""
    temps = {}
    with open(path) as fh:
        lines = fh.readlines()
    headers = [k for k, l in enumerate(lines) if "temperatures" in l]
    if not headers:
        return temps
    i = headers[-1] if last else headers[0]
    for line in lines[i + 1:]:
        parts = line.split()
        if len(parts) == 2:
            try:
                temps[int(parts[0])] = float(parts[1])
            except ValueError:
                if temps:
                    break
        elif parts and temps:
            break
    return temps


@pytest.mark.parametrize(
    "deck,tol,mode,desc",
    THERMAL_CORPUS,
    ids=[d for d, _t, _m, _s in THERMAL_CORPUS],
)
def test_thermal_reference_corpus_matches_calculix(deck, tol, mode, desc):
    """(3.1 / 3.2 / 3.3 / 6.6) Each heat-transfer deck reproduces stock-CalculiX
    .dat.ref nodal temperatures within tolerance (manifest: THERMAL_CORPUS).

    A steady deck is compared against the single .dat block; a transient deck
    marches backward-Euler over the step period and is compared against the FINAL
    increment (the relaxed steady field)."""
    import calculixpp as cx

    inp = os.path.join(CCX_TEST, deck + ".inp")
    ref_path = os.path.join(CCX_TEST, deck + ".dat.ref")
    if not (os.path.exists(inp) and os.path.exists(ref_path)):
        pytest.skip(f"reference deck {deck} not available")
    transient = mode == "transient"
    ref = parse_ref_temperatures(ref_path, last=transient)
    assert ref, f"{deck}: no reference temperatures parsed"

    res = cx.solve(inp)
    expected_proc = ("heat transfer transient" if transient
                     else "heat transfer steady state")
    assert res.get("procedure") == expected_proc, (
        f"{deck}: procedure {res.get('procedure')!r} != {expected_proc!r}")
    ids = res["node_ids"]
    temp = res["temperature"]
    id2t = {int(ids[i]): float(temp[i]) for i in range(len(ids))}

    compared = 0
    max_abs = 0.0
    for nid, tref in ref.items():
        if nid in id2t:
            compared += 1
            max_abs = max(max_abs, abs(id2t[nid] - tref))
    assert compared > 0, f"{deck}: no reference nodes matched the model"
    assert max_abs < tol, f"{deck} ({desc}): max abs temp err {max_abs:.3e} >= {tol:.0e}"


def parse_ref_heat_flux(path):
    """Return {(elem_id, gp): (qx, qy, qz)} from a *EL PRINT HFL 'heat flux' block."""
    hfl = {}
    with open(path) as fh:
        lines = fh.readlines()
    headers = [k for k, l in enumerate(lines) if "heat flux" in l]
    if not headers:
        return hfl
    for line in lines[headers[0] + 1:]:
        parts = line.split()
        if len(parts) == 5:
            try:
                hfl[(int(parts[0]), int(parts[1]))] = [float(x) for x in parts[2:]]
            except ValueError:
                if hfl:
                    break
        elif parts and hfl:
            break
    return hfl


def test_heat_flux_hfl_matches_calculix():
    """(6.2) Integration-point heat flux HFL reproduces stock-CalculiX oneel20cf.dat.ref
    *EL PRINT HFL to machine precision. The C3D20R cube driven by a *CFLUX carries a
    uniform qx = 120 (heat generation / area) and ~1e-14 transverse flux at every one
    of its 8 integration points; we match the reference block per (elem, integ.pnt.)."""
    import calculixpp as cx

    inp = os.path.join(CCX_TEST, "oneel20cf.inp")
    ref_path = os.path.join(CCX_TEST, "oneel20cf.dat.ref")
    if not (os.path.exists(inp) and os.path.exists(ref_path)):
        pytest.skip("reference deck oneel20cf not available")
    ref = parse_ref_heat_flux(ref_path)
    assert ref, "no reference HFL parsed"

    res = cx.solve(inp)
    elem = res["hfl_elem"]
    gp = res["hfl_gp"]
    flux = res["hfl_flux"]
    compared = 0
    max_abs = 0.0
    for k in range(len(elem)):
        key = (int(elem[k]), int(gp[k]))
        if key in ref:
            compared += 1
            for c in range(3):
                max_abs = max(max_abs, abs(float(flux[k][c]) - ref[key][c]))
    assert compared == len(ref), f"matched {compared}/{len(ref)} HFL points"
    assert max_abs < 1e-9, f"HFL max abs err {max_abs:.3e} >= 1e-9"


_COUPLED_BAR = """
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
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*CONDUCTIVITY
50.
*EXPANSION
1.2e-5
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*STEP
*COUPLED TEMPERATURE-DISPLACEMENT
*BOUNDARY
1, 11, 11, 100.
4, 11, 11, 100.
5, 11, 11, 100.
8, 11, 11, 100.
2, 11, 11, 100.
3, 11, 11, 100.
6, 11, 11, 100.
7, 11, 11, 100.
*BOUNDARY
1, 1
4, 1
5, 1
8, 1
2, 1
3, 1
6, 1
7, 1
1, 2
1, 3
2, 3
*END STEP
"""


def test_coupled_temperature_displacement_thermal_stress():
    """(4.1 / 4.3) *COUPLED TEMPERATURE-DISPLACEMENT one-way coupling through the
    Python binding: a bar heated uniformly by dT=100 with u_x fixed on both end faces
    develops the analytical uniaxial thermal stress sigma_xx = -E*alpha*dT = -252, and
    the result dict carries both the temperature and mechanical fields."""
    import calculixpp as cx

    res = cx.solve_text(_COUPLED_BAR)
    assert res.get("procedure") == "coupled temperature-displacement"
    assert "temperature" in res and "stress" in res and "reaction" in res
    # Uniform T=100 field (Tref=0) reproduced on every node.
    for t in res["temperature"]:
        assert abs(float(t) - 100.0) < 1e-6
    # sigma_xx = -E*alpha*dT on every node; lateral components relax to zero.
    sxx_expect = -210000.0 * 1.2e-5 * 100.0
    for row in res["stress"]:
        assert abs(float(row[0]) - sxx_expect) < 1e-2
        assert abs(float(row[1])) < 1e-2
        assert abs(float(row[2])) < 1e-2


def test_coupled_monolithic_equals_staggered():
    """(4.1 / 4.2) The MONOLITHIC 4-DOF/node coupled solve equals the STAGGERED
    (one-way) solve for a pure thermal-stress problem — no mechanical->thermal feedback
    means the coupled tangent is block-triangular, so the single 4-DOF solve reproduces
    the sequential result. Both schemes are reachable from a deck via SOLUTIONS=."""
    import calculixpp as cx

    stag = cx.solve_text(
        _COUPLED_BAR.replace(
            "*COUPLED TEMPERATURE-DISPLACEMENT",
            "*COUPLED TEMPERATURE-DISPLACEMENT, SOLUTIONS=STAGGERED",
        )
    )
    mono = cx.solve_text(
        _COUPLED_BAR.replace(
            "*COUPLED TEMPERATURE-DISPLACEMENT",
            "*COUPLED TEMPERATURE-DISPLACEMENT, SOLUTIONS=MONOLITHIC",
        )
    )
    for ts, tm in zip(stag["temperature"], mono["temperature"]):
        assert abs(float(ts) - float(tm)) < 1e-9
    for rs, rm in zip(stag["stress"], mono["stress"]):
        for cs, cm in zip(rs, rm):
            assert abs(float(cs) - float(cm)) < 1e-6


_KT_BAR = """
*NODE
1, 0., 0., 0.
2, 1., 0., 0.
3, 1., 1., 0.
4, 0., 1., 0.
5, 0., 0., 1.
6, 1., 0., 1.
7, 1., 1., 1.
8, 0., 1., 1.
9, 2., 0., 0.
10, 2., 1., 0.
11, 2., 0., 1.
12, 2., 1., 1.
*ELEMENT, TYPE=C3D8, ELSET=EALL
1, 1, 2, 3, 4, 5, 6, 7, 8
2, 2, 9, 10, 3, 6, 11, 12, 7
*MATERIAL, NAME=EL
*CONDUCTIVITY
50., 0.
100., 100.
*SOLID SECTION, ELSET=EALL, MATERIAL=EL
*STEP
*HEAT TRANSFER, STEADY STATE
*BOUNDARY
1, 11, 11, 0.
4, 11, 11, 0.
5, 11, 11, 0.
8, 11, 11, 0.
9, 11, 11, 100.
10, 11, 11, 100.
11, 11, 11, 100.
12, 11, 11, 100.
*END STEP
"""


def test_temperature_dependent_conductivity_python():
    """(6.1) Temperature-dependent *CONDUCTIVITY k(T)=50+0.5T through the Python solve:
    a 1-D bar (0..2) with fixed end temperatures 0 and 100 has a CONSTANT flux, so the
    Kirchhoff transform q L = int_0^100 (50+0.5T) dT = 7500 fixes the nonlinear profile.
    The mid-plane node (x=1) sits at Tm = -100 + sqrt(25000) ~ 58.11 (NOT 50 as a
    constant k would give), and the hot-face heat-flux reaction sums to q A = 3750."""
    import calculixpp as cx
    import math

    res = cx.solve_text(_KT_BAR)
    assert res.get("procedure") == "heat transfer steady state"
    ids = list(res["node_ids"])
    temp = res["temperature"]
    rfl = res["flux_reaction"]
    coords = res["node_coords"]
    id2i = {int(ids[i]): i for i in range(len(ids))}

    Tm = -100.0 + math.sqrt(25000.0)  # ~58.113883
    assert abs(float(temp[id2i[2]]) - Tm) < 1e-2, float(temp[id2i[2]])
    assert float(temp[id2i[2]]) > 55.0  # meaningfully above the constant-k value 50

    hot = sum(float(rfl[i]) for i in range(len(ids)) if abs(coords[i][0] - 2.0) < 1e-9)
    assert abs(hot - 3750.0) < 1e-1, hot


_CAVITY_PLATES = """
*NODE
1, -0.01, 0., 0.
2, 0., 0., 0.
3, 0., 1., 0.
4, -0.01, 1., 0.
5, -0.01, 0., 1.
6, 0., 0., 1.
7, 0., 1., 1.
8, -0.01, 1., 1.
9, 0.5, 0., 0.
10, 0.51, 0., 0.
11, 0.51, 1., 0.
12, 0.5, 1., 0.
13, 0.5, 0., 1.
14, 0.51, 0., 1.
15, 0.51, 1., 1.
16, 0.5, 1., 1.
*ELEMENT, TYPE=C3D8, ELSET=EALL
1, 1, 2, 3, 4, 5, 6, 7, 8
2, 9, 10, 11, 12, 13, 14, 15, 16
*MATERIAL, NAME=EL
*CONDUCTIVITY
1e6
*SOLID SECTION, ELSET=EALL, MATERIAL=EL
*PHYSICAL CONSTANTS, ABSOLUTE ZERO=0., STEFAN BOLTZMANN=5.67E-8
*STEP
*HEAT TRANSFER, STEADY STATE
*BOUNDARY
1, 11, 11, 800.
4, 11, 11, 800.
5, 11, 11, 800.
8, 11, 11, 800.
10, 11, 11, 500.
11, 11, 11, 500.
14, 11, 11, 500.
15, 11, 11, 500.
*RADIATE
1, R4CR, , 0.8
2, R6CR, , 0.8
*END STEP
"""


def test_cavity_radiation_parallel_plates_python():
    """(3.4) Gray-body cavity radiation (*RADIATE ...,CR) end to end through the Python
    binding. Two thin, high-conductivity unit-square plates face each other across a
    gap; the plate outer faces are held at 800 K and 500 K, so each plate is essentially
    isothermal and the cavity faces sit at the boundary temperatures. The two facing
    patches form a 2-body enclosure (F12 = F21 = 1), so in steady state the heat that
    conducts through the hot plate (its flux reaction) equals the classic parallel-plate
    radiative exchange q = sigma (T1^4 - T2^4)/(1/e1 + 1/e2 - 1) over the unit area."""
    import calculixpp as cx

    res = cx.solve_text(_CAVITY_PLATES)
    assert res.get("procedure") == "heat transfer steady state"
    ids = list(res["node_ids"])
    coords = res["node_coords"]
    rfl = res["flux_reaction"]

    sig, Th, Tc, e = 5.67e-8, 800.0, 500.0, 0.8
    q_expect = sig * (Th ** 4 - Tc ** 4) / (1.0 / e + 1.0 / e - 1.0)  # ~13120.4 W

    # Sum the flux reaction over the hot-plate outer face (x = -0.01): in steady state it
    # carries the full radiative heat flow Q across the cavity.
    hot = sum(
        float(rfl[i]) for i in range(len(ids)) if abs(coords[i][0] + 0.01) < 1e-9
    )
    assert abs(hot - q_expect) < 1e-1, (hot, q_expect)
    # Energy is conserved: the cold-plate outer face absorbs the same flow (opposite sign).
    cold = sum(
        float(rfl[i]) for i in range(len(ids)) if abs(coords[i][0] - 0.51) < 1e-9
    )
    assert abs(cold + q_expect) < 1e-1, (cold, q_expect)


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


# --- Phase-3 thermal API parity (task 6.4) ------------------------------------
# A minimal single-C3D8-cube steady conduction deck: face x=0 held at 0, face x=1
# held at 100 (temperature BCs on DOF 11). The auto-dispatch routes *HEAT TRANSFER
# to the scalar thermal solver; the result carries temperature/flux_reaction.
_HEAT_CUBE = """
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
*MATERIAL, NAME=COND
*CONDUCTIVITY
50.
*SOLID SECTION, ELSET=EALL, MATERIAL=COND
*STEP
*HEAT TRANSFER, STEADY STATE
*BOUNDARY
1, 11, 11, 0.
4, 11, 11, 0.
5, 11, 11, 0.
8, 11, 11, 0.
2, 11, 11, 100.
3, 11, 11, 100.
6, 11, 11, 100.
7, 11, 11, 100.
*END STEP
"""


def test_api_parity_phase3_thermal():
    """(6.4) API parity for the Phase-3 THERMAL track: the heat-transfer solve and
    its nodal-temperature output are reachable from Python, the summary surface
    reports the thermal procedure without solving, and a Phase-3 CONTACT card is
    rejected with an actionable error (never silently mis-solved). Parity means the
    thermal capabilities exposed in C++ each have a Python entry point."""
    import calculixpp as cx

    # solve()/solve_text() auto-dispatch a *HEAT TRANSFER deck to the thermal solver
    # and return the scalar temperature field (NT) and heat-flux reaction (RFL).
    res = cx.solve_text(_HEAT_CUBE)
    assert res.get("procedure") == "heat transfer steady state"
    for key in ("temperature", "flux_reaction", "node_ids", "node_coords",
                "num_nodes", "num_elements", "backend"):
        assert key in res, f"thermal result missing '{key}'"
    ids = res["node_ids"]
    id2t = {int(ids[i]): float(res["temperature"][i]) for i in range(len(ids))}
    # x=0 face held at 0, x=1 face held at 100; interior obeys the linear profile.
    for nid in (1, 4, 5, 8):
        assert abs(id2t[nid] - 0.0) < 1e-9
    for nid in (2, 3, 6, 7):
        assert abs(id2t[nid] - 100.0) < 1e-9
    # The flux reaction at the driven faces balances (sum ~ 0 for a closed system).
    assert abs(sum(float(x) for x in res["flux_reaction"])) < 1e-6

    # summary()/summary_text() report the thermal procedure without solving — the
    # introspection surface distinguishes heat transfer from a mechanical static.
    s = cx.summary_text(_HEAT_CUBE)
    assert s["procedure"] == "heat transfer steady state"
    assert cx.summary_text(C3D4_PRESSURE_DECK)["procedure"] == "static"

    # A Phase-3 CONTACT card must raise a clear, actionable error naming the contact
    # workstream (deferred to the next workflow) — never a silent wrong solve.
    contact_deck = _HEAT_CUBE.replace(
        "*STEP\n", "*SURFACE INTERACTION, NAME=SI\n*SURFACE BEHAVIOR\n1e3\n*STEP\n")
    with pytest.raises(RuntimeError) as ei:
        cx.solve_text(contact_deck)
    msg = str(ei.value).lower()
    assert "contact" in msg and "not implemented" in msg


# --- *MODEL CHANGE element birth-death (tasks 5.1), reachable from Python ------
# Two confined unit C3D8 cubes stacked in y (share the y=1 face). Every node has
# u_y=u_z=0 (so element removal never orphans a DOF) and u_x=0 on x=0 / u_x=delta
# on x=1. Each active cube contributes -M*delta to the fixed-face x-reaction
# (M = confined modulus), so removing a cube halves the reaction and re-adding it
# restores the full model. `{MC}` is spliced into the single step.
_MODEL_CHANGE_BAR = """
*NODE, NSET=NALL
1, 0., 0., 0.
2, 1., 0., 0.
3, 1., 1., 0.
4, 0., 1., 0.
5, 0., 0., 1.
6, 1., 0., 1.
7, 1., 1., 1.
8, 0., 1., 1.
9,  0., 2., 0.
10, 1., 2., 0.
11, 0., 2., 1.
12, 1., 2., 1.
*ELEMENT, TYPE=C3D8, ELSET=CUBE1
1, 1, 2, 3, 4, 5, 6, 7, 8
*ELEMENT, TYPE=C3D8, ELSET=CUBE2
2, 4, 3, 10, 9, 8, 7, 12, 11
*ELSET, ELSET=EALL
1, 2
*MATERIAL, NAME=EL
*ELASTIC
210000., 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=EL
*BOUNDARY
1, 2, 3
2, 2, 3
3, 2, 3
4, 2, 3
5, 2, 3
6, 2, 3
7, 2, 3
8, 2, 3
9, 2, 3
10, 2, 3
11, 2, 3
12, 2, 3
1, 1, 1
4, 1, 1
5, 1, 1
8, 1, 1
9, 1, 1
11, 1, 1
2, 1, 1, 0.01
3, 1, 1, 0.01
6, 1, 1, 0.01
7, 1, 1, 0.01
10, 1, 1, 0.01
12, 1, 1, 0.01
*STEP
*STATIC
{MC}*END STEP
"""


def _fixed_face_rx(res):
    """Sum the x-reaction over the x=0 face nodes (ids 1,4,5,8,9,11)."""
    node_ids = list(res["node_ids"])
    rf = res["reaction"]
    total = 0.0
    for nid in (1, 4, 5, 8, 9, 11):
        k = node_ids.index(nid)
        total += float(rf[k][0])
    return total


def test_model_change_element_birth_death_python():
    """ANALYTICAL (model change): a *MODEL CHANGE, TYPE=ELEMENT deck solves through
    the Python solve_text() entry point, and removing an element halves the
    fixed-face reaction while re-adding it (REMOVE then ADD) reproduces the full
    model exactly (strain-free reactivation). Validates tasks 5.1 end to end from
    Python (spec: model-change — removed elements carry no load / strain-free
    reactivation, invocable from the Python bindings)."""
    import calculixpp as cx

    full = cx.solve_text(_MODEL_CHANGE_BAR.format(MC=""))
    removed = cx.solve_text(
        _MODEL_CHANGE_BAR.format(MC="*MODEL CHANGE, TYPE=ELEMENT, REMOVE\n2\n"))
    readded = cx.solve_text(
        _MODEL_CHANGE_BAR.format(
            MC="*MODEL CHANGE, TYPE=ELEMENT, REMOVE\n2\n"
               "*MODEL CHANGE, TYPE=ELEMENT, ADD\n2\n"))

    rx_full = _fixed_face_rx(full)
    rx_removed = _fixed_face_rx(removed)
    rx_readded = _fixed_face_rx(readded)

    # Removing one of two identical cubes halves the reaction.
    assert abs(rx_removed - 0.5 * rx_full) < 1e-4 * abs(rx_full), (
        f"removed {rx_removed} vs half-full {0.5 * rx_full}")
    # REMOVE then ADD reproduces the full model to solver precision (strain-free).
    assert abs(rx_readded - rx_full) < 1e-7 * abs(rx_full), (
        f"re-added {rx_readded} vs full {rx_full}")
    assert _rel_l2_disp(readded["displacement"], full["displacement"]) < 1e-10


def test_model_change_contact_pair_parses_python():
    """A *MODEL CHANGE, TYPE=CONTACT PAIR deck parses and solves (the contact-pair
    record is stored for the contact workflow; it does not alter the mechanical
    solve here). Confirms the card is accepted through the Python path (tasks 5.2,
    parse-only)."""
    import calculixpp as cx

    with_pair = cx.solve_text(
        _MODEL_CHANGE_BAR.format(
            MC="*MODEL CHANGE, TYPE=CONTACT PAIR, REMOVE\nSURF_A, SURF_B\n"))
    baseline = cx.solve_text(_MODEL_CHANGE_BAR.format(MC=""))
    # Contact-pair model change does not touch the element assembly here.
    assert _rel_l2_disp(with_pair["displacement"], baseline["displacement"]) < 1e-12
