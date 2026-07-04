# Add a damped complex modal eigensolver (`*COMPLEX FREQUENCY`, proportional damping)

## Why

`*COMPLEX FREQUENCY` is currently deferred (BACKLOG item **1.6 / 2.3**, "Complex /
damped modes"). Its two current spec requirements — one in `eigensolution` ("Complex and
damped modes"), one in `modal-and-buckling-analysis` ("Complex frequency analysis") — are
placeholders that gesture vaguely at "damping, Coriolis, or friction effects" with a
"NumPP complex eigensolver". Nothing is implemented; the card fails with a deferral parse
error.

We can ship the **damped complex modes** slice now without any upstream blocker. The
kernel it needs is a **dense** complex non-symmetric eigensolve on a small
`2·nev × 2·nev` matrix (`nev` = number of real modes, typically 5–50), which already
exists as `numpp::linalg::eig()` (verified: `NumPP linalg.hpp`, returns a complex
`EigResult`). The missing SciPP **sparse** complex / non-symmetric eigensolver is **not**
required for this path, because we reduce onto the already-shipped real mass-normalized
`*FREQUENCY` eigenbasis rather than doing a full-scale complex eigensolve.

### Honest scope: this is option (B), NOT CalculiX `CORIOLIS`

This change ships **ABAQUS-style subspace-projected complex modes for proportional
(Rayleigh / modal) damping** — the reduced quadratic eigenproblem
`(λ²M + λC + K)x = 0` projected onto the real modal basis. It is deliberately **not** a
port of CalculiX's `*COMPLEX FREQUENCY, CORIOLIS`, and this proposal does not claim
CalculiX numerical fidelity for it. Verified against the reference source
(`/home/leonardo/work/CalculiX/src`):

- `op_corio.f` header: *"y=A*x for real sparse **antisymmetric** matrices, i.e. the
  transpose is the negative matrix: A^T=-A"* — the CalculiX `CORIOLIS` operator is a
  **skew (gyroscopic)** matrix, **not** a damping matrix.
- `complexfreq.c:607-610`: the reduced matrix `cc` is explicitly **anti-symmetrized**
  (`cc[i*nev+j] = -cc[j*nev+i]`).
- `coriolissolve.f:90-94`: it linearizes `(-ω²I + Λ + i·ω·G_r)q = 0` with an explicit
  imaginary unit `(0,1)` on the ω-linear term — a rotor whirl / Campbell-split problem,
  producing gyroscopically-split real frequency squares, **not** damping ratios.
- `complexfrequencys.f:86-92`: CalculiX **errors** unless `CORIOLIS` or `FLUTTER` is
  present; `CORIOLIS` then consumes a **rotation body load** (rotor speed `om` + axis
  `p2` from `mafillcorio`/`xbody`), **not** a `*DAMPING` card.

CalculixPP has no Coriolis element operator and no rotor-speed/axis parser path, so a
faithful `CORIOLIS` port is out of scope. To avoid the trap of accepting a `CORIOLIS`
deck we would then mis-solve with the wrong physics, the parser **rejects** `CORIOLIS`
and `FLUTTER` with an explicit "not yet implemented" error. The keyword-less
`*COMPLEX FREQUENCY` (or an explicit "proportional damping" acceptance) drives the
option-(B) path. The true gyroscopic `CORIOLIS` operator, `FLUTTER` complex applied
force, and cyclic-symmetry complex eigenproblem remain documented follow-ons.

## What Changes

- **New solver path** — `extract_complex_modes(real_basis, damping, num_modes)`: reduce
  the proportional damping operator onto the real mass-normalized eigenbasis `Φ`
  (`Φᵀ M Φ = I`, `Φᵀ K Φ = Λ = diag(ω_k²)`), form the reduced quadratic modal problem
  `(λ²I + λC_r + Λ)q = 0`, linearize to a real `2·nev × 2·nev` companion pencil, solve
  with `numpp::linalg::eig()`, and extract per-mode complex eigenvalue, damped natural
  frequency, undamped frequency, damping ratio, decay rate, and complex mode shape.
- **Reduced-operator interface built to hold a skew, imaginary coupling** — so a future
  true gyroscopic `G_r` (skew, with `i·ω` coupling) plugs in without redesign, rather
  than assuming the reduced operator is symmetric/diagonal.
- **Analytical validation only** — Rayleigh damping gives a **diagonal** `C_r`
  (`C_r,kk = α + β·ω_k²`), so each mode has the closed form
  `λ_k = -ζ_k ω_k ± i ω_k √(1-ζ_k²)` with `ζ_k = (α/ω_k + β·ω_k)/2`. This is an exact
  correctness proof for the shipped physics — it is **not** a CalculiX-fidelity check.
- **Dense oracle** — a small-problem full `2n` state-space solve
  `[[0,I],[-M⁻¹K,-M⁻¹C]]` (guarded by a max-DOF cap) cross-checks that the modal
  reduction and the linearization agree.
- **Parser** — `*COMPLEX FREQUENCY` sets `Procedure::ComplexFrequency`, parses the mode
  count, requires a preceding `*FREQUENCY` step, and **rejects** `CORIOLIS` / `FLUTTER`.
- **CLI + Python bindings** — a damped-frequency / damping-ratio / stability table and a
  `complex_frequency` result dict (eigenvalues, damped frequencies, damping ratios,
  complex mode shapes).

## Capabilities

### New Capabilities

None. This change modifies existing capabilities.

### Modified Capabilities

- **eigensolution** — replace the vague "Complex and damped modes" requirement with the
  concrete proportional-damping modal-reduction contract, and tighten "Shared basis
  across procedures" so the complex-frequency consumer reuses the preceding `*FREQUENCY`
  basis instead of re-extracting or doing a full-scale complex eigensolve.
- **modal-and-buckling-analysis** — rewrite "Complex frequency analysis" to the shipped
  proportional-damping physics (dropping the false "damping, Coriolis, or friction" and
  "NumPP complex eigensolver at full scale" wording), keep the missing-`*FREQUENCY` error
  requirement, add the `CORIOLIS`/`FLUTTER` rejection requirement, and **replace** the
  CalculiX reference-fidelity requirement for this procedure with an **analytical**
  closed-form fidelity requirement (CalculiX `*COMPLEX FREQUENCY` output is a different
  eigenproblem and would not match).

## Impact

- Code: `include/calculixpp/numerics/eigensolution.hpp`,
  `src/numerics/eigensolution.cpp`, `include/calculixpp/core/model.hpp`,
  `src/io/inp_parser.cpp`, `apps/cli/main.cpp`, `python/bindings.cpp`, plus parser +
  numeric regression tests.
- Behavior: `*COMPLEX FREQUENCY` (no keyword / proportional damping) now solves instead
  of deferring; `CORIOLIS` / `FLUTTER` change from a generic deferral to an explicit
  "not yet implemented" rejection.
- No breaking changes to `*FREQUENCY`, `*MODAL DYNAMIC`, `*STEADY STATE DYNAMICS`, or
  `*DYNAMIC`. GPU backends remain optional.
- BACKLOG **1.6 / 2.3** is **partially** resolved (proportional-damping complex modes);
  CalculiX `CORIOLIS`/`FLUTTER` and cyclic-symmetry (6.x) remain open follow-ons and are
  re-scoped, not marked done.
