# Design — damped complex modal eigensolver (`*COMPLEX FREQUENCY`, proportional damping)

## Context and honest framing

CalculixPP already ships the real symmetric `*FREQUENCY` engine: `extract_modes(K, M, n)`
returns an `EigenBasis` whose columns are mass-normalized (`Φᵀ M Φ = I`,
`Φᵀ K Φ = Λ = diag(ω_k²)`). This change adds the **damped complex modes** step on top of
that basis.

This is **option (B)** in the review: an ABAQUS-style subspace-projected complex-modes
capability for **proportional (Rayleigh / modal) damping**. It is **not** a port of
CalculiX's `*COMPLEX FREQUENCY, CORIOLIS`. That distinction is load-bearing for every
downstream decision (parser contract, validation strategy, spec requirements) and is
grounded in the reference source, not assumed:

| Aspect | CalculiX `CORIOLIS` (verified) | This change (option B) |
|---|---|---|
| Reduced operator | **Skew** `G_r` (`op_corio.f`: `A^T=-A`; `complexfreq.c:607-610` anti-symmetrizes `cc`) | **Symmetric** `C_r` (diagonal for Rayleigh) |
| Linearized problem | `(-ω²I + Λ + i·ω·G_r)q=0` (`coriolissolve.f:90-94`, explicit `(0,1)` on ω-term) | `(λ²I + λC_r + Λ)q=0`, real companion pencil |
| Physics | Rotor whirl / Campbell forward-backward split | Damped decay + damped oscillation |
| Required input | Rotation body load (rotor speed `om`, axis `p2` via `xbody`) | A `*DAMPING` card (Rayleigh / modal ratios) |
| Output | `eig = √(gyroscopic complex freq²)` | ζ (damping ratio), f_d (damped frequency) |

Because the two are genuinely different eigenproblems, this change **does not** claim
CalculiX numerical fidelity and **does not** carry a "match reference CalculiX" spec
requirement for `*COMPLEX FREQUENCY`. Its correctness proof is an **analytical** closed
form (exact for proportional damping). To avoid silently mis-solving, the parser
**rejects** `CORIOLIS` and `FLUTTER`.

## Math

### Reduced quadratic modal eigenproblem

Physical quadratic eigenproblem: `(λ²M + λC + K)x = 0`. Substitute `x = Φq`, premultiply
`Φᵀ` and use `Φᵀ M Φ = I`, `Φᵀ K Φ = Λ`:

```
(λ² I + λ C_r + Λ) q = 0,   C_r = Φᵀ C Φ (nev×nev),   Λ = diag(ω_k²)
```

Reduced damping for the shipped proportional models:

- **Rayleigh** `C = αM + βK` ⇒ `C_r = αI + βΛ` — **diagonal**, `C_r,kk = α + β·ω_k²`.
  Each mode is an uncoupled SDOF with
  `λ_k = -ζ_k ω_k ± i ω_k √(1-ζ_k²)`, `ζ_k = (α/ω_k + β·ω_k)/2`. This equals
  `Damping::ratio(k, ω_k)`; reuse it so the complex path is consistent with the
  modal-dynamics path.
- **`*MODAL DAMPING` ratios** `ζ_k` ⇒ `C_r = diag(2 ζ_k ω_k)` — also diagonal.

Proportional damping keeps `C_r` diagonal, so the modes stay uncoupled and the closed
form above is exact. Non-proportional damping (future contact friction, or the future
skew Coriolis `G_r`) makes the reduced operator full/coupled and requires the general
eigensolve; accuracy then depends on `nev` and is a documented follow-on.

### Linearization to a real companion pencil

Linearize the reduced quadratic to a first-order `2·nev` problem. In mass-normalized
reduced coordinates `M_r = I`, the companion form is a **standard** eigenproblem
`A z = μ z` with `μ = λ`, `z = [q; λq]`:

```
A = [[   0,    I  ],
     [ -Λ,   -C_r ]]        (2·nev × 2·nev, real)
```

Feed `A` (real, non-symmetric) directly to `numpp::linalg::eig(A)`, which returns complex
eigenvalues + eigenvectors.

### Reduced-operator interface (future-proofed for skew Coriolis)

Even though `C_r` is symmetric/diagonal for option (B), the assembly interface is defined
to carry (a) a symmetric damping block **and** (b) an optional **skew** block with an
**imaginary** ω-coupling, so the true gyroscopic `G_r` (skew, `i·ω·G_r`) can populate the
second block later and reuse the same linearization scaffolding. The option-(B) path
leaves the skew/imaginary block empty. This is explicitly to make the "future skew G_r
plugs in without redesign" claim true rather than aspirational.

### Post-processing

For each complex eigenvalue `λ`:

- damped angular frequency `ω_d = |Im(λ)|`; damped frequency `f_d = ω_d / (2π)`
- undamped angular frequency `ω_n = |λ|`
- damping ratio `ζ = -Re(λ) / |λ|` (ζ>0 stable/decaying, ζ<0 unstable/growing)
- decay rate = `Re(λ)`

Keep one representative per complex-conjugate pair (`Im(λ) ≥ 0`); overdamped/real roots
(`ζ ≥ 1`, `ω_d = 0`) are reported, not dropped, so the returned mode count is right. Sort
by ascending `|λ|`. Recombine the physical complex mode from the upper `nev` block of the
eigenvector: `φ_c = Φ q`; store real and imaginary parts via `expand_shape` on each part.

## Files to touch

- **`include/calculixpp/numerics/eigensolution.hpp`** — add `struct ComplexMode`
  (`std::complex<Real> eigenvalue; Real omega_d, omega_n, zeta, frequency, decay_rate;`
  `std::vector<Vec3> shape_real, shape_imag;`) and `struct ComplexEigenBasis`
  (`std::vector<ComplexMode> modes; Index n_free;`). Declare
  `ComplexEigenBasis extract_complex_modes(const EigenBasis& real_basis, const Damping&
  damping, std::size_t num_modes);` with a doc comment stating this is proportional-damping
  option (B) (not CalculiX CORIOLIS) and grounding the math.
- **`src/numerics/eigensolution.cpp`** — implement the C_r assembly (diagonal Rayleigh
  `α + β·ω_k²`, overridden by `2·ζ_k·ω_k` where modal ratios are set, via
  `Damping::ratio`), the `2·nev` companion build, the `numpp::linalg::eig` solve,
  conjugate-pair filtering / sort / post-processing, and the physical-mode recombination.
  Add the guarded small-problem dense `2n` oracle
  `extract_complex_modes_dense(K, C, M)` (cap `kDenseFallbackMaxDof`) for validation only.
- **`include/calculixpp/core/model.hpp`** — add `Procedure::ComplexFrequency` with a doc
  comment; add `int num_complex_modes{0};` and
  `enum class ComplexFrequencyType { Proportional, Coriolis, Flutter } complex_freq_type{Proportional};`
  (only `Proportional` implemented; `Coriolis`/`Flutter` parsed-then-rejected).
- **`src/io/inp_parser.cpp`** — remove `*COMPLEXFREQUENCY` from the deferred registry; add
  `begin_complex_frequency` (accept proportional / no keyword; `CORIOLIS` and `FLUTTER` →
  explicit `ParseError` "not yet implemented"), `complex_frequency_data` (data line = mode
  count → `num_complex_modes`), the data-dispatch hookup, and extend
  `validate_preceding_frequency` to require a preceding `*FREQUENCY` step for
  `Procedure::ComplexFrequency`.
- **`apps/cli/main.cpp`** — add the `Procedure::ComplexFrequency` branch: assemble `K,M`;
  `extract_modes`; build `Damping` from `model.rayleigh` / `model.modal_damping`;
  `extract_complex_modes`; print `MODE | DAMPED FREQ f_d | DAMPING RATIO ζ | DECAY RATE
  Re(λ) | STABILITY`.
- **`python/bindings.cpp`** — add `complex_frequency_result_dict` (numpy arrays:
  `eigenvalues_real`, `eigenvalues_imag`, `damped_frequencies`, `damping_ratios`,
  `omega_n`, and real+imag `mode_shapes`) and the `Procedure::ComplexFrequency` dispatch
  branch in `solve_model`.
- **tests** — parser tests (proportional `*COMPLEX FREQUENCY` sets the procedure + mode
  count; missing preceding `*FREQUENCY` errors; `CORIOLIS` rejected; `FLUTTER` rejected);
  numeric regression tests (Rayleigh SDOF + small cantilever vs the closed form; reduced
  path vs the dense `2n` oracle).

## Decisions

1. **Reduce onto the real basis (option B), not a full-scale complex solve.** Reuses the
   shipped mass-normalized basis, keeps the kernel a tiny dense `2·nev` `eig()`, and
   needs no SciPP sparse complex eigensolver. The problem stays `nev`-sized.
2. **Companion (standard) linearization** `A = [[0,I],[-Λ,-C_r]]`, `B = I`. Simplest form
   that feeds `numpp::linalg::eig` directly; real input, complex output.
3. **Analytical-only validation for this procedure.** CalculiX `*COMPLEX FREQUENCY`
   solves the gyroscopic problem, so its output would not match ζ/f_d from a proportional
   reduction. We therefore validate against the exact closed form and the dense oracle,
   and explicitly drop the CalculiX-fidelity requirement for `*COMPLEX FREQUENCY`.
4. **Reject `CORIOLIS`/`FLUTTER` at parse time.** CalculiX errors without one of them and
   `CORIOLIS` consumes a rotation body load we do not have. Rather than accept a deck we
   would mis-solve, reject with a clear "not yet implemented" message.
5. **Skew/imaginary-capable reduced-operator interface** so the future gyroscopic path is
   a plug-in, honoring the review's requirement that the "future skew G_r" claim be real.

## Risks

- **Modal truncation.** Exact for proportional damping (`C_r` diagonal, modes uncoupled).
  Non-proportional damping couples modes and its accuracy depends on `nev` — out of scope
  and documented as a follow-on.
- **`numpp::linalg::eig` on defective/near-defective pencils** from clustered frequencies
  (repeated `ω` → poor eigenvectors). Mitigate with the dense `2n` oracle cross-check and
  an eigenvector-residual tolerance. Optional NumPP hardening ask if fragile; not
  blocking for typical proportional decks.
- **Overdamped modes** (`ζ ≥ 1`, real negative `λ`, `ω_d = 0`). The `Im ≥ 0` filter and
  `|λ|` sort must still return the right count; handle real roots explicitly so
  `ω_d = 0` modes are reported, not dropped.
- **Single-file (non-multistep) decks.** The CLI re-extracts the real basis inline (as
  `*MODAL DYNAMIC` does) rather than reading a prior step's stored basis; the multi-step
  `validate_preceding_frequency` error still fires for true multi-step jobs.
- **Scope creep** into the Coriolis skew operator (`e_c3d.f:1171-1188`) and `FLUTTER`
  complex-force path — explicitly excluded; `CORIOLIS`/`FLUTTER` rejected at parse time.

## Validation (two tiers — analytical, no CalculiX fidelity claim)

1. **Analytical closed form (correctness proof).** Rayleigh-damped SDOF and a small
   damped cantilever: assert each reduced complex eigenvalue equals
   `λ_k = -ζ_k ω_k ± i ω_k √(1-ζ_k²)`, `ζ_k = (α/ω_k + β·ω_k)/2`, to ~1e-8 relative.
   This is exact (diagonal `C_r`), a genuine correctness proof — not a smoke test. Cross
   check the undamped limit (`α=β=0`) reproduces the real `*FREQUENCY` `ω`s with `ζ=0`.
2. **Internal cross-check.** For small decks (`n_free ≤ kDenseFallbackMaxDof`), compare
   the modal-reduced complex modes against the dense full `2n` state-space oracle to a
   tight relative tolerance, proving the reduction and linearization agree.

No CalculiX `*COMPLEX FREQUENCY` reference-deck comparison is performed for this
procedure: it is a different (gyroscopic) eigenproblem. When the true `CORIOLIS`
operator + rotor-speed input is added as a follow-on, that reference-fidelity requirement
comes back with it.

## Upstream asks

- **None blocking.** The path needs only a **dense** complex non-symmetric eigensolve on
  a `2·nev × 2·nev` real matrix — `numpp::linalg::eig()` already provides it. SciPP's
  missing **sparse** complex / non-symmetric eigensolver is **not** required.
- **Follow-on (not this change).** A SciPP sparse complex generalized eigensolver is
  needed for cyclic symmetry (`*CYCLIC SYMMETRY MODEL`), where the per-nodal-diameter
  pencil is genuinely complex at full FE scale — file against SciPP with that change.
- **Optional.** Confirm `numpp::linalg::eig` eigenvector normalization and behavior on
  repeated eigenvalues from clustered modes; file a NumPP hardening ask only if
  inverse-iteration recovery proves fragile. Not blocking.
