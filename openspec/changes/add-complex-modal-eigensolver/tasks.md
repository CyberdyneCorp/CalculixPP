# Tasks — add-complex-modal-eigensolver

## 1. Model & procedure enum (spec: modal-and-buckling-analysis)

- [x] 1.1 Add `Procedure::ComplexFrequency` to the `Procedure` enum in
  `include/calculixpp/core/model.hpp` with a doc comment stating it is the
  proportional-damping option-(B) path (NOT CalculiX CORIOLIS).
- [x] 1.2 Add model fields `int num_complex_modes{0};` and
  `enum class ComplexFrequencyType { Proportional, Coriolis, Flutter } complex_freq_type{Proportional};`
  to `model.hpp` (only `Proportional` is solved; `Coriolis`/`Flutter` are parsed then
  rejected).
- [x] 1.3 Build the module to enumerate every `switch`/dispatch site that must handle the
  new enum value (CLI procedure-name string, Python procedure-name string). **Blocker for
  §5 and §6.**

## 2. Parser (spec: input-deck-parsing, modal-and-buckling-analysis)

- [x] 2.1 Remove `*COMPLEXFREQUENCY` from the deferred-cards registry in
  `src/io/inp_parser.cpp` (~lines 454-456).
- [x] 2.2 Add card dispatch `else if (card_ == "*COMPLEXFREQUENCY") begin_complex_frequency(line);`.
- [x] 2.3 Implement `begin_complex_frequency`: set `procedure = ComplexFrequency`; accept
  no-keyword / proportional; **reject** `CORIOLIS` and `FLUTTER` with an explicit
  `ParseError` "not yet implemented" (grounded: CalculiX `complexfrequencys.f:86-92`
  errors without one of them, and `CORIOLIS` needs a rotation body load we do not have).
- [x] 2.4 Implement `complex_frequency_data` (data line = number of complex modes →
  `model_.num_complex_modes`) and add its case to the `data_line()` dispatch switch.
- [x] 2.5 Extend `validate_preceding_frequency` to require a preceding `*FREQUENCY` step
  for `Procedure::ComplexFrequency` (mirror the `ModalDynamic`/`SteadyState` check).
- [x] 2.6 Parser regression tests in `tests/test_parser.cpp`: proportional
  `*COMPLEX FREQUENCY` sets `Procedure::ComplexFrequency` and `num_complex_modes`; a
  missing preceding `*FREQUENCY` errors; `*COMPLEX FREQUENCY, CORIOLIS` is rejected;
  `*COMPLEX FREQUENCY, FLUTTER` is rejected.

## 3. Eigensolution header & reduced-operator interface (spec: eigensolution)

- [x] 3.1 Add `struct ComplexMode` (`std::complex<Real> eigenvalue; Real omega_d,
  omega_n, zeta, frequency, decay_rate; std::vector<Vec3> shape_real, shape_imag;`) and
  `struct ComplexEigenBasis { std::vector<ComplexMode> modes; Index n_free; };` to
  `include/calculixpp/numerics/eigensolution.hpp`.
- [x] 3.2 Declare
  `ComplexEigenBasis extract_complex_modes(const EigenBasis& real_basis, const Damping&
  damping, std::size_t num_modes);` with a doc comment: this is proportional-damping
  option (B), reduces onto the real basis, solves a `2·nev` companion via
  `numpp::linalg::eig`, and is explicitly NOT the CalculiX CORIOLIS gyroscopic problem.
- [x] 3.3 Define the reduced-operator assembly interface to carry (a) a symmetric damping
  block and (b) an optional **skew** block with **imaginary** ω-coupling (empty for
  option B), so the future gyroscopic `G_r` plugs in without redesign.

## 4. Eigensolution implementation (spec: eigensolution)

- [x] 4.1 Implement reduced `C_r` assembly in `src/numerics/eigensolution.cpp`: diagonal
  Rayleigh `α + β·ω_k²`, overridden by `2·ζ_k·ω_k` where modal ratios are set — reuse
  `Damping::ratio` for consistency with the modal-dynamics path.
- [x] 4.2 Build the `2·nev × 2·nev` real companion `A = [[0,I],[-Λ,-C_r]]` and solve with
  `numpp::linalg::eig(A)`.
- [x] 4.3 Filter conjugate pairs (`Im(λ) ≥ 0`), handle real/overdamped roots (`ω_d = 0`,
  not dropped), sort by ascending `|λ|`, and extract `ω_d`, `ω_n`, `ζ`, `f_d`, decay rate.
- [x] 4.4 Recombine the physical complex mode `φ_c = Φ q` from the upper `nev` block of
  each eigenvector, using `expand_shape` on the real and imaginary parts; store
  `shape_real` / `shape_imag`.
- [x] 4.5 Implement the guarded small-problem dense `2n` oracle
  `extract_complex_modes_dense(K, C, M)` (`[[0,I],[-M⁻¹K,-M⁻¹C]]`, cap
  `kDenseFallbackMaxDof`) for cross-validation only.

## 5. CLI (spec: modal-and-buckling-analysis)

- [x] 5.1 Add the `Procedure::ComplexFrequency` branch in `apps/cli/main.cpp`: assemble
  `K,M`; `extract_modes`; build `Damping` from `model.rayleigh`/`model.modal_damping`;
  call `extract_complex_modes`; print
  `MODE | DAMPED FREQ f_d | DAMPING RATIO ζ | DECAY RATE Re(λ) | STABILITY (stable if ζ>0)`.
  (Depends on §1.3.)

## 6. Python bindings (spec: python-bindings)

- [x] 6.1 Add `complex_frequency_result_dict` in `python/bindings.cpp` returning numpy
  arrays: `eigenvalues_real`, `eigenvalues_imag`, `damped_frequencies`, `damping_ratios`,
  `omega_n`, and real+imag `mode_shapes`.
- [x] 6.2 Add the `Procedure::ComplexFrequency` dispatch branch in `solve_model` near the
  `Frequency` case. (Depends on §1.3.)

## 7. Numeric regression tests (spec: modal-and-buckling-analysis, eigensolution)

- [x] 7.1 Rayleigh-damped SDOF: assert the reduced complex eigenvalue equals
  `λ = -ζω ± iω√(1-ζ²)` with `ζ = (α/ω + β·ω)/2`, ~1e-8 relative.
- [x] 7.2 Small Rayleigh-damped cantilever (C++ test + a Python pytest deck through the
  bindings): each computed damped frequency and damping ratio matches the closed form
  within the documented tolerance.
- [x] 7.3 Undamped limit (`α=β=0`) reproduces the real `*FREQUENCY` `ω`s with `ζ=0`.
- [x] 7.4 Cross-check the modal-reduced path against the dense `2n` oracle
  (`n_free ≤ kDenseFallbackMaxDof`) to a tight relative tolerance.

## 8. Docs & backlog (spec: modal-and-buckling-analysis)

- [x] 8.1 Update the procedure table and parser supported-cards list in `README.md` to
  add proportional `*COMPLEX FREQUENCY`, and note `CORIOLIS`/`FLUTTER` are rejected.
- [x] 8.2 Update `openspec/BACKLOG.md` item **1.6 / 2.3**: mark the proportional-damping
  complex-modes slice done; re-scope CalculiX `CORIOLIS` (gyroscopic skew operator +
  rotor-speed/axis input), `FLUTTER` (complex applied force), and cyclic symmetry (6.x)
  as explicit follow-ons. Do NOT mark `*COMPLEX FREQUENCY` fully done.
- [x] 8.3 Confirm no regression in the existing `*FREQUENCY` / `*MODAL DYNAMIC` /
  `*STEADY STATE DYNAMICS` / `*DYNAMIC` paths (run their tests, diff against main).
