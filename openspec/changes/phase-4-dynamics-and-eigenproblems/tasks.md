# Phase 4 — Dynamics & Eigenproblems

Implements baseline dynamics/eigen/substructure physics plus the one new capability
(`eigensolution`). Each task cites the capability it satisfies. Depends on
`phase-1-foundation` and `phase-2-nonlinear-statics-and-materials`.

## 1. Eigensolution engine (spec: eigensolution [NEW])

- [x] 1.1 Generalized symmetric eigenpair extraction via NumPP. DENSE generalized path
      implemented (numerics/eigensolution.{hpp,cpp}): Cholesky of M -> standard symmetric
      A = L⁻¹ K L⁻ᵀ -> numpp::linalg::eigh -> back-transform φ = L⁻ᵀ z. Mass matrix
      fem::element_mass (consistent) / element_mass_lumped (row-sum) + fem::assemble_mass
      (same free-DOF numbering / constraint transform as K). Validated: stock CalculiX
      *FREQUENCY deck beam8f (C3D8 cantilever) reproduces all 10 reference eigenvalues to
      < 1e-4 rel, AND an analytical 2-DOF spring-mass chain (λ = (k/m)(3∓√5)/2 exact).
      [~] The SCALABLE shift-invert Lanczos on SciPP's sparse factorization
      ((K-σM)⁻¹ via spsolve) is NOT implemented — NumPP/SciPP expose no sparse GENERALIZED
      eigensolver, so it is a SciPP follow-up (like #10). The dense path is O(n_free³);
      *FREQUENCY runs on it today. See openspec/BACKLOG.md Phase-4 row.
- [x] 1.2 Mass normalization (φᵀ M φ = 1) and ascending-order modal basis. eigh returns
      ascending eigenvalues; each eigenvector is mass-normalized. Verified by the
      generalized-residual test (K φ = λ M φ) and φᵀ M φ = 1 to 1e-9.
- [x] 1.3 Spectral shift + rigid-body handling: extract_modes(..., shift σ) solves
      (K - σ M); eigenvalues returned un-shifted (λ = μ + σ). The dense M-Cholesky path is
      robust for rigid-body / near-zero eigenvalues (M is always PD). Validated on a fully
      unconstrained hex (6 rigid-body modes at ~0 extracted without failure; a nonzero
      shift returns identical un-shifted eigenvalues). [~] Sturm-sequence count check on a
      shifted factorization is deferred with the sparse Lanczos path (dense eigh returns
      the exact spectrum, so no count check is needed for the dense path).
- [x] 1.4 Participation factors and modal effective mass (numerics::participation) for an
      excitation direction: Γ_k = φ_kᵀ M r, m_eff_k = Γ_k². Validated against beam8f's
      reference PARTICIPATION FACTORS + EFFECTIVE MODAL MASS columns and TOTAL (X-direction
      matched to < 1e-4 rel; total effective mass 0.8509573E-07 reproduced).
- [x] 1.5 Modal superposition projection of operators and loads (numerics/modal_dynamics
      .{hpp,cpp}: project_modal_system builds the reduced modal operators — per-mode ω_k,
      ω_k² (unit modal mass), and ζ_k — from the mass-normalized basis + a Damping model;
      ModalSystem::project_load maps a physical free-DOF load to the modal load p_k = φ_kᵀ f).
      Validated: p_k matches the direct dot φ_kᵀ f on a 2-DOF chain (test_modal_dynamics),
      and end-to-end through the bindings (test_regression modal/steady-state).
- [~] 1.6 Complex / damped eigenproblem path (*COMPLEX FREQUENCY). BLOCKED: NumPP/SciPP
      expose a DENSE non-symmetric eigensolver (numpp::linalg::eig) but the complex-damped
      modal path needs the quadratic eigenproblem (M,C,K) recast to a state-space complex
      eigensolve plus the damping/friction operator assembly — not yet built. Not faked.
      The real eigensolution engine it builds on is in place. See openspec/BACKLOG.md
      Phase-4 "Complex / damped modes" row.

## 2. Eigenfrequency, buckling, complex frequency (spec: modal-and-buckling-analysis)

- [x] 2.1 `*FREQUENCY` natural frequencies + mode shapes. Parser: `*FREQUENCY` card sets
      Procedure::Frequency + num_eigenvalues (first data-line field). CLI auto-dispatch
      (apps/cli/main.cpp) prints eigenvalue + cyclic frequency per mode; Python auto-dispatch
      (solve() returns eigenvalue/omega/frequency/mode_shape/participation/effective_mass).
      Reports natural frequencies f = ω/(2π) with ω = sqrt(λ). Validated: beam8f .dat.ref
      (10 modes, freqs to < 1e-4 rel) via C++ (test_eigen) + Python (test_regression.py:
      test_frequency_beam8f_matches_calculix / test_frequency_participation_beam8f).
- [~] 2.2 `*BUCKLE` two-step prestress + geometric-stiffness eigenproblem. BLOCKED on the
      geometric stiffness K_geo, which requires geometric nonlinearity (NLGEOM) — a
      cross-cutting enabler NOT yet built (see Phase-2 BACKLOG). The eigensolution engine
      it would consume (solve (K + λ K_geo) x = 0 via the generalized path with a spectral
      shift for the non-positive spectrum) is in place, so *BUCKLE lands as soon as K_geo
      exists. A fake K_geo is deliberately NOT produced. See openspec/BACKLOG.md Phase-4
      "*BUCKLE" row.
- [~] 2.3 `*COMPLEX FREQUENCY` damped / friction-induced modes. Deferred with the complex
      eigenproblem path (task 1.6) — same blocker (complex quadratic/state-space
      eigensolve + damping/friction operator). Not faked.
- [~] 2.4 Preceding-frequency-step requirement enforced for modal procedures.
      io::validate_preceding_frequency(steps) enforces the spec requirement across a
      MULTI-STEP job: a *MODAL DYNAMIC / *STEADY STATE DYNAMICS step with no preceding
      *FREQUENCY step throws a clear error naming the offending step
      (modal-and-buckling-analysis — "modal-superposition procedures SHALL require a
      preceding *FREQUENCY step in the same job"). Validated in test_parser
      (test_preceding_frequency_required / _satisfied). [~] In the SINGLE-STEP modal decks
      the basis is extracted from the same model at dispatch (the *FREQUENCY extraction is
      implicit, so the requirement is satisfied by construction); wiring a modal step to
      CONSUME a preceding step's STORED basis needs the multi-step dynamics driver (the
      multi-step engine is linear-static only today) — deferred, see openspec/BACKLOG.md
      Phase-4 row (2.4). *COMPLEX FREQUENCY's preceding-*FREQUENCY case lands with the
      complex slice (1.6 / 2.3, deferred).

## 3. Direct-integration dynamics (spec: dynamic-analysis)

- [x] 3.1 Newmark / HHT-α integrator (numerics/direct_dynamics.{hpp,cpp}: direct_dynamic).
      Implicit HHT-α (Hilber-Hughes-Taylor) generalization of Newmark-β integrates
      M a + C v + K u = f(t) in PHYSICAL coordinates. Solving for a_{n+1} directly (Newmark
      correctors substituted into the HHT equilibrium at t_{n+1+α}), the effective dynamic
      operator is A_eff = M + (1+α)γΔt C + (1+α)βΔt² K; the linear path forms it ONCE and
      re-solves each step through numerics::solve_reduced (ComputeBackend/SciPP — no GPU
      needed). Parser: `*DYNAMIC[, ALPHA=][, NLGEOM][, DIRECT]` -> Procedure::Dynamic, data
      line "dt, t_end" (shared with *MODAL DYNAMIC). CLI + Python auto-dispatch (solve()
      returns time/displacement/velocity/acceleration/total_energy history + energy_drift).
      Validated ANALYTICALLY: an undamped SDOF under a suddenly-applied step load reproduces
      u(t)=(F/k)(1-cos ωt) to < 2e-4 (2nd-order, fine step), and free vibration from an
      initial displacement reproduces u0 cos ωt to < 5e-6 (natural period 2π/ω) —
      test_direct_dynamics (C++) + test_dynamic_sdof_step_response (pytest, deck path).
- [x] 3.2 Nonlinear direct dynamics through `nonlinear-solution-control`. The per-step
      Newton loop (step_nonlinear) iterates on a_{n+1}: each iteration assembles the
      material-point tangent + internal force via fem::assemble_material_tangent (and
      fem::add_contact when *CONTACT PAIR present), folds the inertia/damping into the
      effective tangent (K_eff = c_m M + c_c C + k_scale K_tan), forms the HHT dynamic
      residual (inertia + damping + internal force), and solves through solve_reduced. Routed
      by `*DYNAMIC, NLGEOM` or a nonlinear material / contact in the model; material history
      committed per accepted step. Validated: the nonlinear path reproduces the linear path
      on a linear-elastic SDOF to < 1e-9 (test_direct_dynamics
      test_nonlinear_reproduces_linear — the driver's nonlinear-reproduces-linear contract).
- [x] 3.3 Numerical-damping (α) control + ENERGY check. HhtParams::from_alpha maps the α
      knob (α ∈ [-1/3, 0]) to β=(1-α)²/4, γ=½-α; DirectTimePoint carries per-step kinetic
      (½vᵀMv) / strain (½uᵀKu) / total energy, and DirectReport.energy_drift is the peak
      relative total-energy drift over the run. Validated ANALYTICALLY: an UNDAMPED α=0
      (trapezoidal Newmark) free vibration CONSERVES total energy over 20 periods (drift
      < 1e-10, test_energy_conservation_undamped); HHT α=-0.3 DISSIPATES energy over the same
      free run (test_hht_numerical_damping); Rayleigh C=αM+βK decays the free-vibration
      amplitude (test_rayleigh_damping_decay). Deck path exercised by
      test_dynamic_energy_conserved_undamped / test_dynamic_hht_alpha_dissipates (pytest).
      [~] EXPLICIT central-difference *DYNAMIC, EXPLICIT (element wave-speed critical step)
      is NOT implemented — the implicit HHT scheme is the unconditionally-stable default; the
      explicit path is deferred (see openspec/BACKLOG.md Phase-4 row).

## 4. Modal & steady-state dynamics (spec: dynamic-analysis)

- [x] 4.1 `*MODAL DYNAMIC` superposition using the eigensolution basis. numerics::
      modal_dynamic integrates each decoupled modal SDOF q̈_k + 2ζ_kω_k q̇_k + ω_k² q_k =
      p_k(t) with the EXACT piecewise-linear-load recurrence (Nigam-Jennings — zero
      algorithmic damping, natural period reproduced exactly), recombining u = Σ q_k φ_k.
      Parser: `*MODAL DYNAMIC` (dt, t_end) -> Procedure::ModalDynamic. CLI + Python
      auto-dispatch (solve() returns time/displacement history arrays). Validated
      ANALYTICALLY: undamped step response = (F/k)(1-cos ωt) to < 1e-7 over many periods,
      and a damped SDOF step response matches the closed-form under-damped solution
      (test_modal_dynamics + test_regression test_modal_dynamic_sdof_step_response).
- [x] 4.2 `*STEADY STATE DYNAMICS` harmonic response over a frequency sweep. numerics::
      steady_state_response / steady_state_sweep solve the complex modal transfer function
      q_k = p_k/(ω_k²-Ω²+2iζ_kω_kΩ) and recombine to complex physical amplitude ->
      magnitude/phase. Parser: `*STEADY STATE DYNAMICS` (f_lo, f_hi, n_points) ->
      Procedure::SteadyStateDynamics; CLI + Python auto-dispatch (returns frequency/
      amplitude/phase arrays). Validated ANALYTICALLY: resonant peak = (F/k)·Q with
      Q=1/(2ζ) and correct half-power bandwidth Δω ≈ 2ζω (test_modal_dynamics); and
      end-to-end on the real beam8f C3D8 mesh the sweep resonates at a natural frequency of
      the full model (test_regression test_steady_state_sdof_resonance /
      test_steady_state_beam8f_resonates_at_an_eigenfrequency).
- [x] 4.3 Damping (Rayleigh + modal) and base motion. Damping::ratio maps Rayleigh
      C = αM + βK to ζ_k = (α/ω_k + βω_k)/2, with explicit *MODAL DAMPING ratios overriding
      per mode. base_motion_load builds the support-excitation effective load f_eff = -M r
      so the modal forcing is p_k = -Γ_k (participation factor). Parser: `*DAMPING`
      (ALPHA/BETA), `*MODAL DAMPING` (per-mode ratios), `*BASE MOTION` (DOF/magnitude).
      Validated: Rayleigh->modal mapping + override (test_modal_dynamics
      test_rayleigh_mapping); base-motion load projects to -Γ_k
      (test_projection_and_base_motion).
- [ ] 4.4 `*GREEN` Green-function step

## 5. Substructure / superelement (spec: substructure-generation)

- [x] 5.1 Fixed-interface normal modes from eigensolution + constraint modes over retained
      DOFs. numerics/substructure.{hpp,cpp}: reduce_substructure partitions the free DOFs
      into boundary b (retained/master, from *RETAINED NODAL DOFS) and interior i
      (condensed). Constraint modes Ψ = -K_ii⁻¹ K_ib (NumPP dense solve). Fixed-interface
      normal modes Φ come from the eigensolution engine — extract_modes on the
      interior-restricted K_ii x = λ M_ii x (restrict_system builds a fresh identity
      transform on the interior numbering so the reused engine stays self-consistent).
- [x] 5.2 Craig-Bampton reduced stiffness/mass assembly (`*SUBSTRUCTURE GENERATE`,
      `*RETAINED NODAL DOFS`). Guyan static reduction gives the reduced stiffness Schur
      complement K̂ = K_bb - K_bi K_ii⁻¹ K_ib; Craig-Bampton (when mass/modes requested)
      projects K/M onto T = [[I,0],[Ψ,Φ]] -> leading b×b Guyan block + diagonal ω_k² /
      unit-mass modal block + retained/modal mass coupling. Parser: *SUBSTRUCTURE GENERATE
      -> Procedure::Substructure, *RETAINED NODAL DOFS (node/nset, first_dof, last_dof) ->
      model.retained_dofs (dedup, SORTED=NO declaration order). An empty retained set / an
      SPC'd-or-slave retained DOF throws a clear error (spec scenarios). CLI + Python
      auto-dispatch (solve() returns k_reduced/m_reduced + retained-DOF labels/modal_omega).
      Validated: Craig-Bampton reduced eigenfrequencies approximate the FULL model's
      *FREQUENCY to < 2% on a cantilever column (test_substructure
      test_craig_bampton_approximates_full; pytest
      test_substructure_craig_bampton_reduces_full_model); plus an analytical 3-DOF Guyan
      Schur complement and the no-interior identity case.
- [x] 5.3 Reduced / global matrix export (`*SUBSTRUCTURE MATRIX OUTPUT`, `*MATRIX ASSEMBLE`).
      format_substructure_matrix / write_substructure_matrix write the reduced K (and M)
      in the reference CalculiX *MATRIX TYPE=STIFFNESS/MASS lower-triangular row format
      (row i has i+1 entries); the CLI writes <out>.mtx. *SUBSTRUCTURE MATRIX OUTPUT
      (STIFFNESS=/MASS=) selects the exported matrices; *MATRIX ASSEMBLE, MASS=YES also
      requests the reduced mass (the assembled global K/M are reachable via fem::assemble_*).
      VALIDATED against a stock CalculiX .dat.ref: substructure.inp (Guyan, 60 retained
      DOFs) reduced-stiffness *MATRIX block matches to < 1e-6 rel (max entry error 3.5e-11)
      — test_substructure test_reference_substructure_deck (C++) + pytest
      test_substructure_reduced_stiffness_matches_calculix.
      [~] The interior condensation K_ii⁻¹ is the DENSE O(n_interior³) path (same
      scalability note as the eigensolution engine); a sparse Schur-complement is a SciPP
      follow-up. See openspec/BACKLOG.md Phase-4 row.

## 6. Cyclic symmetry (spec: constraints)

- [~] 6.1 Cyclic-symmetry sector complex eigenproblem (`*CYCLIC SYMMETRY MODEL`). BLOCKED
      on two deferred enablers: (a) complex cyclic-symmetry constraints (a `constraints`
      capability not yet built — the per-nodal-diameter phase-shifted tie between the two
      sector cut boundaries), and (b) the complex eigensolve (task 1.6 / 2.3, deferred —
      NumPP/SciPP expose no complex generalized eigensolver). Not faked. See
      openspec/BACKLOG.md Phase-4 "Cyclic symmetry" row.
- [~] 6.2 Nodal-diameter mode selection (`*SELECT CYCLIC SYMMETRY MODES`). Deferred with
      6.1 (it selects which nodal diameters the sector eigenproblem of 6.1 solves).

## 7. Parser, bindings, validation

- [x] 7.1 Extend the deck parser to the Phase-4 card set. `*FREQUENCY` (number of
      eigenvalues), `*DYNAMIC` (dt, t_end; ALPHA=/NLGEOM/DIRECT params for HHT direct
      integration), `*MODAL DYNAMIC` (dt, t_end), `*STEADY STATE DYNAMICS`
      (f_lo, f_hi, n_points), `*DAMPING` (ALPHA/BETA Rayleigh), `*MODAL DAMPING` (per-mode
      ratios), and `*BASE MOTION` (DOF/magnitude) landed; `*SUBSTRUCTURE GENERATE`,
      `*RETAINED NODAL DOFS` (node/nset, first_dof, last_dof), `*SUBSTRUCTURE MATRIX
      OUTPUT` (STIFFNESS=/MASS=), and `*MATRIX ASSEMBLE` (MASS=) landed for the
      substructure slice. The BLOCKED/deferred procedure cards (`*BUCKLE`,
      `*COMPLEX FREQUENCY`, `*CYCLIC SYMMETRY MODEL`, `*SELECT CYCLIC SYMMETRY MODES`,
      `*GREEN`) now raise a clear, ACTIONABLE ParseError naming the deferral and the
      enabler each waits on (reject_card deferred-map, src/io/inp_parser.cpp) — a deck
      using them fails loudly rather than silently mis-solving. Covered by
      test_parser test_blocked_phase4_cards_reject_clearly (C++) + the pytest
      test_api_parity_phase4 blocked-card sweep.
- [x] 7.2 Extend Python bindings + API-parity coverage test. `*FREQUENCY`, `*DYNAMIC`,
      `*MODAL DYNAMIC`, and `*STEADY STATE DYNAMICS` auto-dispatch through
      solve()/solve_text(): frequency returns eigenvalue/omega/frequency/mode_shape/
      participation/effective_mass/total_effective_mass; DIRECT dynamic returns time/
      displacement/velocity/acceleration/total_energy history + hht_alpha/energy_drift/
      newton_iterations; modal dynamic returns time/displacement history + omega/zeta;
      steady-state returns frequency/amplitude/phase sweep arrays; *SUBSTRUCTURE GENERATE
      returns k_reduced/m_reduced + retained-DOF labels (retained_node/retained_comp) +
      modal_omega. API-parity coverage test test_api_parity_phase4 asserts every one of
      the five Phase-4 procedures' documented result keys is present, that summary()
      reports each procedure WITHOUT solving, and that the blocked cards raise the
      actionable deferral error. Remaining (blocked) procedures land with their slices.
- [~] 7.3 Phase-4 reference corpus. `beam8f` (stock *FREQUENCY C3D8 cantilever, 1e-4 rel on
      eigenvalues/frequencies/participation) and `substructure` (stock *SUBSTRUCTURE
      GENERATE Guyan, reduced-stiffness *MATRIX block to < 1e-6 rel) added. Remaining Phase-4 decks
      land with their procedures.
- [x] 7.4 pytest regression. beam8f eigenvalues + participation matched to reference (see
      test_frequency_beam8f_matches_calculix / test_frequency_participation_beam8f). Modal
      superposition dynamics validated ANALYTICALLY end-to-end through the bindings:
      test_modal_dynamic_sdof_step_response (undamped step = (F/k)(1-cos ωt)),
      test_steady_state_sdof_resonance (peak (F/k)·Q, half-power bandwidth), and
      test_steady_state_beam8f_resonates_at_an_eigenfrequency (real C3D8 mesh resonates at a
      natural frequency). DIRECT dynamics validated end-to-end through the bindings:
      test_dynamic_sdof_step_response (HHT step response = (F/k)(1-cos ωt)),
      test_dynamic_energy_conserved_undamped (α=0 energy field bounded), and
      test_dynamic_hht_alpha_dissipates (α=-0.3 decays vs undamped). C++:
      test_modal_dynamics (5 analytical cases) + test_direct_dynamics (6 analytical cases:
      step response, free-vibration period, α=0 energy conservation drift<1e-10, HHT
      dissipation, Rayleigh decay, nonlinear-reproduces-linear). SUBSTRUCTURE validated
      against the stock CalculiX substructure.inp reduced stiffness (.dat.ref *MATRIX block,
      Guyan 60 retained DOFs, < 1e-6 rel) via C++ (test_substructure, 6 cases:
      Guyan 3-DOF analytical, no-interior identity, empty-retained error, export format,
      Craig-Bampton approximates full model < 2%, reference deck) and pytest
      (test_substructure_reduced_stiffness_matches_calculix,
      test_substructure_empty_retained_errors,
      test_substructure_craig_bampton_reduces_full_model). Preceding-*FREQUENCY enforcement
      covered by test_parser (validate_preceding_frequency); blocked-card actionable
      errors covered by test_parser test_blocked_phase4_cards_reject_clearly + pytest
      test_api_parity_phase4. The only regression decks still absent are for the DEFERRED
      procedures (buckle/complex/cyclic/explicit/green) — they land with those slices.
- [x] 7.5 Update README/docs with a modal + dynamics example. Added a Phase-4 README
      section "Python — frequency & substructure" with a runnable *FREQUENCY summary and a
      full Craig-Bampton *SUBSTRUCTURE GENERATE example (solve() returning k_reduced/
      m_reduced + labels), plus a new "Python — dynamics (Phase 4)" section walking a
      runnable *DYNAMIC (HHT-α) SDOF example (returning time/displacement/velocity/
      acceleration/total_energy + hht_alpha/energy_drift) and describing the *MODAL DYNAMIC
      and *STEADY STATE DYNAMICS / *DAMPING variants + their result payloads. Added a
      Phase-4 status blurb to the header (mass matrix + eigensolution engine + the five
      procedures, validated decks, and the deferred cards with actionable errors); the
      roadmap row stays "in progress" (core shipped; buckle/complex/cyclic/explicit/green
      deferred). Mermaid diagrams unchanged (ascii labels, valid).
