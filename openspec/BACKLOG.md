# Deferred work backlog

Implementation items that were intentionally deferred (not dropped). Each is
already specified as a requirement in `openspec/specs/` — only the implementation
is outstanding — so they are tracked here rather than as a delta-less OpenSpec
change. Full context (with reasons) lives in the archived change that deferred
them, under `openspec/changes/archive/`.

## From Phase 2 (`2026-07-03-phase-2-nonlinear-statics-and-materials`)

| Item | What's missing | Enabler needed | Baseline spec |
|---|---|---|---|
| **1.6 `*STATIC, PERTURBATION`** | Linear perturbation about a preloaded base state with stress stiffening | Geometric stiffness `K_geo` / geometric nonlinearity (NLGEOM) | `static-analysis` |
| **2.3 `*CHANGE SOLID SECTION` (cross-step)** ✅ | RESOLVED (2026-07) by the multi-step engine — `*CHANGE SOLID SECTION` rebinds a material at the step boundary (sections accumulate, last-writer-per-element wins) and re-assembles per step. Validated on stock `changesolidsection` (E halved in step 2 -> displacement exactly doubles). `*CHANGE MATERIAL/PLASTIC` cross-step rebind is parsed and applies at the boundary too, but end-to-end plastic multi-step is out of the current linear slice (multi-step nonlinear deferred). | Multi-step analysis ✅ | `loads-and-boundary-conditions` |
| **2.4 `OP=MOD/NEW` (cross-step)** ✅ | RESOLVED (2026-07) by the multi-step engine — loads/BCs carry forward across steps and `OP=NEW` resets a load/BC type once per step. Validated on stock `beampfix` (step 2 `*CLOAD, OP=NEW` + `*BOUNDARY, FIXED`) matching the `.dat.ref` per-step displacements to rel-L2 ~1e-7. | Multi-step analysis ✅ | `loads-and-boundary-conditions` |
| **3.4 Shell / beam / membrane sections** | `*SHELL/*BEAM/*MEMBRANE SECTION` element formulations | Structural (reduced-dimension) element kernels | `element-sections` |
| **3.6 User element / distributions** | `*USER ELEMENT`, `*USER SECTION`, `*NODAL THICKNESS`, `*NORMAL`, `*DISTRIBUTION` | User-element interface + per-element distribution data | `element-sections` |
| **4.3 Hyperelasticity / foam (large strain)** | Large-strain `*HYPERELASTIC`/`*HYPERFOAM` with near-incompressible u/p | Finite-strain kinematics (deformation gradient, PK stresses) | `material-models` |
| **4.4 Creep / visco** | `*CREEP`, `*VISCO`, `*VALUES AT INFINITY` time-dependent integration | Time-dependent constitutive integration in the driver | `material-models` |
| **4.5 Deformation plasticity / Mohr-Coulomb / damage** | `*DEFORMATION PLASTICITY`, `*MOHR COULOMB` (+hardening), `*DAMAGE INITIATION` | (constitutive kernels; damage needs an evolution/coupling model) | `material-models` |

**Cross-cutting enablers** that unblock several rows above, worth planning as
their own work when picked up:
- **Geometric nonlinearity (NLGEOM):** finite-strain kinematics + geometric
  stiffness `K_geo`. Unblocks 1.6 and 4.3, and is a prerequisite for Phase-4
  linear buckling.
- **Multi-step analysis:** ✅ RESOLVED (2026-07). `io::parse_inp_steps` splits a
  deck into one `Model` per `*STEP...*END STEP` block (state/loads/BCs/procedure/
  active-mask carried forward; OP=NEW resets per step; `*CHANGE SOLID SECTION`/
  `*CHANGE MATERIAL` rebind at the boundary; a single-`*STEP` deck returns exactly
  the `parse_inp()` model). `numerics::solve_multistep_static` (src/numerics/
  multistep.cpp) solves each step INCREMENTALLY from the prior converged state, so a
  constant active set reproduces the single-step total-load solve exactly and a
  reactivated element is strain-free relative to the deformed config. `*BOUNDARY,
  FIXED` holds a DOF at its carried value. The single-step path is guarded and
  byte-identical. Unblocked **2.3** (`*CHANGE SOLID SECTION` cross-step) and **2.4**
  (`OP=MOD/NEW` cross-step), and the Phase-3 **5.1** cross-step element birth-death.
  Validated: stock two-step `beampfix`/`changesolidsection` decks + analytical
  superposition and birth-death. Still LINEAR-elastic mechanical only — multi-step
  thermal/coupled and multi-step nonlinear-material are deferred (a step in a
  multi-step list that is thermal/coupled or has a nonlinear material throws a clear
  message rather than silently mis-solving).

Deferred cards currently raise a clear `ParseError` naming the capability, so a
deck that uses them fails loudly rather than solving silently-wrong.

## From Phase 3 (`phase-3-thermal-and-contact`)

| Item | What's missing | Enabler needed | Baseline spec |
|---|---|---|---|
| **5.1 `*MODEL CHANGE, TYPE=ELEMENT` (cross-step)** ✅ | RESOLVED (2026-07) by the multi-step engine. An element removed in a later step contributes nothing thereafter (per-step active mask), and one re-added is strain-free relative to the *deformed* geometry at reactivation (the incremental driver enters it into `K_ff` only for the steps it is active). Validated analytically (two stacked cubes: removed-then-re-added carries only its post-reactivation stress) and via the stock two-step `beampfix` reference. | Multi-step analysis ✅ | `model-change` |
| **5.2 `*MODEL CHANGE, TYPE=CONTACT PAIR` (activation)** ✅ | RESOLVED (2026-07). The multi-step engine is contact-aware: `Model::with_active_contact_pairs()` filters the pairs per step from the accumulated `contact_pair_changes` (a pair active unless its last matching ADD/REMOVE is a REMOVE), and `solve_multistep_contact` solves each step's changed model through the nonlinear contact driver. A pair inactive in step 1 / active in step 2 changes the response (validated: driven two-block stack's base reaction jumps ~0 → large). | Multi-step engine ✅ + contact operator ✅ | `model-change` |
| **1.3 Surface-to-surface MORTAR contact** | Dual/standard-Lagrange segment integration over overlapping master/slave face segments on non-matching meshes (`contactmortar.c`/`multimortar.c`). Node-to-surface penalty contact is complete; mortar is research-level and NOT faked — `build_contact_pairs`/`build_thermal_contact_pairs` reject `TYPE=SURFACE TO SURFACE` with a deferral message. | Mortar segment integration + condensable Lagrange multipliers. | `contact-search`, `contact` |

Model-change note: `*MODEL CHANGE, TYPE=ELEMENT` (5.1) and `TYPE=CONTACT PAIR` (5.2)
are both now FULLY implemented, within-step and cross-step (via the multi-step
engine). The multi-step engine routes a contact deck through the nonlinear driver per
step with each step's active contact-pair set.

Contact note (2026-07): NODE-TO-SURFACE PENALTY CONTACT is now implemented (tasks
1.1, 1.2, 1.4 normal, 1.5, 1.6, 2.1 n2s, 2.2). On top of the geometric core (the
uniform-grid `fem::MasterGrid` proximity search and the `fem::project_onto_face` /
`face_frame_at` node-to-surface projection with signed gap + outward normal/tangent
frame), the penalty operator (`fem/contact.{hpp,cpp}`) evaluates the normal
pressure-overclosure law (`contact_force_weight`: hard/linear/exponential from
`*SURFACE BEHAVIOR`), runs an active-set update, and adds the contact FORCE to the
Newton residual and the consistent penalty `k n⊗n` TANGENT (through the constraint
transform) inside `numerics::solve_nonlinear_static` each iteration — projecting the
DEFORMED slave against the DEFORMED master face (`fem::face_eval_deformed`). The
parser reads `*CONTACT PAIR` (TYPE=NODE TO SURFACE), `*SURFACE INTERACTION`, and
`*SURFACE BEHAVIOR`; a contact deck routes to the nonlinear driver (CLI + Python via
`Model::has_contact()`), and a deck with NO contact is byte-for-byte the pre-contact
path. VALIDATED analytically (two stacked cubes meeting only through a `*CONTACT
PAIR`: the interface transmits the load — base reaction == applied load — with small
penetration ~F/(4*kappa); tests/test_contact.cpp + python) and end to end on the
stock `contact1` n2s deck (converges, right sign/order; its exact NLGEOM+exponential
match is not gated).

TANGENTIAL COULOMB FRICTION is now implemented (tasks 1.4 tangential + 2.3). The
parser reads `*FRICTION` (mu + optional stick stiffness) on the `*SURFACE
INTERACTION` (`Friction{mu,stick,has}`, `Model::contact_friction`). The penalty
operator adds a tangential traction to the same Newton residual+tangent as the normal
contact: an elastic STICK predictor `t_trial = -stick * g_t` (g_t the tangential slip
of the slave vs its master footprint, `I - n⊗n` projected) is return-mapped onto the
Coulomb cone `|t| <= mu*p` (`fem::friction_return_map`) — inside the cone STICKS,
outside SLIPS (traction capped at mu*p along the slip direction). The consistent
regime tangent (`friction_tangent_matrix`: `stick*P` stuck, the small
`mu*p/|g_t|(P - ŝ⊗ŝ)` slip-direction tangent when slipping) is scattered through the
constraint transform, and the auto stick stiffness is sized to the master shear scale
so the stick<->slip transition converges. Only penetrating pairs (g<0) get friction;
a frictionless deck (no `*FRICTION`) is byte-for-byte the normal-only path.
VALIDATED analytically (tests/test_contact.cpp `test_friction_return_map` +
`test_friction_operator_stick_slip`: at the operator level, a slave driven with a
known penetration+slide gives the elastic stick force `-stick*slip` below the cone and
the exact SATURATED `mu*p` traction above it — the defining stick->slip transition,
plus the mu-scaling of the cap; end to end `test_friction_end_to_end_converges` + the
Python `test_contact_friction_stick_slip_python` two-block friction solve converges and
keeps the normal equilibrium).

THERMAL CONTACT is now implemented (task 2.4). `*GAP CONDUCTANCE` (coefficient h_c)
and `*GAP HEAT GENERATION` (generated flux q) parse onto the `*SURFACE INTERACTION`;
the thermal driver reuses the SAME node-to-surface search (`build_thermal_contact_pairs`
+ `add_thermal_contact`) to add the conductance coupling `Q = h_c A (T_slave − Σ N_i
T_master_i)` — a symmetric conductance block over the interface temperature DOFs — plus
the gap heat source, into both the steady and transient thermal operators. A thermal
deck with a `*CONTACT PAIR` routes to `solve_heat_transfer`. VALIDATED analytically:
two cubes touching only through the gap reproduce the series-resistance flux
`Q = ΔT/(2/k + 1/h)` and the interface drop `Q/h` (q = h_c ΔT), with energy
conservation and gap-heat balance (tests/test_thermal_contact.cpp + python).

CONTACT OUTPUT + MODIFIERS is now implemented (task 2.5). `fem::recover_contact`
reports per-slave-node CSTR (status/pressure/gap/tangential traction) filled into
`StaticFields.contact`; the `.dat`/`.frd` writers emit the CSTR block / CONTACT dataset
and the Python result dict a `contact` list. `*CONTACT FILE`/`PRINT`/`OUTPUT` are
accepted; `*CLEARANCE` sets an initial gap the operator adds to the geometric gap
(`g_eff = g_geo + clearance`), opening/pre-loading the interface regardless of the mesh.

PHASE 3 FINALIZED (2026-07): the contact workstream is complete except **1.3
surface-to-surface MORTAR** segment integration (dual/standard Lagrange) — the single
remaining contact piece, research-level, rejected loudly (`TYPE=SURFACE TO SURFACE`),
NOT faked. Node-to-surface penalty contact (spatial search + projection, normal
pressure-overclosure, Coulomb friction, thermal gap conductance/heat generation, CSTR
output + `*CLEARANCE`, and `*MODEL CHANGE` contact-pair activation) is complete,
reachable from Python (`solve`/`summary` + a `contact` CSTR list in the result dict),
and covered by the parametrized `CONTACT_CORPUS` regression plus the C++ operator
suites. Validated against analytical references (two-block equilibrium + penalty
penetration `F/(4κ)`, series-resistance gap flux, stick↔slip return-mapping); the stock
contact decks (NLGEOM + CalculiX's exact exponential law + `*EQUATION` MPCs) are a
harder match and not gated.

Resolved (2026-07): **Coupled temperature-displacement (4.1 monolithic + 4.2
iterated staggered)**. `*COUPLED TEMPERATURE-DISPLACEMENT[, SOLUTIONS=MONOLITHIC|
STAGGERED]` now solves the two fields together. MONOLITHIC (`solve_monolithic_once`)
assembles ONE 4-DOF/node system `[[K_uu,K_uT],[K_Tu,K_TT]][u;T]=[f_ext;q]`, reusing
the reduced mechanical (`assemble_linear_static`) and thermal (`assemble_conduction`)
operators plus the new thermal-strain coupling block `K_uT = fem::element_thermal_
coupling` (∫ Bᵀ D {α,α,α,0,0,0} Nⱼ dV), reduced through the mechanical constraint
transform and thermal Dirichlet map, solved once. STAGGERED (`solve_coupled_
staggered`) is a Gauss-Seidel outer loop with a temperature convergence check. The
mechanical→thermal two-way term is plastic-dissipation (Taylor-Quinney) heating
(`Model::taylor_quinney`, `Model::plastic_dissipation_heat`, fed by the per-element
committed eqplastic now reported from `solve_nonlinear_static`): with it absent the
4-DOF tangent is block-triangular so BOTH schemes equal the sequential one-way solve
EXACTLY (validated: monolithic == staggered, temperature ~1e-9, stress ~1e-6, plus
the analytical bar sigma_xx=-Eα T_bar); with it present both schemes converge to the
same joint (T,u) state and staggered iterates (validated: a stretched-past-yield
element deposits dissipation heat, outer count > 1, temperature rises above the
conducted baseline). See tasks 4.1 / 4.2 (tests/test_thermal_stress.cpp +
python test_coupled_monolithic_equals_staggered).

Resolved (2026-07): **Temp-dependent material tables** (`k(T)`/`c(T)`/`α(T)`
piecewise-linear interpolation with clamping, Picard-iterated in the thermal
kernels; validated vs the Kirchhoff transform) and **Integration-point HFL
output** (`q=-k∇T` recovered at Gauss points, extrapolated to nodes for the
`.frd`, emitted as `*EL PRINT HFL` to the `.dat`; validated vs `oneel20cf.dat.ref`
and the exact Fourier flux). See tasks 6.1 / 6.2.

Resolved (2026-07): **3.4 Cavity radiation** (`*RADIATE ...,CR` gray-body
surface-to-surface exchange). Geometric view factors `F_ij` by direct
double-area quadrature over a refined face sampling (row-normalized to the
`Σ_j F_ij = 1` enclosure summation rule), the gray-body radiosity solve
`(I-(1-eps)F) J = eps sigma T^4`, and the Newton-linearized per-patch net flux
`Q_i(T)` + tangent `dQ_i/dT_k` assembled into the thermal system alongside
surface-to-ambient radiation (steady + transient). Validated analytically: the
summation rule, the known cube view factors (Hottel), the two-plate gray-body
net flux `q = sigma(T1^4-T2^4)/(1/e1+1/e2-1)`, and a Python end-to-end
parallel-plate energy balance. See task 3.4 (`fem::build_cavity` /
`cavity_heat_flow`, `numerics::add_cavity_radiation`).

## From Phase 4 (`phase-4-dynamics-and-eigenproblems`)

| Item | What's missing | Enabler needed | Baseline spec |
|---|---|---|---|
| **1.1 Sparse shift-invert Lanczos** | The scalable eigensolver: a shift-invert Lanczos / subspace iteration driving a SciPP sparse factorization `(K-σM)⁻¹` (via `scipp::sparse::spsolve`) as the Krylov operator, so `*FREQUENCY` scales past the dense `O(n_free³)` path. The DENSE generalized path (Cholesky of M -> `eigh`) is implemented and is the validated oracle; it is what `*FREQUENCY` uses today. | A **sparse GENERALIZED eigensolver** — filed as **[SciPP#12](https://github.com/CyberdyneCorp/SciPP/issues/12)** (`eigsh`, shift-invert Lanczos on the sparse factorization). Also **[NumPP#138](https://github.com/CyberdyneCorp/NumPP/issues/138)**: the DENSE `eigh`/`solve` the current path uses scale **O(n^4)** (measured 0.07s->95s over 132->972 DOF; a 1275-DOF `*FREQUENCY` takes 111s vs numpy's 0.19s) - the eigen analog of #10. Until fixed, the beam8f eigen-validation tests are marked `slow` and excluded from the default gate. | `eigensolution`, `linear-algebra-and-solvers` |
| **1.3 Sturm-sequence count check** | Missed/starved-mode verification via the inertia (Sturm sequence) of the shifted factorization. Not needed for the dense path (`eigh` returns the exact, complete spectrum); required only once the sparse Lanczos path lands. | Sparse Lanczos path (1.1) + inertia from the sparse factorization. | `eigensolution` |
| **1.6 / 2.3 Complex / damped modes** | The complex eigenproblem for `*COMPLEX FREQUENCY` (damping / friction-induced instability). | Complex eigensolve + damping/friction operator assembly. | `eigensolution`, `modal-and-buckling-analysis` |
| **2.2 `*BUCKLE`** | Two-step prestress + geometric-stiffness eigenproblem `(K + λ K_geo) x = 0`. The eigensolution engine it consumes is in place. | Geometric stiffness `K_geo` (see the NLGEOM cross-cutting enabler). | `modal-and-buckling-analysis` |
| **3.x `*DYNAMIC, EXPLICIT` + `*GREEN`** | Explicit central-difference direct dynamics (element wave-speed critical time step) and the `*GREEN` Green-function step. (IMPLICIT direct HHT-α dynamics — 3.1 / 3.2 / 3.3 — are shipped, see the resolved note below.) | Element critical-time-step estimate (explicit); unit-excitation response basis (`*GREEN`). | `dynamic-analysis` |
| **4.3 Base-motion time histories / 4.4 `*GREEN`** | Full `*BASE MOTION` support-excitation TIME HISTORIES driving the transient (the effective-load builder `base_motion_load` and its `p_k = -Γ_k` projection are shipped; wiring a prescribed base-acceleration amplitude into the modal-dynamic march lands with a follow-up) and the `*GREEN` step. | — (builds on the shipped modal engine). | `dynamic-analysis` |
| **6.x Cyclic symmetry** | Cyclic-symmetry sector complex eigenproblem (`*CYCLIC SYMMETRY MODEL`) + nodal-diameter mode selection (`*SELECT CYCLIC SYMMETRY MODES`). | Complex cyclic-symmetry constraints (a deferred `constraints` capability) + the complex eigensolve (1.6 / 2.3). | `constraints`, `modal-and-buckling-analysis` |
| **2.4 Preceding-frequency-step chaining** | The RUNTIME requirement that a `*MODAL DYNAMIC` / `*STEADY STATE DYNAMICS` step reuse the eigenbasis of a PRECEDING `*FREQUENCY` step across steps. Enforced today for the single-deck case (the modal procedures assemble K/M and extract the basis from the same model, so the preceding extraction is implicit); a cross-step "consume the previous step's stored basis" requires the multi-step dynamics driver. | Multi-step dynamics driver that stores + hands the extracted basis to the following step. | `dynamic-analysis`, `modal-and-buckling-analysis` |

Resolved (2026-07): **Phase-4 LEAD slice — mechanical mass matrix + eigensolution
engine + `*FREQUENCY`.** `fem::element_mass` (consistent `M_e = ∫ρ Nᵀ N dV`, the
mechanical analog of `element_capacitance`) + `element_mass_lumped` (row-sum) reuse
the shared shape functions / Gauss rules; `fem::assemble_mass` assembles the global
mass on the SAME free-DOF numbering / constraint transform as `K` (`*MASS` point
masses add a diagonal nodal mass). `numerics::extract_modes` solves the generalized
symmetric eigenproblem `K x = λ M x` for the lowest N modes by the DENSE path
(Cholesky `M = L Lᵀ` -> standard symmetric `A = L⁻¹ K L⁻ᵀ` -> `numpp::linalg::eigh`
-> back-transform + mass-normalize + ascending sort), with a spectral shift for
rigid-body handling and `numerics::participation` for participation factors / modal
effective mass. `*FREQUENCY` parses (number of eigenvalues), auto-dispatches from the
CLI and Python, and reports `f = sqrt(λ)/(2π)` + mode shapes. Validated against the
stock CalculiX `beam8f` `.dat.ref` (10 eigenvalues, natural frequencies, participation
factors, effective + total mass — all < 1e-4 rel) AND an analytical 2-DOF spring-mass
chain (exact). See tasks 1.1-1.4 / 2.1. The scalable sparse shift-invert Lanczos is
deferred (row 1.1 above) pending a sparse generalized eigensolver in NumPP/SciPP.

Resolved (2026-07): **Phase-4 modal superposition dynamics — `*MODAL DYNAMIC` +
`*STEADY STATE DYNAMICS` + damping (tasks 1.5 / 4.1 / 4.2 / 4.3).**
`numerics/modal_dynamics.{hpp,cpp}` projects the mass-normalized eigenbasis into the
decoupled modal SDOFs (`project_modal_system` — per-mode ω_k, ω_k², ζ_k;
`ModalSystem::project_load` maps a physical load to `p_k = φ_kᵀ f`). `*MODAL DYNAMIC`
(`modal_dynamic`) integrates each SDOF `q̈_k + 2ζ_kω_k q̇_k + ω_k² q_k = p_k(t)` with the
EXACT piecewise-linear-load recurrence (Nigam-Jennings — zero algorithmic damping, the
natural period `2π/ω` reproduced exactly) and recombines `u = Σ q_k φ_k`. `*STEADY STATE
DYNAMICS` (`steady_state_response`/`_sweep`) evaluates the complex modal transfer function
`q_k = p_k/(ω_k²-Ω²+2iζ_kω_kΩ)` over a frequency sweep -> per-node amplitude/phase.
Damping: Rayleigh `C = αM + βK` -> `ζ_k = (α/ω_k + βω_k)/2`, `*MODAL DAMPING` per-mode
ratios override; `base_motion_load` builds the support-excitation effective load
`f_eff = -M r` (modal forcing `p_k = -Γ_k`). Parser adds `*MODAL DYNAMIC`, `*STEADY STATE
DYNAMICS`, `*DAMPING`, `*MODAL DAMPING`, `*BASE MOTION`; CLI + Python auto-dispatch (modal
dynamic returns the displacement time history, steady-state the amplitude/phase sweep).
Validated ANALYTICALLY (no clean in-scope reference deck — the stock modal-dynamic decks
use coupling reference nodes / reduced-integration beams out of scope): undamped step
response `(F/k)(1-cos ωt)` to < 1e-7 over many periods, damped step matches the closed-form
under-damped solution, steady-state resonant peak `(F/k)·Q` with `Q=1/(2ζ)` and correct
half-power bandwidth `Δω ≈ 2ζω`; plus an end-to-end check that a steady-state sweep on the
real `beam8f` C3D8 mesh resonates at a natural frequency of the full model. `*GREEN` and
base-motion time histories remain deferred (rows above).

Resolved (2026-07): **Phase-4 direct-integration dynamics — `*DYNAMIC` implicit HHT-α
(tasks 3.1 / 3.2 / 3.3).** `numerics/direct_dynamics.{hpp,cpp}` (`direct_dynamic`)
integrates `M a + C v + K u = f(t)` in PHYSICAL coordinates by the implicit HHT-α
(Hilber-Hughes-Taylor) generalization of Newmark-β. Solving for `a_{n+1}` directly
(Newmark correctors substituted into the HHT equilibrium at `t_{n+1+α}`), the effective
dynamic operator is `A_eff = M + (1+α)γΔt C + (1+α)βΔt² K`; the LINEAR path forms it once
and re-solves each step through `numerics::solve_reduced` (ComputeBackend/SciPP — no GPU),
the NONLINEAR path (`step_nonlinear`, task 3.2) runs a per-step Newton on `a_{n+1}`,
folding the material-point tangent (`fem::assemble_material_tangent` + `fem::add_contact`)
into the effective tangent so inertia + nonlinear material/contact are solved together.
`HhtParams::from_alpha` maps the α knob (`α ∈ [-1/3, 0]`) to `β=(1-α)²/4, γ=½-α`; each
frame carries kinetic `½vᵀMv` / strain `½uᵀKu` / total energy and `DirectReport.energy_drift`
(task 3.3). Parser adds `*DYNAMIC[, ALPHA=][, NLGEOM][, DIRECT]` (data line `dt, t_end`);
CLI + Python auto-dispatch (returns the displacement/velocity/acceleration/total-energy
history + `hht_alpha`/`energy_drift`/`newton_iterations`). Validated ANALYTICALLY (no clean
in-scope reference deck — the stock `*DYNAMIC` decks preload a static step / use out-of-scope
elements): step response `(F/k)(1-cos ωt)` to < 2e-4, free vibration `u0 cos ωt` (natural
period `2π/ω`) to < 5e-6, UNDAMPED α=0 free vibration CONSERVES total energy over 20 periods
(drift < 1e-10), HHT α=-0.3 DISSIPATES energy, Rayleigh damping decays the amplitude, and the
nonlinear path reproduces the linear path to < 1e-9 (`test_direct_dynamics`, 6 cases; pytest
`test_dynamic_*`, deck path). `*DYNAMIC, EXPLICIT` (central difference) and `*GREEN` remain
deferred (rows above).

Resolved (2026-07): **Phase-4 substructure / superelement generation — `*SUBSTRUCTURE
GENERATE` Craig-Bampton / Guyan reduction (tasks 5.1 / 5.2 / 5.3).**
`numerics/substructure.{hpp,cpp}` (`generate_substructure` / `reduce_substructure`)
partitions the free DOFs into boundary `b` (the retained/master DOFs of `*RETAINED NODAL
DOFS`) and interior `i` (condensed). Static (Guyan) reduction forms the reduced stiffness
as the Schur complement `K̂ = K_bb − K_bi K_ii⁻¹ K_ib`; when a mass matrix is requested
(`*SUBSTRUCTURE MATRIX OUTPUT, MASS=YES` / `*MATRIX ASSEMBLE, MASS=YES` / a nonzero mode
count) Craig-Bampton builds the constraint modes `Ψ = −K_ii⁻¹ K_ib` and the fixed-interface
normal modes `Φ` (from the eigensolution engine on the interior-restricted `K_ii x = λ M_ii
x`), then projects K/M onto `T = [[I,0],[Ψ,Φ]]` giving the reduced operators (leading b×b
Guyan block + diagonal `ω_k²` / unit-mass modal block + the retained/modal mass coupling).
The interior condensation `K_ii⁻¹` routes through NumPP dense solve on the ComputeBackend
(no GPU). Parser adds `*SUBSTRUCTURE GENERATE`, `*RETAINED NODAL DOFS` (node/nset, dof
range), `*SUBSTRUCTURE MATRIX OUTPUT` (STIFFNESS=/MASS=), `*MATRIX ASSEMBLE`; CLI writes the
reduced matrices in the reference `*MATRIX TYPE=STIFFNESS/MASS` lower-triangular row format
(`<out>.mtx`), Python `solve()` returns `k_reduced`/`m_reduced` + retained-DOF labels +
`modal_omega`. An empty retained-DOF set / an SPC'd-or-slave retained DOF raises a clear
error (spec scenarios). Validated against the stock CalculiX `substructure.inp` reduced
stiffness — its `.dat.ref` `*MATRIX TYPE=STIFFNESS` block (60 retained DOFs, Guyan) matches
to < 1e-6 rel (C++ `test_substructure`, pytest `test_substructure_reduced_stiffness_matches_calculix`) —
and, for Craig-Bampton, the reduced model's lowest eigenfrequencies approximate the FULL
model's `*FREQUENCY` to < 2% on a cantilever column (`test_craig_bampton_approximates_full`,
pytest `test_substructure_craig_bampton_reduces_full_model`). The DENSE interior condensation
is `O(n_interior³)` (same scalability note as the eigensolution engine — a sparse
Schur-complement path is a SciPP follow-up). Only cyclic symmetry (6.x) and the complex /
buckling paths remain deferred (rows above).

Resolved (2026-07): **Phase-4 FINALIZE — parser blocked-card actionable errors + API
parity + docs (tasks 7.1 / 7.2 / 7.4 / 7.5).** The deferred Phase-4 procedure cards
(`*BUCKLE`, `*COMPLEX FREQUENCY`, `*CYCLIC SYMMETRY MODEL`, `*SELECT CYCLIC SYMMETRY
MODES`, `*GREEN`) now raise a clear, ACTIONABLE `ParseError` naming the deferral and the
enabler each waits on (reject_card deferred-map in src/io/inp_parser.cpp) instead of a
generic "unsupported card" — a deck using them fails loudly rather than silently
mis-solving (regression: test_parser `test_blocked_phase4_cards_reject_clearly` + pytest
`test_api_parity_phase4` blocked-card sweep). A Phase-4 API-parity coverage test
(`test_api_parity_phase4`) asserts every one of the five shipped procedures
(`*FREQUENCY` / `*DYNAMIC` / `*MODAL DYNAMIC` / `*STEADY STATE DYNAMICS` /
`*SUBSTRUCTURE GENERATE`) exposes its documented result payload and that `summary()`
reports each without solving. README gains a "Python — dynamics (Phase 4)" walkthrough
(runnable HHT-α `*DYNAMIC` SDOF example + the modal/steady-state/damping variants) and a
Phase-4 header status blurb. The rows below remain the only outstanding Phase-4 work.
