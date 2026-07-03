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
| **5.2 `*MODEL CHANGE, TYPE=CONTACT PAIR` (activation)** | Actually enabling/disabling a contact pair per step. The card is parsed and each change stored on `Model::contact_pair_changes`, but nothing consumes it yet. | The contact workstream (contact search + assembly) plus multi-step step handling. | `model-change` |

Model-change note: `*MODEL CHANGE, TYPE=ELEMENT` is now FULLY implemented, both
within-step and cross-step (via the multi-step engine — see 5.1 above). Only
`*MODEL CHANGE, TYPE=CONTACT PAIR` (5.2) remains parse-and-store, pending the
contact workstream.

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
