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
