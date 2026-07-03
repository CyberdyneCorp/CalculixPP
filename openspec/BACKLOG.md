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
| **2.3 `*CHANGE MATERIAL/PLASTIC/SOLID SECTION` (cross-step)** | Property redefinition at a step boundary (within-step rebind is done) | Multi-step analysis / step-history loop | `loads-and-boundary-conditions` |
| **2.4 `OP=MOD/NEW` (cross-step)** | Load accumulation/reset *across* steps (within-step reset is done) | Multi-step analysis / step-history loop | `loads-and-boundary-conditions` |
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
- **Multi-step analysis:** a step loop carrying load/state history across `*STEP`
  blocks. Unblocks 2.3 and 2.4.

Deferred cards currently raise a clear `ParseError` naming the capability, so a
deck that uses them fails loudly rather than solving silently-wrong.
