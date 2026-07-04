# Tasks — add-high-cycle-fatigue

## 1. Material & model data (spec: high-cycle-fatigue, material-models)

- [x] 1.1 Add `struct SNCurve { Real a; Real b; bool empty(); }` to
  `include/calculixpp/core/material.hpp` (Basquin `S_a = a·N^b`, `b < 0`) and an
  `std::optional<SNCurve> sn_curve;` field on `Material`.
- [x] 1.2 Add `Procedure::HighCycleFatigue` to the `Procedure` enum in
  `include/calculixpp/core/model.hpp` with a doc comment: stress-life over a preceding
  stress field (NOT the CalculiX crack-growth HCF).
- [x] 1.3 Add an `enum class FatigueCriterion { SignedVonMises, VonMises }` and a
  `FatigueCriterion hcf_criterion{SignedVonMises};` field on `Model`.

## 2. Parser (spec: input-deck-parsing, high-cycle-fatigue)

- [x] 2.1 Dispatch `*HCF` → `begin_hcf(line)`: set `procedure = HighCycleFatigue`; read
  optional `CRITERION=`; accept only `SIGNED-VON-MISES`/`VON-MISES`, reject others with a
  `ParseError`.
- [x] 2.2 Dispatch `*FATIGUE` (material card) → `begin_fatigue` / `fatigue_data`: store the
  data line `a, b` as the current material's `SNCurve`; error if no `*MATERIAL` is active.
- [x] 2.3 Parser regression tests in `tests/test_parser.cpp`: `*HCF` sets the procedure and
  criterion; `*FATIGUE` fills the S-N curve; an unknown `CRITERION=` raises.

## 3. Numerics driver (spec: high-cycle-fatigue)

- [x] 3.1 New header `include/calculixpp/numerics/high_cycle_fatigue.hpp`: `HcfReport`
  (per-node `life`, `amplitude`; worst-case `node_id`, `location`, `amplitude`, `life`) and
  `HcfReport evaluate_hcf(const Model& model);`.
- [x] 3.2 New `src/numerics/high_cycle_fatigue.cpp`: recover the preceding stress field
  (`solve_linear_static`), form the per-node scalar amplitude via the model criterion
  (signed/plain von Mises), invert the material S-N curve `N = (S_a/a)^(1/b)` per node, and
  fill the worst-case (largest amplitude → smallest life).
- [x] 3.3 Diagnostic errors: no `SNCurve` on any material → throw; no `effective_elastic`
  model (nothing to recover stress from) → throw. Message names the missing source.
- [x] 3.4 Register `numerics/high_cycle_fatigue.cpp` in `src/CMakeLists.txt`.

## 4. CLI + Python bindings (spec: high-cycle-fatigue)

- [x] 4.1 CLI dispatch in `apps/cli/main.cpp`: a `*HCF` deck prints the worst-case node,
  amplitude, and cycles-to-failure.
- [x] 4.2 Python `hcf_result_dict` + procedure-name string + `solve_model` branch:
  per-node `life`/`amplitude` arrays, worst-case node/location/life, criterion, `sn_a`/`sn_b`.

## 5. Deferred follow-ons (documented, not implemented)

- [~] 5.1 Multiaxial critical-plane criteria (Findley/Dang Van) — needs a plane search over
  the stress tensor; the scalar uniaxial-equivalent amplitude ships first. Reason: out of
  scope for the first slice; the per-point kernel is the shared substrate.
- [~] 5.2 Mean-stress correction (Goodman/Gerber/Soderberg) — needs the mean stress in
  addition to the amplitude; the amplitude-only Basquin ships first. Reason: requires a
  min/max or mean field the current single-field result does not carry.
- [~] 5.3 Cumulative Palmgren-Miner damage over a load spectrum — needs multiple blocks
  with cycle counts. Reason: single-amplitude life ships first; Miner sums per-block
  `n_i/N_i` on top of the same S-N inversion.
- [~] 5.4 CalculiX crack-growth cumulative HCF fidelity (`hcfs.f`/`combilcfhcf.f`) —
  different physics (crack propagation), not portable in this slice. Reason: replaced by an
  analytical Basquin-inversion fidelity check.

## 6. Validation

- [x] 6.1 `tests/test_high_cycle_fatigue.cpp`: a uniaxial bar under a known stress
  amplitude `S` with an S-N curve `(a, b)` yields the closed-form `N = (S/a)^(1/b)` to
  `< 1e-6` relative; the worst-case node is the highest-amplitude node. Register in
  `tests/CMakeLists.txt`.
- [x] 6.2 Missing-source diagnostics: an `*HCF` deck with no `*FATIGUE` curve, and one with
  no elastic material, each throw.
- [x] 6.3 `ctest --test-dir build` and the Python tests stay green (no regressions).
