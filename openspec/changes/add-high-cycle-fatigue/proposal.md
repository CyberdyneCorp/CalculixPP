# Add stress-based high-cycle-fatigue evaluation (`*HCF`)

## Why

The `high-cycle-fatigue` capability is a Phase-5 placeholder: its four requirements
gesture at "a fatigue/damage criterion over the stored field", a "critical location + life",
a missing-source error, and CalculiX reference fidelity — but nothing is implemented and
the `*HCF` card is unhandled.

The stock-CalculiX HCF path (`src/hcfs.f`, `src/combilcfhcf.f`) is a crack-growth
cumulative-damage engine coupled to the crack-propagation module — far outside a first
implementable slice and not what the baseline requirement needs. We can ship the useful,
self-contained core now: a **stress-amplitude vs. S-N (Basquin) endurance-curve** life
estimate over an already-computed stress field, reporting the worst-case location and its
cycles-to-failure. This is the standard textbook HCF criterion and is analytically
verifiable in closed form.

### Honest scope: stress-life (Basquin S-N), uniaxial-equivalent amplitude

This change ships a **stress-life** evaluation. The preceding step's recovered stress
field is treated as the cyclic **stress amplitude** at each point (the natural
interpretation of a modal / steady-state / static-amplitude result). A scalar
uniaxial-equivalent amplitude `S_a` is formed per point (signed von Mises by default), and
the Basquin S-N curve `S_a = a · N^b` is inverted for the cycles-to-failure
`N = (S_a / a)^(1/b)`. The worst-case location is the point of largest `S_a` (smallest
`N`). Multiaxial critical-plane methods, mean-stress corrections (Goodman/Gerber), and
cumulative Palmgren-Miner damage over a load spectrum are **deferred** (tasks §5, marked
`[~]`) — they layer on top of this same per-point amplitude+curve kernel without redesign.

## What Changes

- **New material data** — `SNCurve { a, b }` (Basquin coefficient + exponent) parsed from a
  `*FATIGUE` material card and stored on `Material`.
- **New procedure** — `Procedure::HighCycleFatigue`, selected by a `*HCF` step card.
- **New driver** — `numerics::evaluate_hcf(model)`: recover the preceding stress field via
  the existing `solve_linear_static` + nodal-stress path, form the per-node stress
  amplitude, invert the S-N curve per node, and return an `HcfReport` (per-node life,
  worst-case node id + location + amplitude + life). It **raises** a clear diagnostic when
  no valid preceding result source exists (no material S-N curve, or no elastic model to
  recover stress from).
- **Parser** — `*HCF` sets the procedure and the optional `CRITERION=` (only
  `SIGNED-VON-MISES` / `VON-MISES` accepted; others rejected); `*FATIGUE` fills the S-N
  curve.
- **CLI + Python bindings** — a worst-case summary line and an `hcf` result dict
  (per-node `life`/`amplitude`, worst-case node/location/life, criterion, S-N `a`/`b`).

## Capabilities

### Modified Capabilities

- **high-cycle-fatigue** — refine the four placeholder requirements to the shipped
  stress-life physics: `*HCF` requests a stress-amplitude/S-N evaluation over the preceding
  stress field; the worst-case location and its Basquin life are reported; a missing valid
  source (no S-N curve or no stress-producing model) is a diagnostic error; and the
  reference-fidelity requirement is replaced with an **analytical** closed-form
  Basquin-inversion fidelity requirement (the CalculiX crack-growth HCF output is a
  different quantity and would not match).

## Impact

- Code: `include/calculixpp/core/material.hpp`, `include/calculixpp/core/model.hpp`,
  `include/calculixpp/numerics/high_cycle_fatigue.hpp` (new),
  `src/numerics/high_cycle_fatigue.cpp` (new), `src/io/inp_parser.cpp`,
  `apps/cli/main.cpp`, `python/bindings.cpp`, `src/CMakeLists.txt`,
  `tests/CMakeLists.txt`, `tests/test_high_cycle_fatigue.cpp` (new).
- Behavior: `*HCF` now evaluates instead of being an unrecognized card; the
  crack-growth cumulative HCF, multiaxial critical-plane, mean-stress, and Miner-spectrum
  variants remain documented follow-ons.
- No breaking changes to any existing procedure. GPU backends remain optional; the
  evaluation runs on the CPU backend.
