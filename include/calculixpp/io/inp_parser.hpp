#pragma once
#include <stdexcept>
#include <string>

#include "calculixpp/core/model.hpp"

// Abaqus-style input-deck parser (spec: input-deck-parsing — Phase-2 card set).
// Supported cards:
//   Model:    *NODE, *ELEMENT (C3D4/C3D10, C3D8/C3D20/C3D6/C3D15 + C3D8R/C3D20R,
//             connector SPRINGA/1/2, MASS, DASHPOTA/1/2), *NSET, *ELSET, *SURFACE.
//   Material: *MATERIAL, *ELASTIC, *DENSITY, *PLASTIC (+HARDENING=), *CYCLIC
//             HARDENING, *HYPERELASTIC (N=1 neo-Hookean), *USER MATERIAL,
//             *DEPVAR, *RATEDEPENDENT.
//   Sections: *SOLID SECTION, *SPRING, *MASS, *DASHPOT.
//   Loads/BC: *BOUNDARY, *CLOAD, *DLOAD (P<face>/GRAV/CENTRIF), *DSLOAD,
//             *AMPLITUDE, OP=MOD|NEW.
//   Step:     *STEP, *STATIC (incl. DIRECT + increment data line, NLGEOM accepted),
//             *CONTROLS, *TIME POINTS, *CHANGE MATERIAL/PLASTIC/SOLID SECTION.
//             MULTIPLE *STEP...*END STEP blocks parse into a per-step model list via
//             parse_inp_steps (multi-step analysis); parse_inp collapses to one model.
//   Constr.:  *EQUATION, *MPC, *RIGID BODY, *COUPLING (+*KINEMATIC/*DISTRIBUTING),
//             *DISTRIBUTING COUPLING, *TIE.
// Output-request cards (*NODE PRINT, *EL PRINT, *NODE FILE, *EL FILE, ...) are
// accepted and ignored. DEFERRED Phase-2 cards (*HYPERFOAM, *CREEP, *VISCO,
// *MOHR COULOMB, *DAMAGE INITIATION, *DEFORMATION PLASTICITY, *SHELL/*BEAM/
// *MEMBRANE SECTION) raise an actionable ParseError naming the deferral rather
// than being silently ignored. Any other keyword raises a generic ParseError.
namespace cxpp::io {

struct ParseError : std::runtime_error {
  ParseError(int line, const std::string& msg)
      : std::runtime_error("input deck line " + std::to_string(line) + ": " + msg),
        line(line) {}
  int line;
};

Model parse_inp(const std::string& text);
Model parse_inp_file(const std::string& path);

// Parse a deck into its per-*STEP models (spec: multi-step analysis). Each
// *STEP...*END STEP block becomes one fully-accumulated Model, carrying the prior
// step's loads/BCs/procedure/controls and *MODEL CHANGE active-element mask forward;
// OP=NEW resets a load/BC type once per step and *CHANGE SOLID SECTION / *CHANGE
// MATERIAL rebinds at the boundary. Global (pre-step) mesh/materials/constraints are
// shared across the returned models. A single-*STEP deck returns exactly one model
// equal to parse_inp()'s, so the single-step path is unchanged. The step-loop driver
// (numerics::solve_multistep_static) solves each step from the previous converged state.
std::vector<Model> parse_inp_steps(const std::string& text);
std::vector<Model> parse_inp_steps_file(const std::string& path);

}  // namespace cxpp::io
