#include "calculixpp/io/inp_parser.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace cxpp::io {
namespace {

std::string trim(std::string s) {
  const auto is_sp = [](unsigned char c) { return std::isspace(c) != 0; };
  while (!s.empty() && is_sp(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
  while (!s.empty() && is_sp(static_cast<unsigned char>(s.back()))) s.pop_back();
  return s;
}

std::string upper(std::string s) {
  for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

// Uppercase and strip spaces — for keyword and parameter-name matching.
std::string normalize(const std::string& s) {
  std::string out;
  for (char c : s)
    if (!std::isspace(static_cast<unsigned char>(c)))
      out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
  return out;
}

std::vector<std::string> split_fields(const std::string& line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string tok;
  while (std::getline(ss, tok, ',')) out.push_back(trim(tok));
  // Drop a single trailing empty field from a trailing comma.
  if (!out.empty() && out.back().empty()) out.pop_back();
  return out;
}

char last_nonspace(const std::string& s) {
  for (auto it = s.rbegin(); it != s.rend(); ++it)
    if (!std::isspace(static_cast<unsigned char>(*it))) return *it;
  return '\0';
}

// A logical record: its source line number and the fields on it.
struct Line {
  int number;
  std::string text;
};

// Merge continuation lines (a data line ending in ',' continues) and drop
// comments/blanks. Keyword lines terminate a continuation.
std::vector<Line> preprocess(const std::string& text) {
  std::vector<std::string> raw;
  std::stringstream ss(text);
  std::string l;
  while (std::getline(ss, l, '\n')) {
    if (!l.empty() && l.back() == '\r') l.pop_back();
    raw.push_back(l);
  }
  std::vector<Line> out;
  for (std::size_t i = 0; i < raw.size();) {
    const std::string t = trim(raw[i]);
    // Skip blanks, comment lines (**), and a bare "*" — a keyword line with no
    // name carries no card and appears as a separator in some stock decks.
    if (t.empty() || t.rfind("**", 0) == 0 || t == "*") {
      ++i;
      continue;
    }
    std::string merged = raw[i];
    const int start = static_cast<int>(i) + 1;
    while (last_nonspace(merged) == ',' && i + 1 < raw.size()) {
      const std::string nxt = trim(raw[i + 1]);
      if (nxt.empty() || nxt.rfind("**", 0) == 0 || nxt.front() == '*') break;
      merged += raw[i + 1];
      ++i;
    }
    out.push_back({start, merged});
    ++i;
  }
  return out;
}

double to_double(const std::string& s, int line) {
  try {
    std::size_t pos = 0;
    const double v = std::stod(s, &pos);
    return v;
  } catch (...) {
    throw ParseError(line, "expected a number, got '" + s + "'");
  }
}

Index to_index(const std::string& s, int line) {
  try {
    return static_cast<Index>(std::stol(s));
  } catch (...) {
    throw ParseError(line, "expected an integer, got '" + s + "'");
  }
}

// Parser state machine.
class Parser {
 public:
  Model run(const std::string& text) {
    const std::vector<Line> lines = preprocess(text);
    for (const Line& ln : lines) {
      const std::vector<std::string> fields = split_fields(ln.text);
      if (fields.empty()) continue;
      if (fields.front().rfind("*", 0) == 0) {
        begin_card(fields, ln.number);
      } else {
        data_line(fields, ln.number);
      }
    }
    flush_sets();
    return std::move(model_);
  }

 public:
  // Parse a deck into its list of *STEP models (spec: multi-step analysis). Each
  // *STEP...*END STEP block is snapshotted at *END STEP as one fully-accumulated
  // Model (loads/BCs/procedure/controls/active-element mask carried forward from the
  // prior step, OP=NEW resetting per step; sections accumulate so *CHANGE SOLID
  // SECTION rebinds at the boundary). Global (pre-step) data — mesh, nsets/elsets/
  // surfaces, materials, amplitudes, physical constants, constraints — is shared, so
  // the returned models differ only in their step-varying state. A single-*STEP deck
  // returns exactly one model that is byte-for-byte the parse_inp() model.
  std::vector<Model> run_steps(const std::string& text) {
    const std::vector<Line> lines = preprocess(text);
    for (const Line& ln : lines) {
      const std::vector<std::string> fields = split_fields(ln.text);
      if (fields.empty()) continue;
      if (fields.front().rfind("*", 0) == 0)
        begin_card(fields, ln.number);
      else
        data_line(fields, ln.number);
    }
    // A deck may omit *END STEP on its last step (or have a single implicit block);
    // if the current step accumulated any step content beyond what was snapshotted,
    // capture it. `in_step_` tracks whether we are inside an unterminated block.
    if (in_step_) snapshot_step();
    if (steps_.empty()) steps_.push_back(model_);  // no *STEP at all -> one model
    // Flush the shared sets/surfaces into every step model's mesh.
    for (Model& m : steps_) flush_sets_into(m);
    return std::move(steps_);
  }

 private:
  Model model_;
  // Per-step snapshots captured at *END STEP (see run_steps). Each is a copy of the
  // accumulated model_ at that boundary; sets are flushed into their meshes at the end.
  std::vector<Model> steps_;
  bool in_step_ = false;  // true between *STEP and *END STEP
  std::string card_;                                        // normalized keyword
  std::unordered_map<std::string, std::string> params_;     // normalized key -> value
  std::string cur_material_;
  std::string cur_amplitude_;
  // Hardening kind carried from the *PLASTIC card's HARDENING= parameter onto the
  // data lines (CalculiX default ISOTROPIC).
  Plastic::Hardening cur_hardening_ = Plastic::Hardening::Isotropic;
  // OP=NEW on a load/BC card resets prior loads of that type; per CalculiX it
  // takes effect on the FIRST card of that type in the step only. These flags
  // record that the first-card reset has already been consumed for the step.
  // (Single-step model: the reset clears loads accumulated earlier in the same
  // deck. Cross-step OP semantics need multi-step step handling — deferred.)
  bool op_seen_cload_ = false;
  bool op_seen_dload_ = false;
  bool op_seen_boundary_ = false;
  bool op_seen_cflux_ = false;
  bool op_seen_dflux_ = false;
  // True once a *HEAT TRANSFER procedure card is seen, so *BOUNDARY dof 11 routes to
  // a prescribed temperature and *DFLUX to a surface heat flux. (spec: heat-transfer.)
  bool heat_step_ = false;
  std::unordered_map<std::string, std::vector<Index>> nsets_, elsets_;
  std::vector<Surface> surfaces_;  // built incrementally; flushed at end
  int cur_surface_ = -1;           // index into surfaces_, or -1 outside a surface

  // *SURFACE INTERACTION / *SURFACE BEHAVIOR / *CONTACT PAIR (spec: contact). The
  // current interaction name (so a following *SURFACE BEHAVIOR attaches to it) and the
  // pair whose slave/master come on its data line. -1 / empty outside those cards.
  std::string cur_interaction_;
  int cur_contact_pair_ = -1;  // index into model_.contact_pairs, or -1

  // *DESIGN RESPONSE / *OBJECTIVE (spec: design-optimization). The name from the
  // current *DESIGN RESPONSE NAME= (attached to the response its data line defines).
  std::string cur_design_response_;

  // Connector (discrete) elements: SPRINGA/SPRING1/SPRING2, MASS, DASHPOTA/1/2.
  // These are not isoparametric solids, so they live outside the Mesh. We keep
  // their connectivity keyed by element id and by elset so a following *SPRING /
  // *MASS / *DASHPOT card (which references the elset) can attach the property.
  enum class ConnKind { SpringA, Spring1, Spring2, Mass, DashpotA, Dashpot1, Dashpot2 };
  struct ConnElem {
    ConnKind kind;
    std::vector<Index> nodes;  // 1 or 2 node ids
  };
  std::unordered_map<Index, ConnElem> conn_elems_;                 // elem id -> def
  std::unordered_map<std::string, std::vector<Index>> conn_elsets_;  // elset -> ids
  // For SPRING1/SPRING2 the DOFs are given on the *SPRING data line before the
  // stiffness; parsed into these while reading the current *SPRING/*DASHPOT card.
  int spring_dof1_ = 0, spring_dof2_ = 0;
  bool spring_dofs_read_ = false;

  // *EQUATION: -1 means "next data line is the term count"; >=0 counts terms still to
  // read for the current equation (which may span continuation lines).
  int eq_terms_remaining_ = -1;
  Equation cur_equation_;
  // *COUPLING: the coupling being built (kind set by a following *KINEMATIC /
  // *DISTRIBUTING modifier). -1 outside a coupling.
  int cur_coupling_ = -1;
  // *TIE: the tie being built (surfaces come on the data line).
  Tie tie_pending_;
  bool tie_has_pending_ = false;

  // *MODEL CHANGE: the current card's TYPE (Element / ContactPair) and ADD/REMOVE
  // sense, carried onto the data lines. (spec: model-change — element/contact birth-death.)
  enum class ModelChangeType { Element, ContactPair };
  ModelChangeType mc_type_ = ModelChangeType::Element;
  bool mc_add_ = true;  // true = ADD (activate), false = REMOVE (deactivate)

  std::string param(const std::string& key) const {
    const auto it = params_.find(key);
    return it == params_.end() ? std::string{} : it->second;
  }

  // Apply OP= on a load/BC card (spec: loads-and-boundary-conditions — load
  // accumulation). OP=MOD (default) keeps accumulating; OP=NEW resets the prior
  // loads of that type, but only for the first such card in the step (`seen`
  // guards the once-per-step rule). Returns whether the reset fired.
  bool apply_op_reset(bool& seen, int line) {
    const std::string op = param("OP");
    if (op.empty() || op == "MOD") return false;
    if (op != "NEW") throw ParseError(line, "OP must be MOD or NEW, got '" + op + "'");
    if (seen) return false;  // NEW takes effect on the first card of the type only
    seen = true;
    return true;
  }

  void begin_card(const std::vector<std::string>& fields, int line) {
    card_ = normalize(fields.front());
    params_.clear();
    for (std::size_t i = 1; i < fields.size(); ++i) {
      const std::string& f = fields[i];
      const auto eq = f.find('=');
      if (eq == std::string::npos) {
        params_[normalize(f)] = "";
      } else {
        params_[normalize(f.substr(0, eq))] = upper(trim(f.substr(eq + 1)));
      }
    }
    // Cards that carry no data or need setup at declaration time.
    static const std::vector<std::string> ignored = {
        "*NODEPRINT",   "*ELPRINT",     "*NODEFILE",    "*ELFILE",
        "*CONTACTFILE", "*CONTACTPRINT", "*CONTACTOUTPUT", "*OUTPUT", "*RESTART",
        "*HEADING",     "*SECTIONPRINT"};
    if (card_ == "*STEP") {
      begin_step();
    } else if (card_ == "*ENDSTEP") {
      end_step();
    } else if (card_ == "*STATIC") {
      select_solver(line);
      begin_static();
    } else if (card_ == "*HEATTRANSFER") {
      begin_heat_transfer(line);
    } else if (card_ == "*FREQUENCY") {
      begin_frequency(line);
    } else if (card_ == "*BUCKLE") {
      begin_buckle(line);
    } else if (card_ == "*DYNAMIC") {
      begin_dynamic(line);
    } else if (card_ == "*MODALDYNAMIC") {
      begin_modal_dynamic(line);
    } else if (card_ == "*STEADYSTATEDYNAMICS") {
      begin_steady_state_dynamics(line);
    } else if (card_ == "*COMPLEXFREQUENCY") {
      begin_complex_frequency(line);
    } else if (card_ == "*DAMPING") {
      begin_damping(line);
    } else if (card_ == "*MODALDAMPING") {
      begin_modal_damping(line);
    } else if (card_ == "*BASEMOTION") {
      begin_base_motion(line);
    } else if (card_ == "*SUBSTRUCTUREGENERATE") {
      begin_substructure_generate(line);
    } else if (card_ == "*RETAINEDNODALDOFS") {
      // Data lines carry "nset_or_node, dof1, dof2" — parsed by retained_nodal_dofs_data.
    } else if (card_ == "*SUBSTRUCTUREMATRIXOUTPUT") {
      begin_substructure_matrix_output(line);
    } else if (card_ == "*MATRIXASSEMBLE") {
      // *MATRIX ASSEMBLE requests export of the assembled GLOBAL K/M (STIFFNESS=/MASS=).
      // The reduced-substructure path already exposes the assembled operators; the flags
      // are accepted and recorded on the substructure output request. Data lines (none).
      begin_matrix_assemble(line);
    } else if (card_ == "*SENSITIVITY") {
      model_.procedure = Procedure::Sensitivity;
    } else if (card_ == "*DESIGNVARIABLES") {
      begin_design_variables(line);
    } else if (card_ == "*DESIGNRESPONSE" || card_ == "*OBJECTIVE") {
      // NAME= labels the response; its type (STRAIN ENERGY / nodal displacement) and
      // scope come on the data line. *OBJECTIVE is the Abaqus spelling of the same.
      cur_design_response_ = param("NAME");
    } else if (card_ == "*CONDUCTIVITY") {
      begin_conductivity(line);
    } else if (card_ == "*SPECIFICHEAT") {
      begin_specific_heat(line);
    } else if (card_ == "*EXPANSION") {
      begin_expansion(line);
    } else if (card_ == "*COUPLEDTEMPERATURE-DISPLACEMENT" ||
               card_ == "*COUPLEDTEMPERATUREDISPLACEMENT") {
      begin_coupled(line);
    } else if (card_ == "*PHYSICALCONSTANTS") {
      begin_physical_constants(line);
    } else if (card_ == "*INITIALCONDITIONS") {
      // TYPE=TEMPERATURE data lines set the transient initial field; other TYPEs
      // (STRESS/PLASTIC STRAIN/...) are accepted and ignored in this slice.
    } else if (card_ == "*FILM") {
      // Data lines carry elem,F<face>,sink,h. No OP reset in this slice.
    } else if (card_ == "*RADIATE") {
      // Data lines carry elem,R<face>,ambient,emissivity.
    } else if (card_ == "*CONTROLS") {
      // Parameters carried on the card; data lines fill the tolerances/limits.
    } else if (card_ == "*TIMEPOINTS") {
      // Data lines list the time points.
    } else if (card_ == "*MATERIAL") {
      cur_material_ = param("NAME");
      model_.materials[cur_material_].name = cur_material_;
    } else if (card_ == "*PLASTIC") {
      begin_plastic(line);
    } else if (card_ == "*CYCLICHARDENING") {
      begin_cyclic_hardening(line);
    } else if (card_ == "*HYPERELASTIC") {
      begin_hyperelastic(line);
    } else if (card_ == "*USERMATERIAL") {
      begin_user_material(line);
    } else if (card_ == "*DEPVAR") {
      begin_depvar(line);
    } else if (card_ == "*RATEDEPENDENT") {
      begin_ratedependent(line);
    } else if (card_ == "*AMPLITUDE") {
      begin_amplitude(line);
    } else if (card_ == "*SURFACE") {
      begin_surface(line);
    } else if (card_ == "*SURFACEINTERACTION") {
      begin_surface_interaction(line);
    } else if (card_ == "*SURFACEBEHAVIOR") {
      begin_surface_behavior(line);
    } else if (card_ == "*FRICTION") {
      begin_friction(line);
    } else if (card_ == "*GAPCONDUCTANCE") {
      begin_gap_conductance(line);
    } else if (card_ == "*GAPHEATGENERATION") {
      begin_gap_heat_generation(line);
    } else if (card_ == "*CONTACTPAIR") {
      begin_contact_pair(line);
    } else if (card_ == "*CLEARANCE") {
      begin_clearance(line);
    } else if (card_ == "*SOLIDSECTION") {
      model_.sections.push_back(SolidSection{param("ELSET"), param("MATERIAL")});
    } else if (card_ == "*CLOAD") {
      if (apply_op_reset(op_seen_cload_, line)) model_.cloads.clear();
    } else if (card_ == "*DLOAD") {
      // OP on *DLOAD resets the whole distributed-load set (pressures + body loads).
      if (apply_op_reset(op_seen_dload_, line)) {
        model_.dloads.clear();
        model_.body_loads.clear();
      }
    } else if (card_ == "*BOUNDARY") {
      if (apply_op_reset(op_seen_boundary_, line)) {
        model_.spcs.clear();
        model_.temp_bcs.clear();
      }
    } else if (card_ == "*CFLUX") {
      if (apply_op_reset(op_seen_cflux_, line)) model_.cfluxes.clear();
    } else if (card_ == "*DFLUX") {
      if (apply_op_reset(op_seen_dflux_, line)) model_.dfluxes.clear();
    } else if (card_ == "*TEMPERATURE") {
      // Data lines carry node,value prescribed temperatures (treated as thermal BCs).
    } else if (card_ == "*SPRING" || card_ == "*DASHPOT") {
      // ELSET= names the connector elements this property attaches to. The
      // stiffness/damping (and, for SPRING1/2, the DOFs) come on the data lines.
      if (param("ELSET").empty())
        throw ParseError(line, std::string(card_ == "*SPRING" ? "*SPRING" : "*DASHPOT") +
                                   " needs ELSET=");
      spring_dofs_read_ = false;
      spring_dof1_ = spring_dof2_ = 0;
    } else if (card_ == "*MASS") {
      if (param("ELSET").empty()) throw ParseError(line, "*MASS needs ELSET=");
    } else if (card_ == "*EQUATION") {
      eq_terms_remaining_ = -1;  // first data line gives the term count
    } else if (card_ == "*MPC") {
      // Type + nodes come on the data line (mpc_data).
    } else if (card_ == "*RIGIDBODY") {
      begin_rigid_body(line);
    } else if (card_ == "*COUPLING") {
      begin_coupling(line);
    } else if (card_ == "*DISTRIBUTINGCOUPLING") {
      begin_distributing_coupling(line);
    } else if (card_ == "*KINEMATIC" || card_ == "*DISTRIBUTING") {
      // Modifier under *COUPLING: sets the coupling kind; DOF ranges follow as data.
      set_coupling_kind(line);
    } else if (card_ == "*TIE") {
      begin_tie(line);
    } else if (card_ == "*MODELCHANGE") {
      begin_model_change(line);
    } else if (card_ == "*CHANGEMATERIAL") {
      begin_change_material(line);
    } else if (card_ == "*CHANGEPLASTIC") {
      // Redefines plastic hardening data of the current material. Plasticity is
      // not implemented yet (workstream 4), so the data lines are accepted and
      // ignored — nothing to change. See tasks.md 2.3 ([~], plasticity-dependent).
    } else if (card_ == "*CHANGESOLIDSECTION") {
      // Re-bind a material to an element set within the step. Appending a section
      // makes it win over earlier ones per element (element_elastic() is
      // last-writer). MATERIAL and ELSET are required.
      if (param("MATERIAL").empty() || param("ELSET").empty())
        throw ParseError(line, "*CHANGE SOLID SECTION needs MATERIAL= and ELSET=");
      model_.sections.push_back(SolidSection{param("ELSET"), param("MATERIAL")});
    } else if (std::find(ignored.begin(), ignored.end(), card_) == ignored.end() &&
               !is_data_card()) {
      reject_card(fields.front(), line);
    }
  }

  // Reject a card the parser does not implement. Capabilities that are
  // deliberately DEFERRED (documented in tasks.md) get a clear, actionable
  // message that names the deferral, so a deck using them fails loudly and
  // specifically rather than being silently mis-solved. Any other keyword is a
  // generic "unsupported card". Either way this THROWS (never returns) — a
  // deferred card must never be silently ignored, which would yield a wrong
  // solve. (spec: input-deck-parsing — graceful handling of unknown cards.)
  [[noreturn]] void reject_card(const std::string& raw, int line) const {
    // Normalized keyword -> the workstream/reason it is deferred (Phase-2
    // material/section families that need extra kinematics or a transient driver).
    static const std::unordered_map<std::string, std::string> deferred = {
        {"*HYPERFOAM", "hyperelastic foam (Ogden) — deferred (tasks 4.3)"},
        {"*CREEP", "creep — deferred, needs a transient driver (tasks 4.4)"},
        {"*VISCO", "viscoelasticity — deferred, needs a transient driver (tasks 4.4)"},
        {"*VALUESATINFINITY",
         "viscoelastic long-term moduli — deferred with *VISCO (tasks 4.4)"},
        {"*DEFORMATIONPLASTICITY",
         "deformation (total-strain) plasticity — deferred (tasks 4.5)"},
        {"*MOHRCOULOMB", "Mohr-Coulomb plasticity — deferred (tasks 4.5)"},
        {"*MOHRCOULOMBHARDENING",
         "Mohr-Coulomb hardening — deferred with *MOHR COULOMB (tasks 4.5)"},
        {"*DAMAGEINITIATION", "damage initiation — deferred (tasks 4.5)"},
        {"*SHELLSECTION",
         "shell sections — deferred, needs shell kinematics (tasks 3.4)"},
        {"*BEAMSECTION",
         "beam sections — deferred, needs beam kinematics (tasks 3.4)"},
        {"*MEMBRANESECTION",
         "membrane sections — deferred, needs shell kinematics (tasks 3.4)"},
        {"*VIEWFACTOR",
         "the *VIEWFACTOR read/write-to-file control is not implemented; cavity "
         "radiation itself (*RADIATE ...,CR) IS supported — the view factors are "
         "computed internally, so drop *VIEWFACTOR (tasks 3.4)"},
        // Phase-4 procedure cards blocked on a missing enabler (documented in
        // openspec/BACKLOG.md Phase-4). The eigensolution/dynamics engines they
        // would consume are in place; each names the enabler it waits on.
        {"*CYCLICSYMMETRYMODEL",
         "cyclic-symmetry sector eigenproblem — deferred; needs complex "
         "cyclic-symmetry constraints and the complex eigensolve (Phase-4 tasks 6.1)"},
        {"*SELECTCYCLICSYMMETRYMODES",
         "nodal-diameter mode selection — deferred with *CYCLIC SYMMETRY MODEL "
         "(Phase-4 tasks 6.2)"},
        {"*GREEN",
         "Green-function step — deferred; builds on the shipped modal engine but "
         "needs the unit-excitation response basis (Phase-4 tasks 4.4)"},
        // Design-optimization layers on top of the shipped sensitivity core, deferred
        // in this slice (openspec add-design-optimization tasks 1.6 / 2.5 / 4.x).
        {"*FILTER",
         "sensitivity filtering — deferred; the raw adjoint gradient (*SENSITIVITY) "
         "IS supported (add-design-optimization tasks 2.5)"},
        {"*FEASIBLEDIRECTION",
         "feasible-direction optimization loop — deferred; consumes the shipped "
         "sensitivity gradient (add-design-optimization tasks 4.1)"},
        {"*ROBUSTDESIGN",
         "robust design — deferred; needs a random-field decomposition "
         "(add-design-optimization tasks 4.2)"},
        {"*RANDOMFIELD",
         "random field — deferred with *ROBUST DESIGN (add-design-optimization tasks 4.2)"},
        {"*CORRELATIONLENGTH",
         "correlation length — deferred with *ROBUST DESIGN (add-design-optimization "
         "tasks 4.2)"},
        {"*GEOMETRICCONSTRAINT",
         "geometric constraint — deferred; part of the optimization loop "
         "(add-design-optimization tasks 4.1)"},
    };
    const auto it = deferred.find(card_);
    if (it != deferred.end())
      throw ParseError(line, "card '" + raw + "' is a recognized capability that is "
                                             "not yet implemented: " +
                                 it->second + ". Remove it or use an implemented "
                                 "material/section; it cannot be silently ignored.");
    // Phase-3 contact + thermal-contact cards (*CONTACT PAIR, *SURFACE INTERACTION,
    // *SURFACE BEHAVIOR, *FRICTION, *GAP CONDUCTANCE, *GAP HEAT GENERATION, *CLEARANCE,
    // *CONTACT FILE/PRINT/OUTPUT) are all handled above; anything reaching here is a
    // genuinely unimplemented keyword.
    throw ParseError(line, "unsupported card '" + raw + "'");
  }

  // Map SOLVER= on *STATIC onto a RequestedSolver (spec 9.2/9.3). Direct family
  // (SPOOLES/PARDISO/PASTIX) -> Direct; ITERATIVE*/CG -> CG; empty -> default
  // Direct. An unrecognized/unavailable name throws with the requested name,
  // mirroring CalculiX's 'solver not available' behavior.
  void select_solver(int line) {
    std::string name = param("SOLVER");
    name.erase(std::remove_if(name.begin(), name.end(),
                              [](unsigned char c) { return std::isspace(c) != 0; }),
               name.end());
    if (name.empty()) {
      model_.solver = RequestedSolver::Auto;  // size-based choice at solve time
    } else if (name == "SPOOLES" || name == "PARDISO" || name == "PASTIX" ||
               name == "DIRECT") {
      model_.solver = RequestedSolver::Direct;
    } else if (name == "CG" || name.rfind("ITERATIVE", 0) == 0) {
      model_.solver = RequestedSolver::CG;
    } else {
      throw ParseError(line, "solver not available: SOLVER=" + param("SOLVER"));
    }
  }

  // *STEP: open a new analysis step (spec: multi-step analysis). Reset the per-step
  // OP=NEW guards so an OP=NEW load/BC card resets once per step (CalculiX's
  // once-per-type-per-step rule), and clear the heat-step flag so each step's
  // procedure card re-selects the routing. All accumulated loads/BCs/state persist
  // (they carry forward from the prior step); OP=NEW inside the step clears them.
  void begin_step() {
    in_step_ = true;
    op_seen_cload_ = op_seen_dload_ = op_seen_boundary_ = false;
    op_seen_cflux_ = op_seen_dflux_ = false;
    heat_step_ = false;
  }

  // *END STEP: snapshot the fully-accumulated model as this step's Model.
  void end_step() {
    snapshot_step();
    in_step_ = false;
  }

  // Capture the current accumulated model_ as one step. Called at *END STEP (and, for
  // an unterminated final block, at end of parse). The snapshot is a plain copy — the
  // shared sets are flushed into its mesh later by flush_sets_into.
  void snapshot_step() { steps_.push_back(model_); }

  // *STATIC card: DIRECT parameter fixes the increment (no auto resizing).
  void begin_static() {
    if (params_.count("DIRECT") != 0) model_.increment.direct = true;
  }

  // *HEAT TRANSFER[, STEADY STATE]: select the thermal procedure. STEADY STATE
  // solves Kt T = q; without it the step is transient (backward-Euler integration of
  // the capacitance term over the step period). (spec: heat-transfer-analysis —
  // steady + transient conduction.) The (unused) `line` keeps the card-handler
  // signature uniform.
  void begin_heat_transfer(int /*line*/) {
    heat_step_ = true;
    model_.procedure = (params_.count("STEADYSTATE") != 0)
                           ? Procedure::HeatTransferSteady
                           : Procedure::HeatTransferTransient;
  }

  // *FREQUENCY: natural-frequency extraction (spec: modal-and-buckling — *FREQUENCY).
  // Selects the eigen procedure; the data line gives the number of eigenvalues
  // requested (first field). The optional second field (a frequency lower bound in
  // CalculiX) is parsed and ignored in this slice — the dense engine returns the
  // lowest N modes regardless. STORAGE/other params are accepted and ignored. `line`
  // keeps the handler signature uniform.
  // A *FREQUENCY card INSIDE a *SUBSTRUCTURE GENERATE step requests the number of
  // fixed-interface normal modes retained by the Craig-Bampton reduction (not a
  // standalone frequency analysis) — keep the Substructure procedure in that case.
  void begin_frequency(int /*line*/) {
    if (model_.procedure != Procedure::Substructure)
      model_.procedure = Procedure::Frequency;
  }

  // *FREQUENCY data line: "N[, ...]" — number of eigenvalues requested. Defaults to 1
  // when the line is blank (CalculiX requires N; we default defensively). Inside a
  // substructure step the count is the number of fixed-interface modes.
  void frequency_data(const std::vector<std::string>& f, int line) {
    const int n = (f.empty() || f[0].empty())
                      ? 1
                      : std::max(1, static_cast<int>(to_double(f[0], line)));
    if (model_.procedure == Procedure::Substructure)
      model_.substructure_modes = n;
    else
      model_.num_eigenvalues = n;
  }

  // *BUCKLE: linear-buckling analysis (spec: input-deck-parsing — *BUCKLE). Selects the
  // buckling procedure; the data line gives the number of buckling factors requested
  // (first field), as *FREQUENCY records num_eigenvalues. The trailing ARPACK-style
  // accuracy / subspace / iteration fields (e.g. the `1.e-2` in `10, 1.e-2`) are
  // recognized but do not affect the model. `line` keeps the handler signature uniform.
  void begin_buckle(int /*line*/) { model_.procedure = Procedure::Buckling; }

  // *BUCKLE data line: "N[, accuracy, ncv, maxiter]" — number of buckling modes (first
  // field). Defaults to 1 when blank. The trailing ARPACK-style fields are accepted and
  // ignored (the dense engine returns the lowest N positive factors regardless).
  void buckle_data(const std::vector<std::string>& f, int line) {
    model_.num_buckling_modes =
        (f.empty() || f[0].empty())
            ? 1
            : std::max(1, static_cast<int>(to_double(f[0], line)));
  }

  // *DESIGN VARIABLES, TYPE=COORDINATE: declare the design field for a *SENSITIVITY
  // step (spec: design-optimization). Only TYPE=COORDINATE (shape) variables are
  // implemented; TYPE=ORIENTATION (sensi_orien) is deferred and rejected here so a deck
  // using it fails loudly rather than being silently mis-solved. The referenced nodes
  // come on the data line(s) as an nset name or node ids.
  void begin_design_variables(int line) {
    const std::string type = param("TYPE");
    if (!type.empty() && type != "COORDINATE")
      throw ParseError(line,
                       "*DESIGN VARIABLES, TYPE=" + type +
                           " is not implemented (only TYPE=COORDINATE shape variables "
                           "are supported this slice); TYPE=ORIENTATION is deferred.");
  }

  // *DESIGN VARIABLES data: an nset name or node id per field. Each referenced node
  // contributes three coordinate design variables (x, y, z), in node-then-component
  // order. Duplicate nodes are de-duplicated so a node listed twice is one variable set.
  void design_variables_data(const std::vector<std::string>& f, int line) {
    for (const std::string& tok : f) {
      if (tok.empty()) continue;
      for (const Index nid : node_refs(tok, line)) {
        bool seen = false;
        for (const DesignVariable& dv : model_.design_variables)
          if (dv.node_id == nid) { seen = true; break; }
        if (seen) continue;
        for (int c = 1; c <= 3; ++c)
          model_.design_variables.push_back(DesignVariable{nid, c});
      }
    }
  }

  // *DESIGN RESPONSE / *OBJECTIVE data: "<type>[, scope]". Supported types:
  //   STRAIN ENERGY / ENERGY / COMPLIANCE  -> the self-adjoint compliance g = fᵀu
  //     (the scope elset, e.g. EALL, is accepted; the response is the whole-model
  //     external work in this slice).
  //   ALL-DISP / U / DISPLACEMENT, <node>, <comp> -> a single nodal displacement DOF.
  // Mass / stress / eigenfrequency responses are out of this slice and rejected.
  void design_response_data(const std::vector<std::string>& f, int line) {
    if (f.empty() || f[0].empty())
      throw ParseError(line, "*DESIGN RESPONSE needs a response type");
    DesignResponse resp;
    resp.name = cur_design_response_;
    const std::string type = normalize(f[0]);  // uppercase, spaces stripped
    if (type == "STRAINENERGY" || type == "ENERGY" || type == "COMPLIANCE") {
      resp.kind = DesignResponse::Kind::Compliance;
    } else if (type == "ALL-DISP" || type == "ALLDISP" || type == "U" ||
               type == "DISPLACEMENT") {
      resp.kind = DesignResponse::Kind::Displacement;
      if (f.size() < 3)
        throw ParseError(line,
                         "*DESIGN RESPONSE displacement needs '<type>, node, comp'");
      resp.node_id = to_index(f[1], line);
      resp.comp = static_cast<int>(to_index(f[2], line));
      if (resp.comp < 1 || resp.comp > 3)
        throw ParseError(line, "*DESIGN RESPONSE displacement comp must be 1..3");
    } else {
      throw ParseError(line, "*DESIGN RESPONSE type '" + type +
                                 "' is not implemented (supported: STRAIN ENERGY, "
                                 "DISPLACEMENT); mass/stress/frequency are deferred.");
    }
    model_.design_responses.push_back(std::move(resp));
  }

  // *DYNAMIC[, DIRECT][, ALPHA=][, NLGEOM]: direct time integration of the equations of
  // motion by the implicit HHT-α scheme (spec: dynamic-analysis — direct dynamics). The
  // data line is "dt, t_end" (shared with *MODAL DYNAMIC). ALPHA= sets the HHT numerical
  // damping (α ∈ [-1/3, 0]); NLGEOM (or a nonlinear material / contact in the model)
  // routes the per-step Newton path. EXPLICIT is accepted but the implicit scheme is used
  // (the stable-implicit default); the parser records it without changing the scheme.
  void begin_dynamic(int line) {
    model_.procedure = Procedure::Dynamic;
    if (params_.count("ALPHA") != 0)
      model_.dynamic_alpha = to_double(param("ALPHA"), line);
    if (params_.count("NLGEOM") != 0) model_.dynamic_nonlinear = true;
    if (params_.count("DIRECT") != 0) model_.increment.direct = true;
  }

  // *MODAL DYNAMIC: transient response by modal superposition over a preceding
  // *FREQUENCY basis (spec: dynamic-analysis — modal dynamic). The data line is
  // "dt, t_end"; both default to a single unit step when blank. SOLVER= is accepted and
  // ignored (modal superposition has no linear factorization). `line` keeps the handler
  // signature uniform.
  void begin_modal_dynamic(int /*line*/) {
    model_.procedure = Procedure::ModalDynamic;
  }
  void modal_dynamic_data(const std::vector<std::string>& f, int line) {
    if (!f.empty() && !f[0].empty()) model_.dynamic_dt = to_double(f[0], line);
    if (f.size() > 1 && !f[1].empty()) model_.dynamic_t_end = to_double(f[1], line);
  }

  // *STEADY STATE DYNAMICS: harmonic response over a frequency sweep by modal
  // superposition (spec: dynamic-analysis — steady-state dynamics). The data line is
  // "f_lo, f_hi, num_points[, bias]"; num_points defaults to 20. `line` uniform.
  void begin_steady_state_dynamics(int /*line*/) {
    model_.procedure = Procedure::SteadyStateDynamics;
  }
  void steady_state_data(const std::vector<std::string>& f, int line) {
    if (!f.empty() && !f[0].empty()) model_.steady_f_lo = to_double(f[0], line);
    if (f.size() > 1 && !f[1].empty()) model_.steady_f_hi = to_double(f[1], line);
    if (f.size() > 2 && !f[2].empty())
      model_.steady_num_points = static_cast<int>(to_double(f[2], line));
    if (model_.steady_num_points <= 0) model_.steady_num_points = 20;
  }

  // *COMPLEX FREQUENCY: damped complex modes over a preceding *FREQUENCY basis (spec:
  // modal-and-buckling-analysis — complex frequency). This slice ships the option-(B)
  // proportional-damping reduction only. CORIOLIS (the gyroscopic skew operator + rotor
  // rotation body load) and FLUTTER (a complex applied force) are a DIFFERENT eigenproblem
  // that needs input this deck does not carry — rather than silently mis-solve them with
  // the wrong physics, they are rejected here with an explicit "not yet implemented"
  // (grounded: CalculiX complexfrequencys.f:86-92 errors unless CORIOLIS or FLUTTER is
  // present, and CORIOLIS consumes a rotation body load via xbody). The keyword-less
  // (or PROPORTIONAL) card drives the shipped path. The data line = number of complex
  // modes. `line` uniform.
  void begin_complex_frequency(int line) {
    model_.procedure = Procedure::ComplexFrequency;
    if (params_.count("CORIOLIS") != 0) {
      model_.complex_freq_type = Model::ComplexFrequencyType::Coriolis;
      throw ParseError(line,
                       "*COMPLEX FREQUENCY, CORIOLIS is not yet implemented: the "
                       "gyroscopic (skew) Coriolis operator and its rotor rotation body "
                       "load are not available; only proportional-damping complex modes "
                       "are supported");
    }
    if (params_.count("FLUTTER") != 0) {
      model_.complex_freq_type = Model::ComplexFrequencyType::Flutter;
      throw ParseError(line,
                       "*COMPLEX FREQUENCY, FLUTTER is not yet implemented: the complex "
                       "applied-force (flutter) path is not available; only "
                       "proportional-damping complex modes are supported");
    }
    model_.complex_freq_type = Model::ComplexFrequencyType::Proportional;
  }
  // *COMPLEX FREQUENCY data line: "N[, ...]" — number of complex modes requested.
  // Defaults to 1 when blank (defensive; CalculiX requires N).
  void complex_frequency_data(const std::vector<std::string>& f, int line) {
    const int n = (f.empty() || f[0].empty())
                      ? 1
                      : std::max(1, static_cast<int>(to_double(f[0], line)));
    model_.num_complex_modes = n;
  }

  // *DAMPING, ALPHA=, BETA=: Rayleigh (proportional) damping C = alpha M + beta K
  // (spec: dynamic-analysis — Rayleigh damping). Coefficients on the card params; a data
  // line (if present) is ignored.
  void begin_damping(int line) {
    model_.rayleigh.active = true;
    if (params_.count("ALPHA") != 0)
      model_.rayleigh.alpha = to_double(param("ALPHA"), line);
    if (params_.count("BETA") != 0)
      model_.rayleigh.beta = to_double(param("BETA"), line);
  }

  // *MODAL DAMPING: explicit per-mode damping ratios (spec: dynamic-analysis — modal
  // damping). Data lines are "mode_lo, mode_hi, zeta" ranges (CalculiX form); each row
  // sets zeta on modes [mode_lo, mode_hi] (1-based). `line` uniform.
  void begin_modal_damping(int /*line*/) {}
  void modal_damping_data(const std::vector<std::string>& f, int line) {
    if (f.size() < 3) return;
    const int lo = static_cast<int>(to_double(f[0], line));
    const int hi = static_cast<int>(to_double(f[1], line));
    const Real zeta = to_double(f[2], line);
    if (lo < 1 || hi < lo) return;
    if (static_cast<int>(model_.modal_damping.size()) < hi)
      model_.modal_damping.resize(static_cast<std::size_t>(hi), 0.0);
    for (int mdx = lo; mdx <= hi; ++mdx)
      model_.modal_damping[static_cast<std::size_t>(mdx - 1)] = zeta;
  }

  // *BASE MOTION, DOF=, AMPLITUDE=: prescribed support excitation for a modal procedure
  // (spec: dynamic-analysis — base motion). DOF= is the excitation direction (1..3); the
  // data line (or MAGNITUDE) gives the amplitude. `line` uniform.
  void begin_base_motion(int line) {
    model_.base_motion.active = true;
    if (params_.count("DOF") != 0)
      model_.base_motion.dof = static_cast<int>(to_double(param("DOF"), line));
  }
  void base_motion_data(const std::vector<std::string>& f, int line) {
    if (!f.empty() && !f[0].empty())
      model_.base_motion.magnitude = to_double(f[0], line);
  }

  // *SUBSTRUCTURE GENERATE: select the substructure (superelement) procedure (spec:
  // substructure-generation — tasks 5.1-5.3). The retained (master) DOFs come from a
  // following *RETAINED NODAL DOFS card; the reduction is Guyan (static) unless a mass
  // matrix and fixed-interface modes are requested (Craig-Bampton). No data line.
  void begin_substructure_generate(int /*line*/) {
    model_.procedure = Procedure::Substructure;
  }

  // *RETAINED NODAL DOFS data: "node_or_nset, first_dof, last_dof" — every DOF in
  // [first, last] (clamped to 1..3 translational) of each node in the set is retained,
  // in declaration order (SORTED=NO). Duplicates are dropped so a node listed twice is
  // retained once. (spec: substructure-generation — retained nodal DOFs definition.)
  void retained_nodal_dofs_data(const std::vector<std::string>& f, int line) {
    if (f.empty() || f[0].empty())
      throw ParseError(line, "*RETAINED NODAL DOFS needs node/nset, first_dof[, last_dof]");
    const int d1 = f.size() > 1 && !f[1].empty()
                       ? static_cast<int>(to_index(f[1], line))
                       : 1;
    const int d2 = f.size() > 2 && !f[2].empty()
                       ? static_cast<int>(to_index(f[2], line))
                       : d1;
    for (const Index nd : node_refs(f[0], line))
      for (int c = std::max(1, d1); c <= std::min(d2, kDofsPerNode); ++c) {
        const bool dup = std::any_of(
            model_.retained_dofs.begin(), model_.retained_dofs.end(),
            [&](const Model::RetainedDof& r) { return r.node_id == nd && r.comp == c; });
        if (!dup) model_.retained_dofs.push_back(Model::RetainedDof{nd, c});
      }
  }

  // *SUBSTRUCTURE MATRIX OUTPUT, STIFFNESS=YES|NO, MASS=YES|NO: select which reduced
  // matrices the superelement exports (spec: substructure-generation — reduced matrix
  // export). MASS=YES (or a nonzero fixed-interface mode count) triggers Craig-Bampton;
  // stiffness defaults to YES. OUTPUT FILE / FILE NAME params are accepted and ignored
  // (the export path is chosen by the driver / bindings). No data line.
  void begin_substructure_matrix_output(int /*line*/) {
    const auto yes = [&](const char* key, bool dflt) {
      if (params_.count(key) == 0) return dflt;
      const std::string& v = params_.at(key);
      return v.empty() || v == "YES";
    };
    model_.substructure_stiffness = yes("STIFFNESS", true);
    model_.substructure_mass = yes("MASS", false);
  }

  // *MATRIX ASSEMBLE, STIFFNESS=YES, MASS=YES: request export of the assembled global
  // operators. In this slice we record the mass request (so a *MATRIX ASSEMBLE, MASS=YES
  // forms the reduced mass too); the assembled global K/M are reachable via the
  // fem::assemble_* API and the substructure export. No data line.
  void begin_matrix_assemble(int /*line*/) {
    if (params_.count("MASS") != 0) {
      const std::string& v = params_.at("MASS");
      if (v.empty() || v == "YES") model_.substructure_mass = true;
    }
  }

  // *CONDUCTIVITY: isotropic scalar conductivity on the current material. The data
  // line is "k[, temperature]"; only the first (isotropic) value is used.
  void begin_conductivity(int line) {
    if (cur_material_.empty())
      throw ParseError(line, "*CONDUCTIVITY without *MATERIAL");
    ensure_thermal();
  }

  // *SPECIFIC HEAT: scalar specific heat on the current material (transient only).
  void begin_specific_heat(int line) {
    if (cur_material_.empty())
      throw ParseError(line, "*SPECIFIC HEAT without *MATERIAL");
    ensure_thermal();
  }

  Thermal& ensure_thermal() {
    auto& mat = model_.materials[cur_material_];
    if (!mat.thermal) mat.thermal = Thermal{};
    return *mat.thermal;
  }

  // *EXPANSION[, ZERO=Tref]: isotropic thermal-expansion coefficient on the current
  // material. ZERO= gives the reference (stress-free) temperature Tref (default 0).
  // The data line is "alpha[, temperature]"; only the first isotropic value is used.
  // (spec: heat-transfer — thermal expansion; material-models — *EXPANSION.)
  void begin_expansion(int line) {
    if (cur_material_.empty())
      throw ParseError(line, "*EXPANSION without *MATERIAL");
    Expansion& ex = ensure_expansion();
    if (params_.count("ZERO") != 0) ex.t_ref = to_double(param("ZERO"), line);
  }

  Expansion& ensure_expansion() {
    auto& mat = model_.materials[cur_material_];
    if (!mat.expansion) mat.expansion = Expansion{};
    return *mat.expansion;
  }

  void expansion_data(const std::vector<std::string>& f, int line) {
    if (cur_material_.empty())
      throw ParseError(line, "*EXPANSION without *MATERIAL");
    append_property_row(ensure_expansion().alpha, f, line, "*EXPANSION");
  }

  // *COUPLED TEMPERATURE-DISPLACEMENT[, SOLUTIONS=MONOLITHIC|STAGGERED]: select the
  // coupled thermomechanical procedure. A heat step for the thermal cards (temperature
  // BCs, fluxes, film) whose solve also computes the displacement/thermal-stress field.
  // SOLUTIONS= picks the scheme (default STAGGERED — its single-pass degeneracy is the
  // historical one-way solve; MONOLITHIC assembles the 4-DOF/node system). Both agree
  // for a pure thermal-stress problem. (spec: heat-transfer — coupled.) `line` keeps the
  // handler signature uniform.
  void begin_coupled(int /*line*/) {
    heat_step_ = true;  // thermal cards (temperature BCs, flux, film) route to thermal
    model_.procedure = Procedure::Coupled;
    const std::string sol = param("SOLUTIONS");
    if (sol == "MONOLITHIC")
      model_.coupled_scheme = CoupledScheme::Monolithic;
    else if (sol == "STAGGERED")
      model_.coupled_scheme = CoupledScheme::Staggered;
  }

  // Append one (value[, temperature]) row to a temperature-dependent property table.
  // A row with only a value (or a bare 0 temperature on the first row) is the constant
  // case; multiple rows form a piecewise-linear k(T)/c(T)/alpha(T). Temperatures must
  // be strictly increasing (CalculiX requirement); an omitted temperature defaults 0.
  void append_property_row(PropertyTable& table, const std::vector<std::string>& f,
                           int line, const char* card) {
    if (f.empty() || f[0].empty())
      throw ParseError(line, std::string(card) + " needs a value");
    const Real value = to_double(f[0], line);
    const Real temperature = (f.size() > 1 && !f[1].empty()) ? to_double(f[1], line) : 0.0;
    if (!table.temp.empty() && temperature <= table.temp.back())
      throw ParseError(line, std::string(card) +
                                 " temperatures must strictly increase");
    table.value.push_back(value);
    table.temp.push_back(temperature);
  }

  void conductivity_data(const std::vector<std::string>& f, int line) {
    if (cur_material_.empty())
      throw ParseError(line, "*CONDUCTIVITY without *MATERIAL");
    append_property_row(ensure_thermal().conductivity, f, line, "*CONDUCTIVITY");
  }

  void specific_heat_data(const std::vector<std::string>& f, int line) {
    if (cur_material_.empty())
      throw ParseError(line, "*SPECIFIC HEAT without *MATERIAL");
    append_property_row(ensure_thermal().specific_heat, f, line, "*SPECIFIC HEAT");
  }

  // *PHYSICAL CONSTANTS: ABSOLUTE ZERO / STEFAN BOLTZMANN carried as card
  // parameters (radiation converts to absolute temperature and needs sigma).
  void begin_physical_constants(int line) {
    if (params_.count("ABSOLUTEZERO") != 0)
      model_.physical.absolute_zero = to_double(param("ABSOLUTEZERO"), line);
    if (params_.count("STEFANBOLTZMANN") != 0)
      model_.physical.sigma = to_double(param("STEFANBOLTZMANN"), line);
  }

  // *INITIAL CONDITIONS, TYPE=TEMPERATURE data: "node_or_nset, value". Sets the
  // transient initial temperature field (a bare "Nall" set covering every node acts
  // as the uniform default). Non-TEMPERATURE types are ignored.
  void initial_conditions_data(const std::vector<std::string>& f, int line) {
    if (param("TYPE") != "TEMPERATURE") return;
    if (f.size() < 2) throw ParseError(line, "*INITIAL CONDITIONS needs node,value");
    const Real val = to_double(f[1], line);
    const std::vector<Index> nodes = node_refs(f[0], line);
    // A set spanning all nodes doubles as the uniform default; still record the
    // per-node values so partial sets win over it.
    if (nodes.size() == model_.mesh.num_nodes()) model_.initial_temperature = val;
    for (const Index nd : nodes) model_.initial_temperature_by_node[nd] = val;
  }

  // Parse a face label like "F1"/"R3"/"S2" -> face number, validating the prefix.
  int face_label(const std::string& raw, char prefix, const char* card, int line) {
    const std::string label = upper(raw);
    if (label.empty() || label[0] != prefix)
      throw ParseError(line, std::string(card) + " face must be " + prefix + "<n>");
    return static_cast<int>(to_index(label.substr(1), line));
  }

  // *FILM data: "elem_or_elset, F<face>, sink_temp, film_coefficient". Applies a
  // convective boundary q = h (T - T_sink) on the face. (Forced-convection network
  // film labels FC/F<face>N are out of scope for this solid-only slice.)
  void film_data(const std::vector<std::string>& f, int line) {
    if (f.size() < 4) throw ParseError(line, "*FILM needs elem,F<face>,sink,h");
    const int face = face_label(f[1], 'F', "*FILM", line);
    const Real sink = to_double(f[2], line);
    const Real h = to_double(f[3], line);
    const std::string amp = param("AMPLITUDE");
    for (const Index eid : elem_refs(f[0], line))
      model_.films.push_back(Film{eid, face, sink, h, amp});
  }

  // *RADIATE data: "elem_or_elset, R<face>[CR], ambient_temp, emissivity". Two forms:
  //   R<face>       surface-to-ambient q = eps sigma (T^4 - T_amb^4) against a fixed
  //                 ambient temperature (field 3);
  //   R<face>CR     gray-body cavity patch — the face joins the radiating enclosure
  //                 whose net flux is computed from the view factors between all CR
  //                 faces (field 3, the environment temperature, is optional/ignored
  //                 in this closed-enclosure slice; field 4 is the emissivity).
  void radiate_data(const std::vector<std::string>& f, int line) {
    if (f.size() < 4) throw ParseError(line, "*RADIATE needs elem,R<face>,ambient,eps");
    std::string label = upper(f[1]);
    const bool cavity = label.size() >= 2 && label.compare(label.size() - 2, 2, "CR") == 0;
    if (cavity) label.erase(label.size() - 2);  // strip the CR suffix -> R<face>
    const int face = face_label(label, 'R', "*RADIATE", line);
    const Real amb = f[2].empty() ? 0.0 : to_double(f[2], line);
    const Real eps = to_double(f[3], line);
    const std::string amp = param("AMPLITUDE");
    for (const Index eid : elem_refs(f[0], line))
      model_.radiates.push_back(Radiate{eid, face, amb, eps, cavity, amp});
  }

  // *CFLUX data: "node_or_nset, dof, value". The dof (0 or 11) selects the
  // temperature DOF and is ignored; the flux is added at the node.
  void cflux_data(const std::vector<std::string>& f, int line) {
    if (f.size() < 3) throw ParseError(line, "*CFLUX needs node,dof,value");
    const Real val = to_double(f[2], line);
    const std::string amp = param("AMPLITUDE");
    for (const Index nd : node_refs(f[0], line))
      model_.cfluxes.push_back(Cflux{nd, val, amp});
  }

  // *DFLUX data: "elem_or_elset, S<face>, magnitude" (distributed surface flux). A
  // BF (body-heat-generation) label is not implemented in this slice.
  void dflux_data(const std::vector<std::string>& f, int line) {
    if (f.size() < 3) throw ParseError(line, "*DFLUX needs element,S<face>,magnitude");
    const std::string label = upper(f[1]);
    if (label == "BF")
      throw ParseError(line,
                       "*DFLUX BF (body heat generation) is not implemented yet; "
                       "only surface flux S<face> is supported");
    if (label.empty() || label[0] != 'S')
      throw ParseError(line, "unsupported *DFLUX type '" + f[1] +
                                 "' (supports surface flux S<face>)");
    const int face = static_cast<int>(to_index(label.substr(1), line));
    const Real flux = to_double(f[2], line);
    const std::string amp = param("AMPLITUDE");
    for (const Index eid : elem_refs(f[0], line))
      model_.dfluxes.push_back(Dflux{eid, face, flux, amp});
  }

  // *TEMPERATURE data: "node_or_nset, value". In a heat-transfer step it prescribes a
  // nodal temperature BC (an alternative to *BOUNDARY dof 11). In a MECHANICAL / coupled
  // step it defines the applied temperature field that drives thermal expansion
  // (eps_th = alpha (T - Tref)); the field lands in model.applied_temperature keyed by
  // node id. (spec: heat-transfer — thermal expansion coupling.)
  void temperature_data(const std::vector<std::string>& f, int line) {
    if (f.size() < 2) throw ParseError(line, "*TEMPERATURE needs node,value");
    const Real val = to_double(f[1], line);
    const std::string amp = param("AMPLITUDE");
    for (const Index nd : node_refs(f[0], line)) {
      if (heat_step_)
        model_.temp_bcs.push_back(TempBc{nd, val, amp});
      else
        model_.applied_temperature[nd] = val;
    }
  }

  // *STATIC data line: initial_inc, total_time, min_inc, max_inc (all optional,
  // CalculiX order). Empty/missing fields keep the current defaults.
  void static_data(const std::vector<std::string>& f, int line) {
    Incrementation& inc = model_.increment;
    if (f.size() > 0 && !f[0].empty()) inc.initial = to_double(f[0], line);
    if (f.size() > 1 && !f[1].empty()) inc.total = to_double(f[1], line);
    if (f.size() > 2 && !f[2].empty()) inc.min = to_double(f[2], line);
    if (f.size() > 3 && !f[3].empty()) inc.max = to_double(f[3], line);
    if (f.size() <= 3) inc.max = inc.total;  // default max = whole step
  }

  // *CONTROLS: PARAMETERS=FIELD gives convergence tolerances (first field is the
  // force-residual tolerance R_n^alpha); PARAMETERS=TIME INCREMENTATION gives the
  // iteration limit (first field I_0). Other rows are accepted and ignored. Absent
  // fields keep the documented defaults on NonlinearControls. (ref: controlss.f)
  void controls_data(const std::vector<std::string>& f, int line) {
    const std::string p = param("PARAMETERS");
    if (p == "FIELD") {
      if (!f.empty() && !f[0].empty())
        model_.controls.force_tol = to_double(f[0], line);
      if (f.size() > 2 && !f[2].empty())
        model_.controls.disp_tol = to_double(f[2], line);
    } else if (p.rfind("TIMEINCREMENTATION", 0) == 0 || p == "TIME INCREMENTATION") {
      if (!f.empty() && !f[0].empty())
        model_.controls.max_iterations = static_cast<int>(to_index(f[0], line));
    }
  }

  // *TIME POINTS: a list of increasing times within the step, used to force
  // increments to land on the given times.
  void time_points_data(const std::vector<std::string>& f, int line) {
    for (const std::string& tok : f)
      if (!tok.empty()) model_.time_points.times.push_back(to_double(tok, line));
  }

  // *AMPLITUDE, NAME=..., DEFINITION=TABULAR|STEP, TIME=..., VALUE=... — the data
  // lines carry (time, value) pairs. Periodic amplitudes are not a NAME-scoped card
  // here; PERIOD (if given as a parameter) sets the wrap length.
  void begin_amplitude(int line) {
    cur_amplitude_ = param("NAME");
    if (cur_amplitude_.empty()) throw ParseError(line, "*AMPLITUDE without NAME");
    Amplitude a;
    a.name = cur_amplitude_;
    const std::string def = param("DEFINITION");
    a.definition = (def == "STEP") ? Amplitude::Definition::Step
                                   : Amplitude::Definition::Tabular;
    if (!param("PERIOD").empty()) a.period = to_double(param("PERIOD"), line);
    model_.amplitudes[cur_amplitude_] = std::move(a);
  }

  // *AMPLITUDE data: a flat list of alternating time, value entries (possibly many
  // pairs per line), appended to the current amplitude curve.
  void amplitude_data(const std::vector<std::string>& f, int line) {
    if (cur_amplitude_.empty()) throw ParseError(line, "*AMPLITUDE data without NAME");
    Amplitude& a = model_.amplitudes[cur_amplitude_];
    for (std::size_t i = 0; i + 1 < f.size(); i += 2) {
      if (f[i].empty() || f[i + 1].empty()) continue;
      a.points.emplace_back(to_double(f[i], line), to_double(f[i + 1], line));
    }
  }

  // *CHANGE MATERIAL, NAME=<mat>: re-open an existing material to redefine its
  // properties within the step (followed by *CHANGE PLASTIC in CalculiX). The
  // material must already exist. Only re-opens the context here; the actual
  // property change is done by the following *CHANGE PLASTIC (deferred: no
  // plasticity yet). (spec: loads-and-boundary-conditions — step property change.)
  void begin_change_material(int line) {
    const std::string name = param("NAME");
    if (name.empty()) throw ParseError(line, "*CHANGE MATERIAL needs NAME=");
    if (model_.materials.find(name) == model_.materials.end())
      throw ParseError(line, "*CHANGE MATERIAL: unknown material '" + name + "'");
    cur_material_ = name;
  }

  // *PLASTIC[, HARDENING=ISOTROPIC|KINEMATIC|COMBINED]: begin the J2 plasticity
  // hardening table on the current material. Data lines are (yield, plastic_strain).
  void begin_plastic(int line) {
    if (cur_material_.empty()) throw ParseError(line, "*PLASTIC without *MATERIAL");
    const std::string h = param("HARDENING");
    cur_hardening_ = Plastic::Hardening::Isotropic;
    if (h == "KINEMATIC") cur_hardening_ = Plastic::Hardening::Kinematic;
    else if (h == "COMBINED") cur_hardening_ = Plastic::Hardening::Combined;
    else if (!h.empty() && h != "ISOTROPIC")
      throw ParseError(line, "*PLASTIC unsupported HARDENING=" + h);
    Plastic& p = ensure_plastic();
    p.hardening = cur_hardening_;
  }

  // *CYCLIC HARDENING: an alternative source of the isotropic hardening table (yield
  // vs. equivalent plastic strain). Reuses the same table on the current material.
  void begin_cyclic_hardening(int line) {
    if (cur_material_.empty())
      throw ParseError(line, "*CYCLIC HARDENING without *MATERIAL");
    ensure_plastic();
  }

  // Get (creating if needed) the Plastic block on the current material.
  Plastic& ensure_plastic() {
    auto& mat = model_.materials[cur_material_];
    if (!mat.plastic) mat.plastic = Plastic{};
    return *mat.plastic;
  }

  // A *PLASTIC / *CYCLIC HARDENING data line: yield_stress, equivalent_plastic_strain
  // [, temperature]. Plastic strain must be non-decreasing (a monotone table); the
  // first row's plastic strain is the yield at zero plastic strain.
  void plastic_data(const std::vector<std::string>& f, int line) {
    if (cur_material_.empty())
      throw ParseError(line, "*PLASTIC/*CYCLIC HARDENING without *MATERIAL");
    if (f.size() < 2)
      throw ParseError(line, "*PLASTIC needs yield, plastic_strain");
    const Real sy = to_double(f[0], line);
    const Real ep = to_double(f[1], line);
    Plastic& p = ensure_plastic();
    if (!p.eqplastic.empty() && ep < p.eqplastic.back())
      throw ParseError(line, "*PLASTIC plastic strain must be non-decreasing");
    p.eqplastic.push_back(ep);
    p.yield.push_back(sy);
    // Kinematic/combined hardening modulus from the table slope (linear kinematic).
    // For a two-point table this is the constant back-stress modulus H_kin; with a
    // single point (perfect plasticity) it stays 0.
    if ((p.hardening == Plastic::Hardening::Kinematic ||
         p.hardening == Plastic::Hardening::Combined) &&
        p.eqplastic.size() >= 2) {
      const std::size_t k = p.eqplastic.size();
      const Real dep = p.eqplastic[k - 1] - p.eqplastic[k - 2];
      if (dep > 0.0) p.kinematic_modulus = (p.yield[k - 1] - p.yield[k - 2]) / dep;
    }
  }

  // The connector kind of the elements in an elset (they must be homogeneous).
  // Throws if the elset is unknown or empty.
  ConnKind conn_elset_kind(const std::string& elset, int line) const {
    const auto it = conn_elsets_.find(elset);
    if (it == conn_elsets_.end() || it->second.empty())
      throw ParseError(line, "connector ELSET '" + elset + "' has no elements");
    return conn_elems_.at(it->second.front()).kind;
  }

  // *SPRING data. Depending on the referenced element type:
  //   SPRINGA: one line, the (linear) stiffness.
  //   SPRING1: first line the acting DOF, second line the stiffness.
  //   SPRING2: first line the two DOFs, second line the stiffness.
  // (NONLINEAR force-displacement tables are not supported yet; linear only.)
  void spring_data(const std::vector<std::string>& f, int line) {
    const std::string elset = param("ELSET");
    const ConnKind kind = conn_elset_kind(elset, line);
    if (kind == ConnKind::SpringA) {
      add_springs(elset, to_double(f.at(0), line), 0, 0, line);
    } else if (kind == ConnKind::Spring1) {
      if (!spring_dofs_read_) {
        spring_dof1_ = static_cast<int>(to_index(f.at(0), line));
        spring_dofs_read_ = true;
      } else {
        add_springs(elset, to_double(f.at(0), line), spring_dof1_, 0, line);
      }
    } else if (kind == ConnKind::Spring2) {
      if (!spring_dofs_read_) {
        spring_dof1_ = static_cast<int>(to_index(f.at(0), line));
        spring_dof2_ = static_cast<int>(to_index(f.at(1), line));
        spring_dofs_read_ = true;
      } else {
        add_springs(elset, to_double(f.at(0), line), spring_dof1_, spring_dof2_, line);
      }
    } else {
      throw ParseError(line, "*SPRING elset '" + elset + "' is not a spring element");
    }
  }

  // Emit a Spring for every connector element in the elset with the given stiffness
  // and (for SPRING1/2) DOFs.
  void add_springs(const std::string& elset, Real k, int dof1, int dof2, int line) {
    for (const Index eid : conn_elsets_.at(elset)) {
      const ConnElem& ce = conn_elems_.at(eid);
      Spring sp;
      sp.stiffness = k;
      if (ce.kind == ConnKind::SpringA) {
        sp.kind = Spring::Kind::Axial;
        sp.node1 = ce.nodes.at(0);
        sp.node2 = ce.nodes.at(1);
      } else if (ce.kind == ConnKind::Spring1) {
        sp.kind = Spring::Kind::Grounded;
        sp.node1 = ce.nodes.at(0);
        sp.dof1 = dof1;
      } else {  // Spring2
        sp.kind = Spring::Kind::Dof;
        sp.node1 = ce.nodes.at(0);
        sp.node2 = ce.nodes.at(1);
        sp.dof1 = dof1;
        sp.dof2 = dof2;
      }
      model_.springs.push_back(sp);
    }
    (void)line;
  }

  // *DASHPOT data (viscous damping coefficient), stored for dynamics. Same
  // DOF-line convention as *SPRING for DASHPOT1/DASHPOT2.
  void dashpot_data(const std::vector<std::string>& f, int line) {
    const std::string elset = param("ELSET");
    const ConnKind kind = conn_elset_kind(elset, line);
    Real c = 0.0;
    if (kind == ConnKind::DashpotA) {
      c = to_double(f.at(0), line);
    } else if (kind == ConnKind::Dashpot1 || kind == ConnKind::Dashpot2) {
      if (!spring_dofs_read_) {
        spring_dof1_ = static_cast<int>(to_index(f.at(0), line));
        if (kind == ConnKind::Dashpot2) spring_dof2_ = static_cast<int>(to_index(f.at(1), line));
        spring_dofs_read_ = true;
        return;
      }
      c = to_double(f.at(0), line);
    } else {
      throw ParseError(line, "*DASHPOT elset '" + elset + "' is not a dashpot element");
    }
    for (const Index eid : conn_elsets_.at(elset)) {
      const ConnElem& ce = conn_elems_.at(eid);
      Dashpot d;
      d.coefficient = c;
      d.node1 = ce.nodes.at(0);
      if (ce.kind == ConnKind::DashpotA) {
        d.kind = Dashpot::Kind::Axial;
        d.node2 = ce.nodes.at(1);
      } else if (ce.kind == ConnKind::Dashpot1) {
        d.kind = Dashpot::Kind::Grounded;
        d.dof1 = spring_dof1_;
      } else {
        d.kind = Dashpot::Kind::Dof;
        d.node2 = ce.nodes.at(1);
        d.dof1 = spring_dof1_;
        d.dof2 = spring_dof2_;
      }
      model_.dashpots.push_back(d);
    }
  }

  // *MASS data (point mass magnitude), stored for dynamics.
  void mass_data(const std::vector<std::string>& f, int line) {
    const std::string elset = param("ELSET");
    const Real mass = to_double(f.at(0), line);
    const auto it = conn_elsets_.find(elset);
    if (it == conn_elsets_.end())
      throw ParseError(line, "*MASS ELSET '" + elset + "' has no elements");
    for (const Index eid : it->second) {
      const ConnElem& ce = conn_elems_.at(eid);
      model_.point_masses.push_back(PointMass{ce.nodes.at(0), mass});
    }
  }

  // *EQUATION data. The first data line after the card is the term count N; the
  // subsequent line(s) carry N (node, dof, coeff) triplets (possibly continued). The
  // first triplet is the dependent DOF. (spec: constraints — Linear equations.)
  void equation_data(const std::vector<std::string>& f, int line) {
    if (eq_terms_remaining_ < 0) {
      if (f.empty()) throw ParseError(line, "*EQUATION needs a term count");
      eq_terms_remaining_ = static_cast<int>(to_index(f[0], line));
      cur_equation_ = Equation{};
      cur_equation_.origin = "*EQUATION";
      return;
    }
    for (std::size_t i = 0; i + 3 <= f.size() && eq_terms_remaining_ > 0; i += 3) {
      const Index nd = to_index(f[i], line);
      const int dof = static_cast<int>(to_index(f[i + 1], line));
      const Real c = to_double(f[i + 2], line);
      cur_equation_.terms.push_back(EquationTerm{nd, dof, c});
      --eq_terms_remaining_;
    }
    if (eq_terms_remaining_ == 0) {
      model_.equations.push_back(std::move(cur_equation_));
      eq_terms_remaining_ = -1;  // next data line starts a new equation
    }
  }

  // *MPC data: "MPCTYPE, node, node, ...". First token is the type (BEAM/PLANE/
  // STRAIGHT or a user name); the rest are node ids / node sets in card order (first
  // is dependent). (spec: constraints — Multi-point constraints.)
  void mpc_data(const std::vector<std::string>& f, int line) {
    if (f.empty()) throw ParseError(line, "*MPC needs a type and nodes");
    Mpc m;
    const std::string type = upper(f[0]);
    if (type == "BEAM") m.kind = Mpc::Kind::Beam;
    else if (type == "PLANE") m.kind = Mpc::Kind::Plane;
    else if (type == "STRAIGHT") m.kind = Mpc::Kind::Straight;
    else { m.kind = Mpc::Kind::User; m.user_name = type; }
    for (std::size_t i = 1; i < f.size(); ++i)
      if (!f[i].empty())
        for (const Index nd : node_refs(f[i], line)) m.nodes.push_back(nd);
    model_.mpcs.push_back(std::move(m));
  }

  // *RIGID BODY, NSET=|ELSET=, REF NODE=, ROT NODE=. Ties the set to the reference
  // node (translation) and optional rotation node. (spec: constraints — Rigid bodies.)
  void begin_rigid_body(int line) {
    RigidBody rb;
    const std::string nset = param("NSET");
    const std::string elset = param("ELSET");
    if (!nset.empty()) rb.nset = nset;
    else if (!elset.empty()) rb.nodes = elset_nodes(elset, line);
    else throw ParseError(line, "*RIGID BODY needs NSET= or ELSET=");
    const std::string ref = param("REFNODE");
    if (ref.empty()) throw ParseError(line, "*RIGID BODY needs REF NODE=");
    rb.ref_node = to_index(ref, line);
    const std::string rot = param("ROTNODE");
    if (!rot.empty()) rb.rot_node = to_index(rot, line);
    model_.rigid_bodies.push_back(std::move(rb));
  }

  // Collect the node ids of an element set (union of the elements' nodes), for a
  // *RIGID BODY on an ELSET.
  std::vector<Index> elset_nodes(const std::string& elset, int line) {
    const auto it = elsets_.find(upper(elset));
    if (it == elsets_.end()) throw ParseError(line, "unknown element set '" + elset + "'");
    std::vector<Index> out;
    std::unordered_map<Index, bool> seen;
    for (const Index eid : it->second) {
      const Index ei = model_.mesh.element_index(eid);
      if (ei < 0) continue;
      for (const Index nd : model_.mesh.elements()[static_cast<std::size_t>(ei)].nodes)
        if (!seen[nd]) { seen[nd] = true; out.push_back(nd); }
    }
    return out;
  }

  // *COUPLING, REF NODE=, SURFACE=. The kind (KINEMATIC / DISTRIBUTING) comes from a
  // following modifier card. (spec: constraints — couplings.)
  void begin_coupling(int line) {
    Coupling cp;
    const std::string ref = param("REFNODE");
    if (ref.empty()) throw ParseError(line, "*COUPLING needs REF NODE=");
    cp.ref_node = to_index(ref, line);
    cp.surface = param("SURFACE");
    if (cp.surface.empty()) throw ParseError(line, "*COUPLING needs SURFACE=");
    cp.nodes = surface_nodes(cp.surface, line);
    cp.dofs.clear();  // filled by the modifier's data lines (default all 3 if none)
    cur_coupling_ = static_cast<int>(model_.couplings.size());
    model_.couplings.push_back(std::move(cp));
  }

  // Set the kind of the coupling currently being built from a *KINEMATIC /
  // *DISTRIBUTING modifier card.
  void set_coupling_kind(int line) {
    if (cur_coupling_ < 0)
      throw ParseError(line, std::string(card_ == "*KINEMATIC" ? "*KINEMATIC" : "*DISTRIBUTING") +
                                 " outside a *COUPLING");
    model_.couplings[static_cast<std::size_t>(cur_coupling_)].kind =
        card_ == "*KINEMATIC" ? Coupling::Kind::Kinematic : Coupling::Kind::Distributing;
  }

  // Coupling DOF-range data line: "first_dof, last_dof" (or a single dof). Appends the
  // constrained DOFs to the current coupling. Applies to *COUPLING / *KINEMATIC /
  // *DISTRIBUTING data lines.
  void coupling_data(const std::vector<std::string>& f, int line) {
    if (cur_coupling_ < 0 || f.empty()) return;
    Coupling& cp = model_.couplings[static_cast<std::size_t>(cur_coupling_)];
    const int d1 = static_cast<int>(to_index(f[0], line));
    const int d2 = f.size() > 1 && !f[1].empty() ? static_cast<int>(to_index(f[1], line)) : d1;
    for (int c = d1; c <= d2 && c <= kDofsPerNode; ++c) cp.dofs.push_back(c);
  }

  // *DISTRIBUTING COUPLING, ELSET=|REF NODE=... The data lines carry the coupled node
  // ids and weights: "node, weight". Simplified to equal-weight averaging here.
  void begin_distributing_coupling(int line) {
    Coupling cp;
    cp.kind = Coupling::Kind::Distributing;
    const std::string ref = param("REFNODE");
    if (!ref.empty()) cp.ref_node = to_index(ref, line);
    cp.dofs = {1, 2, 3};
    cur_coupling_ = static_cast<int>(model_.couplings.size());
    model_.couplings.push_back(std::move(cp));
  }

  void distributing_coupling_data(const std::vector<std::string>& f, int line) {
    if (cur_coupling_ < 0 || f.empty()) return;
    Coupling& cp = model_.couplings[static_cast<std::size_t>(cur_coupling_)];
    // "node_or_nset[, weight]" — the reference node is the first if REF NODE= absent.
    for (const Index nd : node_refs(f[0], line)) {
      if (cp.ref_node == 0) { cp.ref_node = nd; continue; }
      cp.nodes.push_back(nd);
    }
  }

  // *TIE, NAME=, POSITION TOLERANCE=. The single data line carries "slave, master"
  // surface names. (spec: constraints — Surface ties.)
  void begin_tie(int line) {
    tie_pending_ = Tie{};
    tie_pending_.name = param("NAME");
    const std::string tol = param("POSITIONTOLERANCE");
    if (!tol.empty()) tie_pending_.tolerance = to_double(tol, line);
    tie_has_pending_ = true;
  }

  // *TIE data line: "slave_surface, master_surface". Surfaces are resolved to node
  // lists; matching-mesh ties pair coincident nodes at expansion time.
  void tie_data(const std::vector<std::string>& f, int line) {
    if (!tie_has_pending_) throw ParseError(line, "*TIE data without a *TIE card");
    if (f.size() < 2) throw ParseError(line, "*TIE needs slave, master surfaces");
    tie_pending_.slave_surface = f[0];
    tie_pending_.master_surface = f[1];
    tie_pending_.slave_nodes = surface_nodes(f[0], line);
    tie_pending_.master_nodes = surface_nodes(f[1], line);
    model_.ties.push_back(std::move(tie_pending_));
    tie_has_pending_ = false;
  }

  // *MODEL CHANGE card: TYPE=ELEMENT (default) or TYPE=CONTACT PAIR, with ADD or
  // REMOVE selecting activation/deactivation. The affected element ids / contact-pair
  // surfaces come on the data lines. (spec: model-change — element/contact-pair
  // birth-death.) Within the current single-step model, a REMOVE deactivates the
  // element for the whole step; an ADD (re)activates it strain-free (the solve starts
  // undeformed). True cross-step birth-death needs the multi-step engine (tasks 5.1).
  void begin_model_change(int line) {
    const std::string type = param("TYPE");
    if (type.empty() || type == "ELEMENT") {
      mc_type_ = ModelChangeType::Element;
    } else if (type == "CONTACTPAIR" || type == "CONTACT PAIR") {
      mc_type_ = ModelChangeType::ContactPair;
    } else {
      throw ParseError(line, "*MODEL CHANGE TYPE must be ELEMENT or CONTACT PAIR, got '" +
                                 type + "'");
    }
    // ADD / REMOVE given as bare parameters (CalculiX). Exactly one is required.
    const bool has_add = params_.count("ADD") != 0;
    const bool has_remove = params_.count("REMOVE") != 0;
    if (has_add == has_remove)
      throw ParseError(line, "*MODEL CHANGE needs exactly one of ADD or REMOVE");
    mc_add_ = has_add;
  }

  // *MODEL CHANGE data line. For TYPE=ELEMENT: an element id or elset name per token;
  // REMOVE deactivates them, ADD reactivates (drops them from the deactivated set).
  // For TYPE=CONTACT PAIR: "surface_a, surface_b" recorded for the contact workflow.
  void model_change_data(const std::vector<std::string>& f, int line) {
    if (f.empty() || f[0].empty()) return;
    if (mc_type_ == ModelChangeType::ContactPair) {
      if (f.size() < 2)
        throw ParseError(line, "*MODEL CHANGE, TYPE=CONTACT PAIR needs two surfaces");
      model_.contact_pair_changes.push_back(
          ContactPairChange{f[0], f[1], mc_add_});
      return;
    }
    for (const std::string& tok : f) {
      if (tok.empty()) continue;
      for (const Index eid : elem_refs(tok, line)) {
        if (mc_add_)
          std::erase(model_.deactivated_elements, eid);
        else if (std::find(model_.deactivated_elements.begin(),
                           model_.deactivated_elements.end(),
                           eid) == model_.deactivated_elements.end())
          model_.deactivated_elements.push_back(eid);
      }
    }
  }

  // Resolve a *SURFACE (or nset) name to its node ids. Element surfaces contribute the
  // nodes of their listed faces; node surfaces contribute their node list; an nset
  // name contributes its nodes.
  std::vector<Index> surface_nodes(const std::string& name, int line) {
    const std::string key = upper(name);
    // A node set with this name.
    const auto ns = nsets_.find(key);
    if (ns != nsets_.end()) return ns->second;
    // A *SURFACE built so far.
    for (const Surface& s : surfaces_)
      if (upper(s.name) == key) {
        if (s.type == Surface::Type::Node) return s.nodes;
        std::vector<Index> out;
        std::unordered_map<Index, bool> seen;
        for (const auto& [eid, face] : s.faces) {
          const Index ei = model_.mesh.element_index(eid);
          if (ei < 0) continue;
          for (const Index nd :
               model_.mesh.elements()[static_cast<std::size_t>(ei)].nodes)
            if (!seen[nd]) { seen[nd] = true; out.push_back(nd); }
        }
        return out;
      }
    throw ParseError(line, "constraint references unknown surface/nset '" + name + "'");
  }

  bool is_data_card() const {
    static const std::vector<std::string> data_cards = {
        "*NODE",    "*ELEMENT",  "*NSET",         "*ELSET",     "*MATERIAL",
        "*ELASTIC", "*DENSITY",  "*SOLIDSECTION", "*BOUNDARY",  "*CLOAD",
        "*DLOAD",   "*DSLOAD",   "*SURFACE",      "*STATIC",    "*CONTROLS",
        "*TIMEPOINTS", "*AMPLITUDE", "*CHANGEPLASTIC", "*SPRING", "*MASS",
        "*DASHPOT", "*PLASTIC", "*CYCLICHARDENING", "*HYPERELASTIC",
        "*USERMATERIAL", "*DEPVAR", "*RATEDEPENDENT", "*EQUATION", "*MPC",
        "*RIGIDBODY", "*COUPLING", "*DISTRIBUTINGCOUPLING", "*KINEMATIC",
        "*DISTRIBUTING", "*TIE", "*MODELCHANGE",
        // Contact (Phase 3) data cards with a data line.
        "*GAPCONDUCTANCE", "*GAPHEATGENERATION",
        // Thermal (Phase 3) data cards.
        "*HEATTRANSFER", "*CONDUCTIVITY", "*SPECIFICHEAT", "*EXPANSION", "*CFLUX",
        "*DFLUX", "*TEMPERATURE", "*FILM", "*RADIATE", "*INITIALCONDITIONS",
        "*FREQUENCY",
        // Dynamics (Phase 4) data cards with a data line.
        "*DYNAMIC", "*MODALDYNAMIC", "*STEADYSTATEDYNAMICS", "*COMPLEXFREQUENCY",
        "*MODALDAMPING", "*BASEMOTION", "*RETAINEDNODALDOFS",
        "*COUPLEDTEMPERATURE-DISPLACEMENT", "*COUPLEDTEMPERATUREDISPLACEMENT"};
    return std::find(data_cards.begin(), data_cards.end(), card_) != data_cards.end();
  }

  void data_line(const std::vector<std::string>& f, int line) {
    if (card_ == "*NODE") return node(f, line);
    if (card_ == "*ELEMENT") return element(f, line);
    if (card_ == "*NSET") return set_data(f, line, nsets_, param("NSET"));
    if (card_ == "*ELSET") return set_data(f, line, elsets_, param("ELSET"));
    if (card_ == "*ELASTIC") return elastic(f, line);
    if (card_ == "*PLASTIC" || card_ == "*CYCLICHARDENING")
      return plastic_data(f, line);
    if (card_ == "*HYPERELASTIC") return hyperelastic_data(f, line);
    if (card_ == "*USERMATERIAL") return user_material_data(f, line);
    if (card_ == "*DEPVAR") return depvar_data(f, line);
    if (card_ == "*RATEDEPENDENT") return ratedependent_data(f, line);
    if (card_ == "*DENSITY") return density(f, line);
    if (card_ == "*BOUNDARY") return boundary(f, line);
    if (card_ == "*CLOAD") return cload(f, line);
    if (card_ == "*DLOAD") return dload(f, line);
    if (card_ == "*DSLOAD") return dsload(f, line);
    if (card_ == "*SURFACE") return surface_data(f, line);
    if (card_ == "*SURFACEBEHAVIOR") return surface_behavior_data(f, line);
    if (card_ == "*FRICTION") return friction_data(f, line);
    if (card_ == "*GAPCONDUCTANCE") return gap_conductance_data(f, line);
    if (card_ == "*GAPHEATGENERATION") return gap_heat_generation_data(f, line);
    if (card_ == "*CONTACTPAIR") return contact_pair_data(f, line);
    if (card_ == "*STATIC") return static_data(f, line);
    if (card_ == "*HEATTRANSFER") return static_data(f, line);  // same inc data line
    if (card_ == "*FREQUENCY") return frequency_data(f, line);
    if (card_ == "*BUCKLE") return buckle_data(f, line);
    if (card_ == "*DYNAMIC") return modal_dynamic_data(f, line);  // "dt, t_end"
    if (card_ == "*MODALDYNAMIC") return modal_dynamic_data(f, line);
    if (card_ == "*STEADYSTATEDYNAMICS") return steady_state_data(f, line);
    if (card_ == "*COMPLEXFREQUENCY") return complex_frequency_data(f, line);
    if (card_ == "*MODALDAMPING") return modal_damping_data(f, line);
    if (card_ == "*BASEMOTION") return base_motion_data(f, line);
    if (card_ == "*RETAINEDNODALDOFS") return retained_nodal_dofs_data(f, line);
    if (card_ == "*DESIGNVARIABLES") return design_variables_data(f, line);
    if (card_ == "*DESIGNRESPONSE" || card_ == "*OBJECTIVE")
      return design_response_data(f, line);
    if (card_ == "*COUPLEDTEMPERATURE-DISPLACEMENT" ||
        card_ == "*COUPLEDTEMPERATUREDISPLACEMENT")
      return static_data(f, line);  // same increment data line as *STATIC
    if (card_ == "*CONDUCTIVITY") return conductivity_data(f, line);
    if (card_ == "*SPECIFICHEAT") return specific_heat_data(f, line);
    if (card_ == "*EXPANSION") return expansion_data(f, line);
    if (card_ == "*CFLUX") return cflux_data(f, line);
    if (card_ == "*DFLUX") return dflux_data(f, line);
    if (card_ == "*FILM") return film_data(f, line);
    if (card_ == "*RADIATE") return radiate_data(f, line);
    if (card_ == "*INITIALCONDITIONS") return initial_conditions_data(f, line);
    if (card_ == "*TEMPERATURE") return temperature_data(f, line);
    if (card_ == "*CONTROLS") return controls_data(f, line);
    if (card_ == "*TIMEPOINTS") return time_points_data(f, line);
    if (card_ == "*AMPLITUDE") return amplitude_data(f, line);
    if (card_ == "*SPRING") return spring_data(f, line);
    if (card_ == "*DASHPOT") return dashpot_data(f, line);
    if (card_ == "*MASS") return mass_data(f, line);
    if (card_ == "*EQUATION") return equation_data(f, line);
    if (card_ == "*MPC") return mpc_data(f, line);
    if (card_ == "*RIGIDBODY") return;  // all info on the card
    if (card_ == "*COUPLING" || card_ == "*KINEMATIC" || card_ == "*DISTRIBUTING")
      return coupling_data(f, line);
    if (card_ == "*DISTRIBUTINGCOUPLING") return distributing_coupling_data(f, line);
    if (card_ == "*TIE") return tie_data(f, line);
    if (card_ == "*MODELCHANGE") return model_change_data(f, line);
    // *SOLID SECTION optional thickness line, output requests, amplitudes:
    // accepted and ignored in Phase 1.
  }

  void node(const std::vector<std::string>& f, int line) {
    const Index id = to_index(f.at(0), line);
    Vec3 x{0, 0, 0};
    for (std::size_t k = 1; k < f.size() && k <= 3; ++k)
      x[k - 1] = to_double(f[k], line);
    model_.mesh.add_node(id, x);
    if (!param("NSET").empty()) nsets_[param("NSET")].push_back(id);
  }

  // Map a connector element TYPE= to its kind and expected node count, or return
  // false if TYPE is not a connector element.
  static bool connector_kind(const std::string& type, ConnKind& kind, int& nnodes) {
    if (type == "SPRINGA") { kind = ConnKind::SpringA; nnodes = 2; return true; }
    if (type == "SPRING1") { kind = ConnKind::Spring1; nnodes = 1; return true; }
    if (type == "SPRING2") { kind = ConnKind::Spring2; nnodes = 2; return true; }
    if (type == "MASS") { kind = ConnKind::Mass; nnodes = 1; return true; }
    if (type == "DASHPOTA") { kind = ConnKind::DashpotA; nnodes = 2; return true; }
    if (type == "DASHPOT1") { kind = ConnKind::Dashpot1; nnodes = 1; return true; }
    if (type == "DASHPOT2") { kind = ConnKind::Dashpot2; nnodes = 2; return true; }
    return false;
  }

  void element(const std::vector<std::string>& f, int line) {
    const std::string type = param("TYPE");
    ConnKind ck{};
    int cn = 0;
    if (connector_kind(type, ck, cn)) return connector_element(f, line, ck, cn);

    ElementType et{};
    if (!parse_element_type(type, et))
      throw ParseError(line, "unsupported element TYPE=" + type);
    const Index id = to_index(f.at(0), line);
    std::vector<Index> conn;
    for (std::size_t k = 1; k < f.size(); ++k) conn.push_back(to_index(f[k], line));
    if (static_cast<int>(conn.size()) != nodes_per_element(et))
      throw ParseError(line, "element " + std::to_string(id) + ": expected " +
                                 std::to_string(nodes_per_element(et)) +
                                 " nodes, got " + std::to_string(conn.size()));
    model_.mesh.add_element(id, et, conn);
    if (!param("ELSET").empty()) elsets_[param("ELSET")].push_back(id);
  }

  // A connector *ELEMENT line: "elem_id, node[, node2]". Stored aside (not in the
  // Mesh) until the matching *SPRING/*MASS/*DASHPOT card attaches its property.
  void connector_element(const std::vector<std::string>& f, int line, ConnKind kind,
                         int nnodes) {
    const Index id = to_index(f.at(0), line);
    std::vector<Index> conn;
    for (std::size_t k = 1; k < f.size(); ++k) conn.push_back(to_index(f[k], line));
    if (static_cast<int>(conn.size()) != nnodes)
      throw ParseError(line, "connector element " + std::to_string(id) +
                                 ": expected " + std::to_string(nnodes) +
                                 " node(s), got " + std::to_string(conn.size()));
    conn_elems_[id] = ConnElem{kind, conn};
    if (!param("ELSET").empty()) conn_elsets_[param("ELSET")].push_back(id);
  }

  void set_data(const std::vector<std::string>& f, int line,
                std::unordered_map<std::string, std::vector<Index>>& sets,
                const std::string& name) {
    if (name.empty()) throw ParseError(line, "set card without NSET/ELSET name");
    auto& set = sets[name];
    if (params_.count("GENERATE") != 0) {
      const Index a = to_index(f.at(0), line), b = to_index(f.at(1), line);
      const Index inc = f.size() > 2 ? to_index(f[2], line) : 1;
      for (Index v = a; v <= b; v += inc) set.push_back(v);
      return;
    }
    for (const std::string& tok : f) {
      // Integer id, or the name of a previously-defined set to expand.
      if (!tok.empty() && (std::isdigit(static_cast<unsigned char>(tok.front())) ||
                           tok.front() == '-' || tok.front() == '+')) {
        set.push_back(to_index(tok, line));
      } else {
        const std::string ref = upper(tok);
        const auto it = sets.find(ref);
        if (it == sets.end()) throw ParseError(line, "unknown set '" + tok + "'");
        set.insert(set.end(), it->second.begin(), it->second.end());
      }
    }
  }

  // *HYPERELASTIC[, N=1]: compressible neo-Hookean (spec: material-models 4.3). The
  // data line carries the coefficients; for N=1 (the only order implemented) that is
  // C10, D1 (CalculiX order). mu=2*C10, kappa=2/D1 are derived at model build time.
  void begin_hyperelastic(int line) {
    if (cur_material_.empty())
      throw ParseError(line, "*HYPERELASTIC without *MATERIAL");
    const std::string n = param("N");
    if (!n.empty() && n != "1")
      throw ParseError(line, "*HYPERELASTIC only N=1 (neo-Hookean) is implemented");
    auto& mat = model_.materials[cur_material_];
    if (!mat.hyperelastic) mat.hyperelastic = Hyperelastic{};
  }

  // *HYPERELASTIC data line: C10, D1.
  void hyperelastic_data(const std::vector<std::string>& f, int line) {
    if (cur_material_.empty())
      throw ParseError(line, "*HYPERELASTIC without *MATERIAL");
    if (f.empty()) throw ParseError(line, "*HYPERELASTIC needs C10[, D1]");
    auto& mat = model_.materials[cur_material_];
    if (!mat.hyperelastic) mat.hyperelastic = Hyperelastic{};
    Hyperelastic& h = *mat.hyperelastic;
    h.model = Hyperelastic::Model::NeoHookean;
    h.c10 = to_double(f[0], line);
    h.d1 = f.size() > 1 && !f[1].empty() ? to_double(f[1], line) : 0.0;
    h.mu = 2.0 * h.c10;
    h.kappa = h.d1 != 0.0 ? 2.0 / h.d1 : 0.0;
    if (h.d1 == 0.0)
      throw ParseError(line, "*HYPERELASTIC needs D1>0 (compressible); "
                             "near-incompressible u/p is deferred");
  }

  // *USER MATERIAL[, CONSTANTS=n]: name a registered C++ user-material factory
  // (spec: material-models 4.6). The registry key is NAME= if given, else the
  // material name. Data lines carry the `n` constants.
  void begin_user_material(int line) {
    if (cur_material_.empty())
      throw ParseError(line, "*USER MATERIAL without *MATERIAL");
    auto& mat = model_.materials[cur_material_];
    if (!mat.user) mat.user = UserMaterial{};
    // Registry key: explicit NAME= wins, else the material name (upper-cased).
    const std::string name = param("NAME");
    mat.user->name = name.empty() ? upper(cur_material_) : name;
  }

  // *USER MATERIAL data line: the material constants, appended in order.
  void user_material_data(const std::vector<std::string>& f, int line) {
    if (cur_material_.empty())
      throw ParseError(line, "*USER MATERIAL data without *MATERIAL");
    auto& mat = model_.materials[cur_material_];
    if (!mat.user) mat.user = UserMaterial{};
    for (const std::string& tok : f)
      if (!tok.empty()) mat.user->constants.push_back(to_double(tok, line));
  }

  // *DEPVAR: the number of solution-dependent state variables per integration point,
  // read from the single data line (spec: material-models 4.6).
  void begin_depvar(int line) {
    if (cur_material_.empty()) throw ParseError(line, "*DEPVAR without *MATERIAL");
    auto& mat = model_.materials[cur_material_];
    if (!mat.user) mat.user = UserMaterial{};
  }
  void depvar_data(const std::vector<std::string>& f, int line) {
    if (cur_material_.empty()) throw ParseError(line, "*DEPVAR without *MATERIAL");
    auto& mat = model_.materials[cur_material_];
    if (!mat.user) mat.user = UserMaterial{};
    if (!f.empty() && !f[0].empty())
      mat.user->ndepvar = static_cast<int>(to_index(f[0], line));
  }

  // *RATEDEPENDENT: a rate-scaling factor passed to the user material
  // (spec: material-models 4.6). The first data value is stored as the scale.
  void begin_ratedependent(int line) {
    if (cur_material_.empty())
      throw ParseError(line, "*RATEDEPENDENT without *MATERIAL");
    auto& mat = model_.materials[cur_material_];
    if (!mat.user) mat.user = UserMaterial{};
  }
  void ratedependent_data(const std::vector<std::string>& f, int line) {
    if (cur_material_.empty())
      throw ParseError(line, "*RATEDEPENDENT without *MATERIAL");
    auto& mat = model_.materials[cur_material_];
    if (!mat.user) mat.user = UserMaterial{};
    if (!f.empty() && !f[0].empty())
      mat.user->rate_scale = to_double(f[0], line);
  }

  void elastic(const std::vector<std::string>& f, int line) {
    if (cur_material_.empty()) throw ParseError(line, "*ELASTIC without *MATERIAL");
    const Real E = to_double(f.at(0), line);
    const Real nu = f.size() > 1 ? to_double(f[1], line) : 0.0;
    model_.materials[cur_material_].elastic = ElasticIso{E, nu};
  }

  void density(const std::vector<std::string>& f, int line) {
    if (cur_material_.empty()) throw ParseError(line, "*DENSITY without *MATERIAL");
    model_.materials[cur_material_].density = to_double(f.at(0), line);
  }

  // Expand a node reference (integer id or nset name) to node ids.
  std::vector<Index> node_refs(const std::string& tok, int line) {
    if (!tok.empty() && (std::isdigit(static_cast<unsigned char>(tok.front())) ||
                         tok.front() == '-' || tok.front() == '+'))
      return {to_index(tok, line)};
    const auto it = nsets_.find(upper(tok));
    if (it == nsets_.end()) throw ParseError(line, "unknown node set '" + tok + "'");
    return it->second;
  }

  void boundary(const std::vector<std::string>& f, int line) {
    if (f.size() < 2) throw ParseError(line, "*BOUNDARY needs node,dof[,dof2[,value]]");
    const int d1 = static_cast<int>(to_index(f[1], line));
    const int d2 = f.size() > 2 && !f[2].empty() ? static_cast<int>(to_index(f[2], line)) : d1;
    const Real val = f.size() > 3 ? to_double(f[3], line) : 0.0;
    const std::string amp = param("AMPLITUDE");
    // *BOUNDARY, FIXED holds each listed DOF at its currently-attained (deformed)
    // value; the multi-step driver resolves that at solve time (see Spc::fixed).
    const bool fixed = params_.count("FIXED") != 0;
    // Temperature DOF: CalculiX numbers the nodal temperature as DOF 11 (and accepts
    // 0 as the temperature DOF in *CFLUX). In a heat step, *BOUNDARY on DOF 11 (or 0)
    // prescribes a temperature; mechanical DOFs 1..3 still prescribe displacements.
    const bool temp_dof = (d1 == 11) || (heat_step_ && d1 == 0);
    for (const Index nd : node_refs(f[0], line)) {
      if (temp_dof) {
        model_.temp_bcs.push_back(TempBc{nd, val, amp});
        continue;
      }
      for (int c = d1; c <= d2 && c <= kDofsPerNode; ++c)
        model_.spcs.push_back(Spc{nd, c, val, amp, fixed});
    }
  }

  void cload(const std::vector<std::string>& f, int line) {
    if (f.size() < 3) throw ParseError(line, "*CLOAD needs node,dof,value");
    const int dof = static_cast<int>(to_index(f[1], line));
    const Real val = to_double(f[2], line);
    const std::string amp = param("AMPLITUDE");
    for (const Index nd : node_refs(f[0], line))
      model_.cloads.push_back(Cload{nd, dof, val, amp});
  }

  // Expand an element reference (integer id or elset name) to element ids.
  std::vector<Index> elem_refs(const std::string& tok, int line) {
    if (!tok.empty() && (std::isdigit(static_cast<unsigned char>(tok.front())) ||
                         tok.front() == '-' || tok.front() == '+'))
      return {to_index(tok, line)};
    const auto it = elsets_.find(upper(tok));
    if (it == elsets_.end()) throw ParseError(line, "unknown element set '" + tok + "'");
    return it->second;
  }

  // *DLOAD: pressure "elem, P<face>, magnitude"; body loads "GRAV, g, nx,ny,nz"
  // and "elset, CENTRIF, omega2, px,py,pz, ax,ay,az".
  void dload(const std::vector<std::string>& f, int line) {
    if (f.size() < 2) throw ParseError(line, "*DLOAD needs at least a label");
    const std::string label = upper(f[1]);
    const std::string amp = param("AMPLITUDE");
    if (label == "GRAV") return gravity_load(f, line, amp);
    if (label == "CENTRIF") return centrif_load(f, line, amp);
    if (label.empty() || label[0] != 'P')
      throw ParseError(line, "unsupported *DLOAD type '" + f[1] +
                                 "' (supports P<face>, GRAV, CENTRIF)");
    if (f.size() < 3) throw ParseError(line, "*DLOAD needs element,P<face>,magnitude");
    const int face = static_cast<int>(to_index(label.substr(1), line));
    const Real p = to_double(f[2], line);
    for (const Index eid : elem_refs(f[0], line))
      model_.dloads.push_back(Dload{eid, face, p, amp});
  }

  // *DLOAD GRAV: "elset, GRAV, g, nx, ny, nz". A blank/numeric elset -> all elements.
  void gravity_load(const std::vector<std::string>& f, int line,
                    const std::string& amp) {
    if (f.size() < 6) throw ParseError(line, "*DLOAD GRAV needs elset, g, nx, ny, nz");
    BodyLoad bl;
    bl.kind = BodyLoad::Kind::Gravity;
    bl.elset = elset_name(f[0]);
    bl.magnitude = to_double(f[2], line);
    bl.dir = {to_double(f[3], line), to_double(f[4], line), to_double(f[5], line)};
    bl.amplitude = amp;
    model_.body_loads.push_back(std::move(bl));
  }

  // *DLOAD CENTRIF: "elset, CENTRIF, omega2, px, py, pz, ax, ay, az".
  void centrif_load(const std::vector<std::string>& f, int line,
                    const std::string& amp) {
    if (f.size() < 9)
      throw ParseError(line, "*DLOAD CENTRIF needs omega2, px,py,pz, ax,ay,az");
    BodyLoad bl;
    bl.kind = BodyLoad::Kind::Centrifugal;
    bl.elset = elset_name(f[0]);
    bl.magnitude = to_double(f[2], line);  // omega^2
    bl.point = {to_double(f[3], line), to_double(f[4], line), to_double(f[5], line)};
    bl.dir = {to_double(f[6], line), to_double(f[7], line), to_double(f[8], line)};
    bl.amplitude = amp;
    model_.body_loads.push_back(std::move(bl));
  }

  // Interpret a *DLOAD/CENTRIF first field as an elset name (empty -> all elements).
  // A bare integer is treated as a singleton — but body loads are set-scoped, so an
  // integer element id is uncommon; keep the name form (empty means all).
  std::string elset_name(const std::string& tok) const {
    const std::string t = upper(trim(tok));
    if (t.empty() || std::isdigit(static_cast<unsigned char>(t.front())) != 0)
      return {};  // blank or numeric -> apply to all elements
    return t;
  }

  // *DSLOAD: distributed surface load "elem, P<face>, magnitude" — same consistent
  // pressure-face machinery as *DLOAD P<face>. (Surface-name form deferred.)
  void dsload(const std::vector<std::string>& f, int line) {
    if (f.size() < 3) throw ParseError(line, "*DSLOAD needs element,P<face>,magnitude");
    const std::string label = upper(f[1]);
    if (label.empty() || label[0] != 'P')
      throw ParseError(line, "unsupported *DSLOAD type '" + f[1] +
                                 "' (supports pressure P<face>)");
    const int face = static_cast<int>(to_index(label.substr(1), line));
    const Real p = to_double(f[2], line);
    const std::string amp = param("AMPLITUDE");
    for (const Index eid : elem_refs(f[0], line))
      model_.dloads.push_back(Dload{eid, face, p, amp});
  }

  // Start a new named *SURFACE (TYPE=ELEMENT default, or TYPE=NODE).
  void begin_surface(int line) {
    const std::string name = param("NAME");
    if (name.empty()) throw ParseError(line, "*SURFACE without NAME");
    const std::string type = param("TYPE");
    Surface s;
    s.name = name;
    if (type.empty() || type == "ELEMENT") {
      s.type = Surface::Type::Element;
    } else if (type == "NODE") {
      s.type = Surface::Type::Node;
    } else {
      throw ParseError(line, "unsupported *SURFACE TYPE=" + type +
                                 " (expected ELEMENT or NODE)");
    }
    cur_surface_ = static_cast<int>(surfaces_.size());
    surfaces_.push_back(std::move(s));
  }

  // Parse a face id 'S<n>' -> n (1..4 tet, 1..5 wedge, 1..6 hex); reject anything else.
  int surface_face(const std::string& tok, int line) {
    const std::string t = upper(tok);
    if (t.size() < 2 || t[0] != 'S')
      throw ParseError(line, "*SURFACE face must be S<n>, got '" + tok + "'");
    const int face = static_cast<int>(to_index(t.substr(1), line));
    if (face < 1 || face > 6)
      throw ParseError(line, "*SURFACE face out of range (1..6): '" + tok + "'");
    return face;
  }

  void surface_data(const std::vector<std::string>& f, int line) {
    if (cur_surface_ < 0) throw ParseError(line, "*SURFACE data without a surface");
    Surface& s = surfaces_[static_cast<std::size_t>(cur_surface_)];
    if (s.type == Surface::Type::Node) {
      if (f.empty()) throw ParseError(line, "*SURFACE (NODE) needs a node id or nset");
      for (const std::string& tok : f)
        for (const Index nd : node_refs(tok, line)) s.nodes.push_back(nd);
      return;
    }
    // TYPE=ELEMENT: "elset_or_elem, S<face>".
    if (f.size() < 2)
      throw ParseError(line, "*SURFACE (ELEMENT) needs 'elset_or_elem, S<face>'");
    const int face = surface_face(f[1], line);
    for (const Index eid : elem_refs(f[0], line))
      s.faces.emplace_back(eid, face);
  }

  // *SURFACE INTERACTION, NAME=<name>: open a named interaction that a *SURFACE BEHAVIOR
  // fills and a *CONTACT PAIR references. (spec: contact — surface interaction.)
  void begin_surface_interaction(int line) {
    const std::string name = param("NAME");
    if (name.empty()) throw ParseError(line, "*SURFACE INTERACTION without NAME");
    cur_interaction_ = name;
    SurfaceInteraction si;
    si.name = name;
    model_.surface_interactions[name] = si;
  }

  // *SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=HARD|LINEAR|EXPONENTIAL. The law is set from
  // the parameter here; its parameters (slope / p0,c0) come on the data line. Attaches to
  // the most recent *SURFACE INTERACTION. (spec: contact — surface behavior normal.)
  void begin_surface_behavior(int line) {
    if (cur_interaction_.empty())
      throw ParseError(line, "*SURFACE BEHAVIOR outside a *SURFACE INTERACTION");
    SurfaceInteraction& si = model_.surface_interactions[cur_interaction_];
    std::string po = param("PRESSURE-OVERCLOSURE");
    if (po.empty()) po = param("PRESSUREOVERCLOSURE");
    if (po.empty() || po == "HARD") si.behavior.law = SurfaceBehavior::Law::Hard;
    else if (po == "LINEAR") si.behavior.law = SurfaceBehavior::Law::Linear;
    else if (po == "EXPONENTIAL") si.behavior.law = SurfaceBehavior::Law::Exponential;
    else if (po == "TIED") si.behavior.law = SurfaceBehavior::Law::Hard;  // treat tied as hard
    else throw ParseError(line, "unsupported *SURFACE BEHAVIOR PRESSURE-OVERCLOSURE=" + po);
    si.has_behavior = true;
  }

  // *SURFACE BEHAVIOR data line: LINEAR -> "slope"; EXPONENTIAL -> "c0, p0" (CalculiX
  // orders the exponential card as c0=pressure-at-zero-clearance, then the clearance);
  // HARD -> no data. (spec: contact — surface behavior parameters.)
  void surface_behavior_data(const std::vector<std::string>& f, int line) {
    if (cur_interaction_.empty()) return;
    SurfaceInteraction& si = model_.surface_interactions[cur_interaction_];
    if (si.behavior.law == SurfaceBehavior::Law::Linear) {
      if (f.empty()) throw ParseError(line, "*SURFACE BEHAVIOR LINEAR needs a slope");
      si.behavior.k = to_double(f[0], line);
    } else if (si.behavior.law == SurfaceBehavior::Law::Exponential) {
      if (f.size() < 2)
        throw ParseError(line, "*SURFACE BEHAVIOR EXPONENTIAL needs 'c0, p0'");
      // CalculiX EXPONENTIAL card: field 1 = c0 (pressure at zero clearance p0 in our
      // struct), field 2 = the clearance at which pressure vanishes (c0 in our struct).
      si.behavior.p0 = to_double(f[0], line);
      si.behavior.c0 = to_double(f[1], line);
    }
    // HARD: no parameters.
  }

  // *FRICTION: open Coulomb friction on the most recent *SURFACE INTERACTION; mu and the
  // optional stick stiffness come on the data line. (spec: contact — Coulomb friction.)
  void begin_friction(int line) {
    if (cur_interaction_.empty())
      throw ParseError(line, "*FRICTION outside a *SURFACE INTERACTION");
    model_.surface_interactions[cur_interaction_].friction.has = true;
  }

  // *FRICTION data line: "mu [, stick_stiffness]". mu is the Coulomb coefficient; the
  // optional second field is the tangential stick (lambda) stiffness (0/absent -> the
  // contact operator auto-sizes it from the normal penalty). Attaches to the most recent
  // *SURFACE INTERACTION. (spec: contact — tangential contact / stick-slip.)
  void friction_data(const std::vector<std::string>& f, int line) {
    if (cur_interaction_.empty()) return;
    Friction& fr = model_.surface_interactions[cur_interaction_].friction;
    if (f.empty()) throw ParseError(line, "*FRICTION needs a friction coefficient mu");
    fr.mu = to_double(f[0], line);
    if (f.size() >= 2 && !f[1].empty()) fr.stick = to_double(f[1], line);
    fr.has = true;
  }

  // *GAP CONDUCTANCE: open thermal contact conductance on the most recent *SURFACE
  // INTERACTION; the coefficient comes on the data line. (spec: contact — thermal contact.)
  void begin_gap_conductance(int line) {
    if (cur_interaction_.empty())
      throw ParseError(line, "*GAP CONDUCTANCE outside a *SURFACE INTERACTION");
    model_.surface_interactions[cur_interaction_].conductance.has = true;
  }

  // *GAP CONDUCTANCE data line: "h [, ...]". The first field is the conductance
  // coefficient (heat per unit area per unit temperature drop). CalculiX allows a
  // pressure/temperature table; this slice consumes the constant first coefficient.
  // (spec: contact — thermal contact conductance.)
  void gap_conductance_data(const std::vector<std::string>& f, int line) {
    if (cur_interaction_.empty()) return;
    if (f.empty() || f[0].empty())
      throw ParseError(line, "*GAP CONDUCTANCE needs a conductance coefficient");
    GapConductance& gc = model_.surface_interactions[cur_interaction_].conductance;
    gc.h = to_double(f[0], line);
    gc.has = true;
  }

  // *GAP HEAT GENERATION: open interface heat generation on the most recent *SURFACE
  // INTERACTION; the generated flux comes on the data line. (spec: contact — gap heat.)
  void begin_gap_heat_generation(int line) {
    if (cur_interaction_.empty())
      throw ParseError(line, "*GAP HEAT GENERATION outside a *SURFACE INTERACTION");
    model_.surface_interactions[cur_interaction_].heat_generation.has = true;
  }

  // *GAP HEAT GENERATION data line: "q [, ...]". The first field is the generated heat
  // flux per unit contacting area (total; split evenly between the two surfaces).
  // (spec: contact — gap heat generation.)
  void gap_heat_generation_data(const std::vector<std::string>& f, int line) {
    if (cur_interaction_.empty()) return;
    if (f.empty() || f[0].empty())
      throw ParseError(line, "*GAP HEAT GENERATION needs a generated-flux value");
    GapHeatGeneration& gh = model_.surface_interactions[cur_interaction_].heat_generation;
    gh.q = to_double(f[0], line);
    gh.has = true;
  }

  // *CONTACT PAIR, INTERACTION=<name>, TYPE=NODE TO SURFACE|SURFACE TO SURFACE. The pair
  // is opened here; slave/master surfaces come on the data line. (spec: contact — pairs.)
  void begin_contact_pair(int line) {
    ContactPair cp;
    cp.interaction = param("INTERACTION");
    if (cp.interaction.empty())
      throw ParseError(line, "*CONTACT PAIR needs INTERACTION=");
    std::string type = normalize(param("TYPE"));  // strip spaces: "NODE TO SURFACE"
    if (type.empty() || type == "NODETOSURFACE")
      cp.formulation = ContactPair::Formulation::NodeToSurface;
    else if (type == "SURFACETOSURFACE")
      cp.formulation = ContactPair::Formulation::SurfaceToSurface;
    else
      throw ParseError(line, "*CONTACT PAIR TYPE must be NODE TO SURFACE or SURFACE TO "
                             "SURFACE, got '" + type + "'");
    if (!param("ADJUST").empty()) cp.search = to_double(param("ADJUST"), line);
    cur_contact_pair_ = static_cast<int>(model_.contact_pairs.size());
    model_.contact_pairs.push_back(cp);
  }

  // *CONTACT PAIR data line: "slave_surface, master_surface". (spec: contact — pairs.)
  void contact_pair_data(const std::vector<std::string>& f, int line) {
    if (cur_contact_pair_ < 0)
      throw ParseError(line, "*CONTACT PAIR data without a pair");
    if (f.size() < 2)
      throw ParseError(line, "*CONTACT PAIR needs 'slave_surface, master_surface'");
    ContactPair& cp = model_.contact_pairs[static_cast<std::size_t>(cur_contact_pair_)];
    // Surface names are stored uppercased (like *SURFACE NAME=); match that here so a
    // mixed-case reference on the data line resolves. (Mesh surface keys are uppercase.)
    cp.slave_surface = upper(f[0]);
    cp.master_surface = upper(f[1]);
  }

  // *CLEARANCE, MASTER=<surf>, SLAVE=<surf>, VALUE=<v>: set a uniform initial clearance on
  // the matching *CONTACT PAIR (spec: contact — contact modifiers). The clearance overrides
  // the geometric gap so the interface starts at exactly `VALUE` regardless of the mesh.
  // The card carries no data line. (spec: contact — *CLEARANCE initial gap.)
  void begin_clearance(int line) {
    const std::string master = upper(param("MASTER"));
    const std::string slave = upper(param("SLAVE"));
    if (master.empty() || slave.empty())
      throw ParseError(line, "*CLEARANCE needs MASTER= and SLAVE= surface names");
    if (param("VALUE").empty())
      throw ParseError(line, "*CLEARANCE needs VALUE=");
    const Real value = to_double(param("VALUE"), line);
    bool matched = false;
    for (ContactPair& cp : model_.contact_pairs) {
      const bool same = (cp.master_surface == master && cp.slave_surface == slave) ||
                        (cp.master_surface == slave && cp.slave_surface == master);
      if (same) {
        cp.clearance = value;
        cp.has_clearance = true;
        matched = true;
      }
    }
    if (!matched)
      throw ParseError(line, "*CLEARANCE MASTER=/SLAVE= match no *CONTACT PAIR");
  }

  void flush_sets() { flush_sets_into(model_); }

  // Add the parsed nsets/elsets/surfaces into `m`'s mesh. Used for the single-model
  // path (into model_) and for each step model in run_steps. Surfaces are copied (not
  // moved) so several step meshes can each receive them.
  void flush_sets_into(Model& m) {
    for (auto& [name, ids] : nsets_) m.mesh.add_nset(name, ids);
    for (auto& [name, ids] : elsets_) m.mesh.add_elset(name, ids);
    for (const auto& s : surfaces_) m.mesh.add_surface(s);
  }
};

}  // namespace

Model parse_inp(const std::string& text) { return Parser{}.run(text); }

std::vector<Model> parse_inp_steps(const std::string& text) {
  return Parser{}.run_steps(text);
}

Model parse_inp_file(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ParseError(0, "cannot open file '" + path + "'");
  std::stringstream ss;
  ss << in.rdbuf();
  return parse_inp(ss.str());
}

std::vector<Model> parse_inp_steps_file(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ParseError(0, "cannot open file '" + path + "'");
  std::stringstream ss;
  ss << in.rdbuf();
  return parse_inp_steps(ss.str());
}

void validate_preceding_frequency(const std::vector<Model>& steps) {
  bool saw_frequency = false;
  for (std::size_t i = 0; i < steps.size(); ++i) {
    const Procedure p = steps[i].procedure;
    const bool needs_basis = p == Procedure::ModalDynamic ||
                             p == Procedure::SteadyStateDynamics ||
                             p == Procedure::ComplexFrequency;
    if (needs_basis && !saw_frequency)
      throw std::runtime_error(
          "step " + std::to_string(i + 1) +
          " (*MODAL DYNAMIC / *STEADY STATE DYNAMICS / *COMPLEX FREQUENCY) requires a "
          "preceding *FREQUENCY step in the same job that produced the eigenmodes (none "
          "found)");
    if (p == Procedure::Frequency) saw_frequency = true;
  }
}

}  // namespace cxpp::io
