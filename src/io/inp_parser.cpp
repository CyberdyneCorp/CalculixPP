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
    if (t.empty() || t.rfind("**", 0) == 0) {
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

 private:
  Model model_;
  std::string card_;                                        // normalized keyword
  std::unordered_map<std::string, std::string> params_;     // normalized key -> value
  std::string cur_material_;
  std::string cur_amplitude_;
  // OP=NEW on a load/BC card resets prior loads of that type; per CalculiX it
  // takes effect on the FIRST card of that type in the step only. These flags
  // record that the first-card reset has already been consumed for the step.
  // (Single-step model: the reset clears loads accumulated earlier in the same
  // deck. Cross-step OP semantics need multi-step step handling — deferred.)
  bool op_seen_cload_ = false;
  bool op_seen_dload_ = false;
  bool op_seen_boundary_ = false;
  std::unordered_map<std::string, std::vector<Index>> nsets_, elsets_;
  std::vector<Surface> surfaces_;  // built incrementally; flushed at end
  int cur_surface_ = -1;           // index into surfaces_, or -1 outside a surface

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
        "*NODEPRINT", "*ELPRINT",  "*NODEFILE", "*ELFILE", "*CONTACTFILE",
        "*OUTPUT",    "*RESTART",  "*STEP",     "*ENDSTEP",
        "*HEADING"};
    if (card_ == "*STATIC") {
      select_solver(line);
      begin_static();
    } else if (card_ == "*CONTROLS") {
      // Parameters carried on the card; data lines fill the tolerances/limits.
    } else if (card_ == "*TIMEPOINTS") {
      // Data lines list the time points.
    } else if (card_ == "*MATERIAL") {
      cur_material_ = param("NAME");
      model_.materials[cur_material_].name = cur_material_;
    } else if (card_ == "*AMPLITUDE") {
      begin_amplitude(line);
    } else if (card_ == "*SURFACE") {
      begin_surface(line);
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
      if (apply_op_reset(op_seen_boundary_, line)) model_.spcs.clear();
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
      throw ParseError(line, "unsupported card '" + fields.front() + "'");
    }
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

  // *STATIC card: DIRECT parameter fixes the increment (no auto resizing).
  void begin_static() {
    if (params_.count("DIRECT") != 0) model_.increment.direct = true;
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

  bool is_data_card() const {
    static const std::vector<std::string> data_cards = {
        "*NODE",    "*ELEMENT",  "*NSET",         "*ELSET",     "*MATERIAL",
        "*ELASTIC", "*DENSITY",  "*SOLIDSECTION", "*BOUNDARY",  "*CLOAD",
        "*DLOAD",   "*DSLOAD",   "*SURFACE",      "*STATIC",    "*CONTROLS",
        "*TIMEPOINTS", "*AMPLITUDE", "*CHANGEPLASTIC"};
    return std::find(data_cards.begin(), data_cards.end(), card_) != data_cards.end();
  }

  void data_line(const std::vector<std::string>& f, int line) {
    if (card_ == "*NODE") return node(f, line);
    if (card_ == "*ELEMENT") return element(f, line);
    if (card_ == "*NSET") return set_data(f, line, nsets_, param("NSET"));
    if (card_ == "*ELSET") return set_data(f, line, elsets_, param("ELSET"));
    if (card_ == "*ELASTIC") return elastic(f, line);
    if (card_ == "*DENSITY") return density(f, line);
    if (card_ == "*BOUNDARY") return boundary(f, line);
    if (card_ == "*CLOAD") return cload(f, line);
    if (card_ == "*DLOAD") return dload(f, line);
    if (card_ == "*DSLOAD") return dsload(f, line);
    if (card_ == "*SURFACE") return surface_data(f, line);
    if (card_ == "*STATIC") return static_data(f, line);
    if (card_ == "*CONTROLS") return controls_data(f, line);
    if (card_ == "*TIMEPOINTS") return time_points_data(f, line);
    if (card_ == "*AMPLITUDE") return amplitude_data(f, line);
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

  void element(const std::vector<std::string>& f, int line) {
    ElementType type{};
    if (!parse_element_type(param("TYPE"), type))
      throw ParseError(line, "unsupported element TYPE=" + param("TYPE"));
    const Index id = to_index(f.at(0), line);
    std::vector<Index> conn;
    for (std::size_t k = 1; k < f.size(); ++k) conn.push_back(to_index(f[k], line));
    if (static_cast<int>(conn.size()) != nodes_per_element(type))
      throw ParseError(line, "element " + std::to_string(id) + ": expected " +
                                 std::to_string(nodes_per_element(type)) +
                                 " nodes, got " + std::to_string(conn.size()));
    model_.mesh.add_element(id, type, conn);
    if (!param("ELSET").empty()) elsets_[param("ELSET")].push_back(id);
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
    for (const Index nd : node_refs(f[0], line))
      for (int c = d1; c <= d2 && c <= kDofsPerNode; ++c)
        model_.spcs.push_back(Spc{nd, c, val, amp});
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

  // Parse a face id 'S<n>' -> n (1..4 for tets); reject anything else.
  int surface_face(const std::string& tok, int line) {
    const std::string t = upper(tok);
    if (t.size() < 2 || t[0] != 'S')
      throw ParseError(line, "*SURFACE face must be S<n>, got '" + tok + "'");
    const int face = static_cast<int>(to_index(t.substr(1), line));
    if (face < 1 || face > 4)
      throw ParseError(line, "*SURFACE face out of range (1..4): '" + tok + "'");
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

  void flush_sets() {
    for (auto& [name, ids] : nsets_) model_.mesh.add_nset(name, ids);
    for (auto& [name, ids] : elsets_) model_.mesh.add_elset(name, ids);
    for (auto& s : surfaces_) model_.mesh.add_surface(std::move(s));
  }
};

}  // namespace

Model parse_inp(const std::string& text) { return Parser{}.run(text); }

Model parse_inp_file(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ParseError(0, "cannot open file '" + path + "'");
  std::stringstream ss;
  ss << in.rdbuf();
  return parse_inp(ss.str());
}

}  // namespace cxpp::io
