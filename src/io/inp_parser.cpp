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
  std::unordered_map<std::string, std::vector<Index>> nsets_, elsets_;
  std::vector<Surface> surfaces_;  // built incrementally; flushed at end
  int cur_surface_ = -1;           // index into surfaces_, or -1 outside a surface

  std::string param(const std::string& key) const {
    const auto it = params_.find(key);
    return it == params_.end() ? std::string{} : it->second;
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
        "*OUTPUT",    "*RESTART",  "*STEP",     "*STATIC", "*ENDSTEP",
        "*HEADING",   "*AMPLITUDE"};
    if (card_ == "*STATIC") {
      select_solver(line);
    } else if (card_ == "*MATERIAL") {
      cur_material_ = param("NAME");
      model_.materials[cur_material_].name = cur_material_;
    } else if (card_ == "*SURFACE") {
      begin_surface(line);
    } else if (card_ == "*SOLIDSECTION") {
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

  bool is_data_card() const {
    static const std::vector<std::string> data_cards = {
        "*NODE",    "*ELEMENT",     "*NSET",     "*ELSET", "*MATERIAL",
        "*ELASTIC", "*DENSITY",     "*SOLIDSECTION", "*BOUNDARY", "*CLOAD",
        "*DLOAD",   "*SURFACE"};
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
    if (card_ == "*SURFACE") return surface_data(f, line);
    // *SOLID SECTION optional thickness line, output requests, step/static
    // control lines, amplitudes: accepted and ignored in Phase 1.
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
    for (const Index nd : node_refs(f[0], line))
      for (int c = d1; c <= d2 && c <= kDofsPerNode; ++c)
        model_.spcs.push_back(Spc{nd, c, val});
  }

  void cload(const std::vector<std::string>& f, int line) {
    if (f.size() < 3) throw ParseError(line, "*CLOAD needs node,dof,value");
    const int dof = static_cast<int>(to_index(f[1], line));
    const Real val = to_double(f[2], line);
    for (const Index nd : node_refs(f[0], line))
      model_.cloads.push_back(Cload{nd, dof, val});
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

  void dload(const std::vector<std::string>& f, int line) {
    if (f.size() < 3) throw ParseError(line, "*DLOAD needs element,P<face>,magnitude");
    const std::string label = upper(f[1]);
    if (label.empty() || label[0] != 'P')
      throw ParseError(line, "unsupported *DLOAD type '" + f[1] +
                                 "' (Phase 1 supports pressure P<face> only)");
    const int face = static_cast<int>(to_index(label.substr(1), line));
    const Real p = to_double(f[2], line);
    for (const Index eid : elem_refs(f[0], line))
      model_.dloads.push_back(Dload{eid, face, p});
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
