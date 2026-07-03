#include "calculixpp/fem/constraint_transform.hpp"

#include <stdexcept>
#include <unordered_map>

#include "calculixpp/core/types.hpp"

namespace cxpp::fem {
namespace {

// Global DOF index for (node index, comp 0..2).
inline std::size_t gdof(Index ni, int comp) {
  return static_cast<std::size_t>(ni) * kDofsPerNode + static_cast<std::size_t>(comp);
}

// The slave DOF's defining relation, stored as u_slave = constant + sum (coeff*master).
struct SlaveDef {
  std::size_t dof{};
  std::vector<std::pair<std::size_t, Real>> masters;  // (global master dof, coeff)
  std::string origin;
};

}  // namespace

ConstraintTransform build_constraint_transform(
    const Model& model, const std::vector<bool>& spc_constrained,
    const std::vector<Real>& prescribed_value) {
  const Mesh& mesh = model.mesh;
  const std::size_t n_dof = mesh.num_nodes() * kDofsPerNode;
  ConstraintTransform tf;
  tf.dof_eq.assign(n_dof, -1);
  tf.is_slave.assign(n_dof, false);
  tf.expansion.assign(n_dof, DofExpansion{});
  tf.spc_terms.assign(n_dof, {});

  const std::vector<Equation> eqs = model.expand_constraints();

  // Resolve each equation's dependent (first) term to a SlaveDef, checking for
  // over-constraints as we go.
  std::unordered_map<std::size_t, SlaveDef> slaves;  // slave dof -> definition
  for (const Equation& eq : eqs) {
    if (eq.terms.empty()) continue;
    const EquationTerm& dep = eq.terms.front();
    const Index dni = mesh.node_index(dep.node_id);
    if (dni < 0)
      throw std::runtime_error(eq.origin + ": unknown node " + std::to_string(dep.node_id));
    if (dep.comp < 1 || dep.comp > kDofsPerNode)
      throw std::runtime_error(eq.origin + ": dependent DOF out of range");
    const std::size_t sdof = gdof(dni, dep.comp - 1);
    if (dep.coeff == 0.0)
      throw std::runtime_error(eq.origin + ": dependent term has zero coefficient");
    if (spc_constrained[sdof])
      throw std::runtime_error(
          "over-constraint: dependent DOF of " + eq.origin +
          " (node " + std::to_string(dep.node_id) + ") is also a *BOUNDARY SPC");
    if (slaves.count(sdof) != 0)
      throw std::runtime_error(
          "over-constraint: DOF (node " + std::to_string(dep.node_id) +
          ") is the dependent term of more than one constraint (" +
          slaves.at(sdof).origin + " and " + eq.origin + ")");

    SlaveDef def;
    def.dof = sdof;
    def.origin = eq.origin;
    for (std::size_t k = 1; k < eq.terms.size(); ++k) {
      const EquationTerm& t = eq.terms[k];
      const Index ni = mesh.node_index(t.node_id);
      if (ni < 0)
        throw std::runtime_error(eq.origin + ": unknown node " + std::to_string(t.node_id));
      if (t.comp < 1 || t.comp > kDofsPerNode)
        throw std::runtime_error(eq.origin + ": DOF out of range");
      def.masters.emplace_back(gdof(ni, t.comp - 1), -t.coeff / dep.coeff);
    }
    slaves.emplace(sdof, std::move(def));
    tf.is_slave[sdof] = true;
  }

  // Number the free DOFs: everything that is neither SPC nor slave.
  Index n_free = 0;
  for (std::size_t g = 0; g < n_dof; ++g)
    if (!spc_constrained[g] && !tf.is_slave[g]) tf.dof_eq[g] = n_free++;
  tf.n_free = n_free;

  // Recursively expand a global DOF into free-equation contributions (+ constant from
  // the full-magnitude prescribed values) and, separately, the SPC-master terms
  // (global spc dof, weight) so the driver can rebuild slave motion under a scaled
  // load factor. `visiting` guards against a dependency cycle.
  struct Expanded {
    DofExpansion free;                 // free-equation terms + full-magnitude constant
    std::vector<DofTerm> spc;          // SPC masters as (global spc dof, weight)
  };
  std::unordered_map<std::size_t, Expanded> memo;
  std::vector<bool> visiting(n_dof, false);
  auto expand = [&](std::size_t g, auto&& self) -> Expanded {
    if (tf.dof_eq[g] >= 0)
      return Expanded{DofExpansion{{{tf.dof_eq[g], 1.0}}, 0.0}, {}};
    if (!tf.is_slave[g])  // SPC: constant (full magnitude) + one SPC master term
      return Expanded{DofExpansion{{}, prescribed_value[g]},
                      {{static_cast<Index>(g), 1.0}}};
    const auto it = memo.find(g);
    if (it != memo.end()) return it->second;
    if (visiting[g])
      throw std::runtime_error("over-constraint: cyclic constraint dependency at DOF " +
                               std::to_string(g));
    visiting[g] = true;
    Expanded acc;
    std::unordered_map<Index, Real> term_by_eq;
    std::unordered_map<Index, Real> spc_by_dof;
    const SlaveDef& def = slaves.at(g);
    for (const auto& [mdof, coeff] : def.masters) {
      const Expanded me = self(mdof, self);
      acc.free.constant += coeff * me.free.constant;
      for (const DofTerm& dt : me.free.terms) term_by_eq[dt.eq] += coeff * dt.coeff;
      for (const DofTerm& st : me.spc) spc_by_dof[st.eq] += coeff * st.coeff;
    }
    for (const auto& [eq, c] : term_by_eq)
      if (c != 0.0) acc.free.terms.push_back({eq, c});
    for (const auto& [dof, c] : spc_by_dof)
      if (c != 0.0) acc.spc.push_back({dof, c});
    visiting[g] = false;
    memo.emplace(g, acc);
    return acc;
  };

  for (std::size_t g = 0; g < n_dof; ++g) {
    Expanded e = expand(g, expand);
    tf.expansion[g] = std::move(e.free);
    tf.spc_terms[g] = std::move(e.spc);
  }
  return tf;
}

}  // namespace cxpp::fem
