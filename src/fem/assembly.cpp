#include "calculixpp/fem/assembly.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>

#include "calculixpp/fem/element.hpp"
#include "calculixpp/fem/loads.hpp"

namespace cxpp::fem {
namespace {

// Set up the DOF map (single-point constraints applied, MPC slaves eliminated, free
// DOFs numbered) and the prescribed / amplitude arrays. Shared by the linear and
// material-point assembly. Builds sys.transform from the model's constraints; over-
// constraints throw here (spec: constraints — dependent-DOF elimination + detection).
void build_dof_map(const Model& model, LinearSystem& sys) {
  const Mesh& mesh = model.mesh;
  const std::size_t n_dof = mesh.num_nodes() * kDofsPerNode;
  sys.prescribed.assign(n_dof, 0.0);
  sys.prescribed_amp.assign(n_dof, std::string{});
  std::vector<bool> constrained(n_dof, false);

  for (const Spc& spc : model.spcs) {
    const Index ni = mesh.node_index(spc.node_id);
    if (ni < 0) throw std::runtime_error("*BOUNDARY references unknown node");
    if (spc.comp < 1 || spc.comp > kDofsPerNode)
      throw std::runtime_error("*BOUNDARY dof out of range");
    const std::size_t g = static_cast<std::size_t>(ni) * kDofsPerNode +
                          static_cast<std::size_t>(spc.comp - 1);
    constrained[g] = true;
    sys.prescribed[g] = spc.value;
    sys.prescribed_amp[g] = spc.amplitude;
  }

  sys.transform = build_constraint_transform(model, constrained, sys.prescribed);
  sys.dof_eq = sys.transform.dof_eq;
  sys.n_free = sys.transform.n_free;
  sys.rhs.assign(static_cast<std::size_t>(sys.n_free), 0.0);
}

// Gather an element's node indices and coordinates (and, when `u` is non-null, the
// current element displacements ue).
void gather(const Mesh& mesh, const Element& elem, int n,
            const std::vector<Vec3>* u, std::vector<Index>& nidx,
            std::vector<Vec3>& coords, std::vector<Vec3>& ue) {
  nidx.resize(static_cast<std::size_t>(n));
  coords.resize(static_cast<std::size_t>(n));
  if (u) ue.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(i)]);
    if (ni < 0) throw std::runtime_error("element references unknown node");
    nidx[static_cast<std::size_t>(i)] = ni;
    coords[static_cast<std::size_t>(i)] = mesh.nodes()[static_cast<std::size_t>(ni)].x;
    if (u) ue[static_cast<std::size_t>(i)] = (*u)[static_cast<std::size_t>(ni)];
  }
}

// Global DOF index of element-local dof `a` for node indices `nidx`.
inline std::size_t global_dof(const std::vector<Index>& nidx, int a) {
  return static_cast<std::size_t>(nidx[static_cast<std::size_t>(a / kDofsPerNode)]) *
             kDofsPerNode +
         static_cast<std::size_t>(a % kDofsPerNode);
}

// Scatter a dense symmetric block into the reduced free/free COO map, distributing
// each DOF through the constraint transform (T^T K T). For DOFs a,b with block value
// kab and expansions Ea = {(eq, cA)} + constA, Eb = {(eq, cB)} + constB:
//   matrix:  cA * kab * cB  at (eqA, eqB)
//   rhs:    -cA * kab * constB  at eqA   (eliminated SPC/slave constant on column b)
// This is the congruence transform, so symmetry / positive-definiteness are
// preserved. `gdofs` maps local index -> global DOF. Shared by element tangents and
// spring blocks. (spec: constraints — SPD-preserving dependent-DOF elimination.)
void scatter_block(const std::vector<Real>& K, int m,
                   const std::vector<std::size_t>& gdofs, const LinearSystem& sys,
                   Index n_free, std::unordered_map<std::int64_t, Real>& kmap,
                   std::vector<Real>& rhs) {
  const auto key = [n_free](Index r, Index c) {
    return static_cast<std::int64_t>(r) * n_free + c;
  };
  const ConstraintTransform& tf = sys.transform;
  for (int a = 0; a < m; ++a) {
    const DofExpansion& Ea = tf.expansion[gdofs[static_cast<std::size_t>(a)]];
    if (Ea.terms.empty() && Ea.constant == 0.0) continue;
    for (int b = 0; b < m; ++b) {
      const Real kab = K[static_cast<std::size_t>(a) * static_cast<std::size_t>(m) +
                         static_cast<std::size_t>(b)];
      if (kab == 0.0) continue;
      const DofExpansion& Eb = tf.expansion[gdofs[static_cast<std::size_t>(b)]];
      for (const DofTerm& ta : Ea.terms) {
        for (const DofTerm& tb : Eb.terms)
          kmap[key(ta.eq, tb.eq)] += ta.coeff * kab * tb.coeff;
        if (Eb.constant != 0.0)
          rhs[static_cast<std::size_t>(ta.eq)] -= ta.coeff * kab * Eb.constant;
      }
    }
  }
}

// Element-tangent adapter: build the global-DOF list for the element and scatter.
void scatter_tangent(const std::vector<Real>& Ke, int ndof,
                     const std::vector<Index>& nidx, const LinearSystem& sys,
                     Index n_free, std::unordered_map<std::int64_t, Real>& kmap,
                     std::vector<Real>& rhs) {
  std::vector<std::size_t> gdofs(static_cast<std::size_t>(ndof));
  for (int a = 0; a < ndof; ++a) gdofs[static_cast<std::size_t>(a)] = global_dof(nidx, a);
  scatter_block(Ke, ndof, gdofs, sys, n_free, kmap, rhs);
}

// Global DOF (node_index*3 + comp) for a spring endpoint (node id, comp 1..3).
std::size_t spring_global_dof(const Mesh& mesh, Index node_id, int comp) {
  const Index ni = mesh.node_index(node_id);
  if (ni < 0) throw std::runtime_error("*SPRING references unknown node");
  if (comp < 1 || comp > kDofsPerNode)
    throw std::runtime_error("*SPRING dof out of range (1..3)");
  return static_cast<std::size_t>(ni) * kDofsPerNode + static_cast<std::size_t>(comp - 1);
}

// Scatter a spring's local stiffness (given as a list of global DOFs `gdof` and a
// dense symmetric block `kloc`) into the reduced free/free COO map, moving
// constrained columns to the RHS exactly like scatter_tangent. Shared by all three
// spring kinds. (spec: element-sections — *SPRING stiffness contribution.)
// `u` and `f_int` (both full nodal DOF length) are optional: when non-null, the
// spring's internal force K_loc * u_loc is scattered into f_int (used by the
// nonlinear material-point residual; springs are linear so K u is exact).
void scatter_spring_block(const std::vector<std::size_t>& gdof,
                          const std::vector<Real>& kloc, const LinearSystem& sys,
                          Index n_free, std::unordered_map<std::int64_t, Real>& kmap,
                          std::vector<Real>& rhs, const std::vector<Vec3>* u,
                          std::vector<Real>* f_int) {
  const int m = static_cast<int>(gdof.size());
  // Internal force K_loc * u_loc into the full-DOF f_int (reduced through the
  // transform later by the driver). Springs are linear so K u is exact.
  if (u && f_int) {
    const auto dof_val = [&](std::size_t g) {
      return (*u)[g / kDofsPerNode][g % kDofsPerNode];
    };
    for (int a = 0; a < m; ++a)
      for (int b = 0; b < m; ++b)
        (*f_int)[gdof[static_cast<std::size_t>(a)]] +=
            kloc[static_cast<std::size_t>(a) * static_cast<std::size_t>(m) +
                 static_cast<std::size_t>(b)] * dof_val(gdof[static_cast<std::size_t>(b)]);
  }
  // Tangent through the constraint transform (same congruence path as elements).
  scatter_block(kloc, m, gdof, sys, n_free, kmap, rhs);
}

// Assemble all *SPRING connectors into the reduced tangent. SPRING1 (grounded) is a
// scalar k on one DOF; SPRING2 (dof) couples two DOFs with [[k,-k],[-k,k]]; SPRINGA
// (axial) builds the 6x6 block k*[[n nᵀ,-n nᵀ],[-n nᵀ,n nᵀ]] along the node line.
void scatter_springs(const Model& model, const LinearSystem& sys, Index n_free,
                     std::unordered_map<std::int64_t, Real>& kmap,
                     std::vector<Real>& rhs, const std::vector<Vec3>* u = nullptr,
                     std::vector<Real>* f_int = nullptr) {
  const Mesh& mesh = model.mesh;
  for (const Spring& sp : model.springs) {
    const Real k = sp.stiffness;
    if (sp.kind == Spring::Kind::Grounded) {
      const std::vector<std::size_t> g{spring_global_dof(mesh, sp.node1, sp.dof1)};
      scatter_spring_block(g, {k}, sys, n_free, kmap, rhs, u, f_int);
    } else if (sp.kind == Spring::Kind::Dof) {
      const std::vector<std::size_t> g{spring_global_dof(mesh, sp.node1, sp.dof1),
                                       spring_global_dof(mesh, sp.node2, sp.dof2)};
      scatter_spring_block(g, {k, -k, -k, k}, sys, n_free, kmap, rhs, u, f_int);
    } else {  // Axial (SPRINGA): direction from node1 -> node2.
      const Index i1 = mesh.node_index(sp.node1), i2 = mesh.node_index(sp.node2);
      if (i1 < 0 || i2 < 0) throw std::runtime_error("*SPRING references unknown node");
      const Vec3 x1 = mesh.nodes()[static_cast<std::size_t>(i1)].x;
      const Vec3 x2 = mesh.nodes()[static_cast<std::size_t>(i2)].x;
      std::array<Real, 3> n{x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
      const Real len = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
      if (len <= 0.0) throw std::runtime_error("SPRINGA has coincident nodes");
      for (Real& c : n) c /= len;
      // 6 global DOFs: node1 x,y,z then node2 x,y,z.
      std::vector<std::size_t> g(6);
      for (int c = 0; c < 3; ++c) {
        g[static_cast<std::size_t>(c)] = spring_global_dof(mesh, sp.node1, c + 1);
        g[static_cast<std::size_t>(3 + c)] = spring_global_dof(mesh, sp.node2, c + 1);
      }
      // kloc[a][b] = k * s_a s_b n_(a%3) n_(b%3), s = +1 for node1 block, -1 node2.
      std::vector<Real> kloc(36, 0.0);
      for (int a = 0; a < 6; ++a)
        for (int b = 0; b < 6; ++b) {
          const Real sa = a < 3 ? 1.0 : -1.0, sbb = b < 3 ? 1.0 : -1.0;
          kloc[static_cast<std::size_t>(a * 6 + b)] =
              k * sa * sbb * n[static_cast<std::size_t>(a % 3)] * n[static_cast<std::size_t>(b % 3)];
        }
      scatter_spring_block(g, kloc, sys, n_free, kmap, rhs, u, f_int);
    }
  }
}

// Add a full-DOF external load vector into the reduced rhs through the constraint
// transform: a load on a free DOF lands on its equation; a load on an MPC slave DOF
// is distributed to its masters by the equation coefficients (T^T f). Loads on SPC
// DOFs go to the reactions, not the free rhs. (spec: constraints — consistent load
// transfer through the elimination transform.)
void add_external_load(const LinearSystem& sys, const std::vector<Real>& fext,
                       std::vector<Real>& rhs) {
  const ConstraintTransform& tf = sys.transform;
  for (std::size_t g = 0; g < fext.size(); ++g) {
    if (fext[g] == 0.0) continue;
    for (const DofTerm& t : tf.expansion[g].terms)
      rhs[static_cast<std::size_t>(t.eq)] += t.coeff * fext[g];
  }
}

// Flush the accumulated free/free COO map into the LinearSystem triplet arrays.
void flush_coo(const std::unordered_map<std::int64_t, Real>& kmap, Index n_free,
               LinearSystem& sys) {
  sys.rows.reserve(kmap.size());
  sys.cols.reserve(kmap.size());
  sys.vals.reserve(kmap.size());
  for (const auto& [k, v] : kmap) {
    sys.rows.push_back(static_cast<Index>(k / n_free));
    sys.cols.push_back(static_cast<Index>(k % n_free));
    sys.vals.push_back(v);
  }
}

}  // namespace

LinearSystem assemble_linear_static(const Model& model) {
  const Mesh& mesh = model.mesh;
  LinearSystem sys;
  build_dof_map(model, sys);
  const Index n_free = sys.n_free;

  const std::vector<ElasticIso> elastic = model.element_elastic();
  const std::vector<bool> active = model.element_active_mask();
  std::unordered_map<std::int64_t, Real> kmap;
  std::vector<Index> nidx;
  std::vector<Vec3> coords, ue;

  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    if (!active[e]) continue;  // *MODEL CHANGE, REMOVE: element carries no stiffness
    const Element& elem = mesh.elements()[e];
    const int n = nodes_per_element(elem.type);
    gather(mesh, elem, n, nullptr, nidx, coords, ue);
    const std::vector<Real> Ke = element_stiffness(elem.type, coords, elastic[e]);
    scatter_tangent(Ke, n * kDofsPerNode, nidx, sys, n_free, kmap, sys.rhs);
  }

  // Discrete spring connectors add their stiffness directly to the reduced system.
  scatter_springs(model, sys, n_free, kmap, sys.rhs);

  // External loads, distributed through the constraint transform.
  add_external_load(sys, external_load_vector(model), sys.rhs);
  // Thermal-strain equivalent load (zero unless *EXPANSION + an applied temperature
  // field are present, so the pure-mechanical path is untouched).
  add_external_load(sys, thermal_strain_load_vector(model), sys.rhs);

  flush_coo(kmap, n_free, sys);
  return sys;
}

LinearSystem assemble_mass(const Model& model, bool lumped) {
  const Mesh& mesh = model.mesh;
  LinearSystem sys;
  build_dof_map(model, sys);
  const Index n_free = sys.n_free;

  const std::vector<Real> density = model.element_density();
  const std::vector<bool> active = model.element_active_mask();
  std::unordered_map<std::int64_t, Real> mmap;
  std::vector<Index> nidx;
  std::vector<Vec3> coords, ue;
  std::vector<Real> dummy_rhs(static_cast<std::size_t>(n_free), 0.0);

  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    if (!active[e]) continue;  // *MODEL CHANGE, REMOVE: element carries no mass
    const Element& elem = mesh.elements()[e];
    const int n = nodes_per_element(elem.type);
    gather(mesh, elem, n, nullptr, nidx, coords, ue);
    const std::vector<Real> Me =
        lumped ? element_mass_lumped(elem.type, coords, density[e])
               : element_mass(elem.type, coords, density[e]);
    // Scatter through the same congruence transform as the stiffness; the RHS side
    // effect (eliminated-constant columns) is discarded — mass never drives a load.
    scatter_tangent(Me, n * kDofsPerNode, nidx, sys, n_free, mmap, dummy_rhs);
  }

  // *MASS point-mass connectors: a diagonal nodal mass on the 3 translational DOFs.
  for (const PointMass& pm : model.point_masses) {
    if (pm.mass == 0.0) continue;
    for (int c = 0; c < kDofsPerNode; ++c) {
      const std::size_t g = spring_global_dof(mesh, pm.node, c + 1);
      const std::vector<std::size_t> gd{g};
      scatter_block({pm.mass}, 1, gd, sys, n_free, mmap, dummy_rhs);
    }
  }

  flush_coo(mmap, n_free, sys);
  return sys;
}

// Build the MaterialModel for element `e` from the resolved per-element material data,
// dispatching by constitutive kind (user material > plastic > hyperelastic > elastic).
// `ndepvar` is set to the *DEPVAR count so the caller can size the state vectors.
static std::unique_ptr<MaterialModel> make_element_model(
    const ElasticIso& elastic, const std::optional<Plastic>& plastic,
    const std::optional<Hyperelastic>& hyper,
    const std::optional<UserMaterial>& user, int& ndepvar) {
  ndepvar = 0;
  if (user && !user->empty()) {
    ndepvar = user->ndepvar;
    const UserMaterialFactory* factory =
        UserMaterialRegistry::instance().find(user->name);
    if (factory == nullptr)
      throw std::runtime_error("*USER MATERIAL '" + user->name +
                               "' is not registered (no factory)");
    return (*factory)(*user, elastic);
  }
  if (plastic && !plastic->empty())
    return std::make_unique<J2PlasticMaterial>(elastic, *plastic);
  if (hyper && !hyper->empty())
    return std::make_unique<HyperelasticNeoHookean>(*hyper);
  return std::make_unique<ElasticIsoMaterial>(elastic);
}

MaterialPoints make_material_points(const Model& model) {
  const Mesh& mesh = model.mesh;
  const std::vector<ElasticIso> elastic = model.element_elastic();
  const std::vector<std::optional<Plastic>> plastic = model.element_plastic();
  const std::vector<std::optional<Hyperelastic>> hyper =
      model.element_hyperelastic();
  const std::vector<std::optional<UserMaterial>> user =
      model.element_user_material();
  MaterialPoints mp;
  mp.models.reserve(mesh.num_elements());
  mp.state.reserve(mesh.num_elements());
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    const Element& elem = mesh.elements()[e];
    int ndepvar = 0;
    mp.models.push_back(
        make_element_model(elastic[e], plastic[e], hyper[e], user[e], ndepvar));
    const std::size_t ngp = gauss_rule(elem.type).size();
    std::vector<MaterialState> gp(ngp);
    if (ndepvar > 0)
      for (MaterialState& s : gp) {
        s.committed_depvar.assign(static_cast<std::size_t>(ndepvar), 0.0);
        s.depvar.assign(static_cast<std::size_t>(ndepvar), 0.0);
      }
    mp.state.push_back(std::move(gp));
  }
  return mp;
}

LinearSystem assemble_material_tangent(const Model& model,
                                       const std::vector<Vec3>& u,
                                       MaterialPoints& mp,
                                       std::vector<Real>& f_int) {
  const Mesh& mesh = model.mesh;
  LinearSystem sys;
  build_dof_map(model, sys);
  const Index n_free = sys.n_free;

  f_int.assign(mesh.num_nodes() * kDofsPerNode, 0.0);
  const std::vector<bool> active = model.element_active_mask();
  std::unordered_map<std::int64_t, Real> kmap;
  std::vector<Index> nidx;
  std::vector<Vec3> coords, ue;

  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    if (!active[e]) continue;  // *MODEL CHANGE, REMOVE: no tangent / internal force
    const Element& elem = mesh.elements()[e];
    const int n = nodes_per_element(elem.type);
    const int ndof = n * kDofsPerNode;
    gather(mesh, elem, n, &u, nidx, coords, ue);
    const ElementResponse resp = element_tangent_force(
        elem.type, coords, ue, *mp.models[e], mp.state[e]);
    scatter_tangent(resp.Ke, ndof, nidx, sys, n_free, kmap, sys.rhs);
    for (int a = 0; a < ndof; ++a)
      f_int[global_dof(nidx, a)] += resp.fe[static_cast<std::size_t>(a)];
  }

  // Spring connectors: add tangent to the reduced system and K_spring*u to f_int.
  scatter_springs(model, sys, n_free, kmap, sys.rhs, &u, &f_int);

  // External loads at full magnitude, distributed through the constraint transform
  // (driver overwrites rhs with its residual; this keeps rhs meaningful for callers
  // that want K u = f directly).
  add_external_load(sys, external_load_vector(model), sys.rhs);

  flush_coo(kmap, n_free, sys);
  return sys;
}

}  // namespace cxpp::fem
