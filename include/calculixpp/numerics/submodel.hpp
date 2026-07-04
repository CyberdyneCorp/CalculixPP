#pragma once
#include <optional>
#include <span>
#include <vector>

#include "calculixpp/compute/backend.hpp"
#include "calculixpp/core/model.hpp"
#include "calculixpp/core/results.hpp"
#include "calculixpp/core/submodel.hpp"
#include "calculixpp/core/types.hpp"

// Submodeling driver (spec: submodeling). Displacement-driven slice: a finer local
// model whose cut-boundary DOFs are prescribed to the global displacement field
// interpolated at each boundary node's location, then solved as an ordinary linear
// static step. The interpolation locates the host global element containing the node
// and evaluates the global nodal displacement with that element's shape functions —
// reusing the fem element machinery, no new element math. (ref: src/submodels.f,
// src/interpolsubmodel.f, src/attach_3d.f.)
namespace cxpp::numerics {

// The natural coordinates of a physical point inside a global element, with the
// residual distance |x(ξ) - X| of the converged inverse map. `inside` is true when ξ
// lies within the element's reference domain (within tolerance).
struct NaturalCoords {
  Real xi{0}, et{0}, ze{0};
  Real distance{0};   // |x(ξ) - X| at the converged ξ
  Real overshoot{0};  // how far ξ lies outside the reference domain (0 when inside)
  bool inside{false};
};

// Newton inverse-isoparametric map: find the natural coordinates ξ of physical point
// `X` inside a global element of the given `type` with node coordinates `coords`
// (element node order). Iterates x(ξ) = Σ N_a(ξ) x_a toward X, solving J dξ = X - x(ξ)
// with J = dx/dξ each step. Returns the converged ξ, the residual distance, and whether
// ξ is inside the reference domain. Robust for a point outside the element (distance is
// then the map's closest attainable image).
NaturalCoords natural_coords(ElementType type, std::span<const Vec3> coords,
                             const Vec3& X);

// Interpolate the global displacement field at physical point `X`: locate the host
// global element (the one containing X, restricted to `global_elset_mask` when
// non-empty) and return Σ_a N_a(ξ) U_a over its nodes. `global_elset_mask`, when
// non-empty, is a per-global-element bool selecting the searched subset (a *SUBMODEL's
// global element set). Throws std::runtime_error when no host element contains X.
Vec3 interpolate_global_displacement(const GlobalSolution& global, const Vec3& X,
                                     const std::vector<bool>& global_elset_mask = {});

// Solve a submodel: fill every DRIVEN *BOUNDARY, SUBMODEL DOF with the global
// displacement interpolated at the boundary node's undeformed location, then run the
// ordinary linear-static solve. The global element set of the first *SUBMODEL card (or
// all global elements when it names none) bounds the host-element search. Throws when a
// driven node lies outside the global element set.
StaticFields solve_submodel(const Model& model, const GlobalSolution& global,
                            std::optional<compute::SolverKind> forced = std::nullopt);

}  // namespace cxpp::numerics
