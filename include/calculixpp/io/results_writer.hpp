#pragma once
#include <string>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/results.hpp"

// Result writers (spec: results-output). `.frd` targets the CalculiX/CGX format
// (nodes, elements, DISP/STRESS/FORC datasets); `.dat` is a plain tabular dump.
namespace cxpp::io {

void write_frd(const std::string& path, const Model& model, const StaticFields& f);
void write_dat(const std::string& path, const Model& model, const StaticFields& f);

// Thermal (heat-transfer) result writers (spec: results-output — NT field). `.frd`
// writes the NDTEMP dataset (CGX temperature field); `.dat` writes the nodal
// temperature and heat-flux-reaction tables matching CalculiX's *NODE PRINT NT/RFL.
void write_frd(const std::string& path, const Model& model, const ThermalFields& t);
void write_dat(const std::string& path, const Model& model, const ThermalFields& t);

}  // namespace cxpp::io
