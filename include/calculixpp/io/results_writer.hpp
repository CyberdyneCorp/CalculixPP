#pragma once
#include <string>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/results.hpp"

// Result writers (spec: results-output). `.frd` targets the CalculiX/CGX format
// (nodes, elements, DISP/STRESS/FORC datasets); `.dat` is a plain tabular dump.
namespace cxpp::io {

void write_frd(const std::string& path, const Model& model, const StaticFields& f);
void write_dat(const std::string& path, const Model& model, const StaticFields& f);

}  // namespace cxpp::io
