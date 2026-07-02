#pragma once
#include <stdexcept>
#include <string>

#include "calculixpp/core/model.hpp"

// Abaqus-style input-deck parser (spec: input-deck-parsing — Phase-1 subset).
// Supported cards: *NODE, *ELEMENT (C3D4/C3D10), *NSET, *ELSET, *SURFACE,
// *MATERIAL, *ELASTIC, *DENSITY, *SOLID SECTION, *BOUNDARY, *CLOAD, *STEP,
// *STATIC.
// Output-request cards (*NODE PRINT, *EL PRINT, *NODE FILE, *EL FILE, ...) are
// accepted and ignored. Unsupported cards raise ParseError.
namespace cxpp::io {

struct ParseError : std::runtime_error {
  ParseError(int line, const std::string& msg)
      : std::runtime_error("input deck line " + std::to_string(line) + ": " + msg),
        line(line) {}
  int line;
};

Model parse_inp(const std::string& text);
Model parse_inp_file(const std::string& path);

}  // namespace cxpp::io
