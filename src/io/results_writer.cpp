#include "calculixpp/io/results_writer.hpp"

#include <cstdio>
#include <fstream>
#include <stdexcept>

namespace cxpp::io {
namespace {

// CGX/.frd element type codes (cgx element type numbering).
int frd_elem_type(ElementType t) {
  switch (t) {
    case ElementType::C3D8:
    case ElementType::C3D8R:
      return 1;
    case ElementType::C3D6:
      return 2;
    case ElementType::C3D4:
      return 3;
    case ElementType::C3D20:
    case ElementType::C3D20R:
      return 4;
    case ElementType::C3D15:
      return 5;
    case ElementType::C3D10:
      return 6;
  }
  return 3;
}

std::string e12_5(Real v) {
  char buf[24];
  std::snprintf(buf, sizeof(buf), "%12.5E", v);
  return buf;
}
std::string i10(long v) {
  char buf[24];
  std::snprintf(buf, sizeof(buf), "%10ld", v);
  return buf;
}

// .frd order is XX,YY,ZZ,XY,YZ,ZX; ours is xx,yy,zz,xy,xz,yz (swap last two).
Voigt6 to_frd_order(const Voigt6& s) {
  return {s[0], s[1], s[2], s[3], s[5], s[4]};
}

// Write the .frd header + node + element mesh blocks (shared by static and thermal
// result writers). Leaves the stream positioned to append result datasets.
void write_frd_mesh(std::ofstream& out, const Mesh& mesh) {
  out << "    1C\n";
  out << "    1UUSER   CalculiX++\n";

  // Nodes.
  out << "    2C" << i10(static_cast<long>(mesh.num_nodes())) << "                  2\n";
  for (const Node& nd : mesh.nodes())
    out << " -1" << i10(nd.id) << e12_5(nd.x[0]) << e12_5(nd.x[1]) << e12_5(nd.x[2]) << "\n";
  out << " -3\n";

  // Elements.
  out << "    3C" << i10(static_cast<long>(mesh.num_elements())) << "                  2\n";
  for (const Element& el : mesh.elements()) {
    out << " -1" << i10(el.id) << i10(frd_elem_type(el.type)).substr(5) << "    0    1\n";
    out << " -2";
    for (std::size_t k = 0; k < el.nodes.size(); ++k) {
      out << i10(el.nodes[k]);
      if ((k + 1) % 10 == 0 && k + 1 < el.nodes.size()) out << "\n -2";
    }
    out << "\n";
  }
  out << " -3\n";
}

}  // namespace

void write_frd(const std::string& path, const Model& model, const StaticFields& f) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot open .frd for writing: " + path);
  const Mesh& mesh = model.mesh;
  write_frd_mesh(out, mesh);

  auto dataset_header = [&](const char* name, int ncomp) {
    out << "  100CL  101" << e12_5(1.0) << i10(static_cast<long>(mesh.num_nodes()))
        << "                     0    1           1\n";
    out << " -4  " << name;
    for (int i = static_cast<int>(std::string(name).size()); i < 8; ++i) out << ' ';
    out << i10(ncomp).substr(5) << "    1\n";
  };

  // DISP (U).
  dataset_header("DISP", 4);
  for (const char* c : {"D1", "D2", "D3"})
    out << " -5  " << c << "          1    2    " << (c[1]) << "    0\n";
  for (std::size_t i = 0; i < mesh.num_nodes(); ++i)
    out << " -1" << i10(mesh.nodes()[i].id) << e12_5(f.displacement[i][0])
        << e12_5(f.displacement[i][1]) << e12_5(f.displacement[i][2]) << "\n";
  out << " -3\n";

  // STRESS (S).
  dataset_header("STRESS", 6);
  for (const char* c : {"SXX", "SYY", "SZZ", "SXY", "SYZ", "SZX"})
    out << " -5  " << c << "         1    4    1    1\n";
  for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
    const Voigt6 s = to_frd_order(f.stress[i]);
    out << " -1" << i10(mesh.nodes()[i].id);
    for (int c = 0; c < 6; ++c) out << e12_5(s[static_cast<std::size_t>(c)]);
    out << "\n";
  }
  out << " -3\n";

  // TOSTRAIN (E).
  dataset_header("TOSTRAIN", 6);
  for (const char* c : {"EXX", "EYY", "EZZ", "EXY", "EYZ", "EZX"})
    out << " -5  " << c << "         1    4    1    1\n";
  for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
    const Voigt6 e = to_frd_order(f.strain[i]);
    out << " -1" << i10(mesh.nodes()[i].id);
    for (int c = 0; c < 6; ++c) out << e12_5(e[static_cast<std::size_t>(c)]);
    out << "\n";
  }
  out << " -3\n";

  // FORC (RF).
  dataset_header("FORC", 4);
  for (const char* c : {"F1", "F2", "F3"})
    out << " -5  " << c << "          1    2    " << (c[1]) << "    0\n";
  for (std::size_t i = 0; i < mesh.num_nodes(); ++i)
    out << " -1" << i10(mesh.nodes()[i].id) << e12_5(f.reaction[i][0])
        << e12_5(f.reaction[i][1]) << e12_5(f.reaction[i][2]) << "\n";
  out << " -3\n";
  out << " 9999\n";
}

void write_dat(const std::string& path, const Model& model, const StaticFields& f) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot open .dat for writing: " + path);
  const Mesh& mesh = model.mesh;

  out << "\n displacements (vx,vy,vz)\n\n";
  for (std::size_t i = 0; i < mesh.num_nodes(); ++i)
    out << i10(mesh.nodes()[i].id) << e12_5(f.displacement[i][0])
        << e12_5(f.displacement[i][1]) << e12_5(f.displacement[i][2]) << "\n";

  out << "\n stresses (sxx,syy,szz,sxy,sxz,syz)\n\n";
  for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
    out << i10(mesh.nodes()[i].id);
    for (int c = 0; c < 6; ++c) out << e12_5(f.stress[i][static_cast<std::size_t>(c)]);
    out << "\n";
  }

  out << "\n strains (exx,eyy,ezz,exy,exz,eyz)\n\n";
  for (std::size_t i = 0; i < mesh.num_nodes(); ++i) {
    out << i10(mesh.nodes()[i].id);
    for (int c = 0; c < 6; ++c) out << e12_5(f.strain[i][static_cast<std::size_t>(c)]);
    out << "\n";
  }

  out << "\n forces (fx,fy,fz)\n\n";
  for (std::size_t i = 0; i < mesh.num_nodes(); ++i)
    out << i10(mesh.nodes()[i].id) << e12_5(f.reaction[i][0])
        << e12_5(f.reaction[i][1]) << e12_5(f.reaction[i][2]) << "\n";
}

void write_frd(const std::string& path, const Model& model, const ThermalFields& t) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot open .frd for writing: " + path);
  const Mesh& mesh = model.mesh;
  write_frd_mesh(out, mesh);

  // NDTEMP scalar temperature dataset (CGX field name NDTEMP / component T1).
  out << "  100CL  101" << e12_5(1.0) << i10(static_cast<long>(mesh.num_nodes()))
      << "                     0    1           1\n";
  out << " -4  NDTEMP      1    1\n";
  out << " -5  T           1    1    0    0\n";
  for (std::size_t i = 0; i < mesh.num_nodes(); ++i)
    out << " -1" << i10(mesh.nodes()[i].id) << e12_5(t.temperature[i]) << "\n";
  out << " -3\n";

  // FLUX vector dataset (nodal, extrapolated heat flux HFL; CGX field FLUX / F1..F3).
  if (!t.heat_flux.empty()) {
    out << "  100CL  101" << e12_5(1.0) << i10(static_cast<long>(mesh.num_nodes()))
        << "                     0    1           1\n";
    out << " -4  FLUX        4    1\n";
    for (const char* c : {"F1", "F2", "F3"})
      out << " -5  " << c << "          1    2    " << (c[1]) << "    0\n";
    for (std::size_t i = 0; i < mesh.num_nodes(); ++i)
      out << " -1" << i10(mesh.nodes()[i].id) << e12_5(t.heat_flux[i][0])
          << e12_5(t.heat_flux[i][1]) << e12_5(t.heat_flux[i][2]) << "\n";
    out << " -3\n";
  }
  out << " 9999\n";
}

void write_dat(const std::string& path, const Model& model, const ThermalFields& t) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot open .dat for writing: " + path);
  const Mesh& mesh = model.mesh;

  out << "\n temperatures (t)\n\n";
  for (std::size_t i = 0; i < mesh.num_nodes(); ++i)
    out << i10(mesh.nodes()[i].id) << e12_5(t.temperature[i]) << "\n";

  out << "\n heat generation (rfl)\n\n";
  for (std::size_t i = 0; i < mesh.num_nodes(); ++i)
    out << i10(mesh.nodes()[i].id) << e12_5(t.flux_reaction[i]) << "\n";

  // Integration-point heat flux HFL (*EL PRINT HFL), matching CalculiX's block layout
  // "heat flux (elem, integ.pnt.,qx,qy,qz)": one row per element per Gauss point.
  if (!t.hfl_points.empty()) {
    out << "\n heat flux (elem, integ.pnt.,qx,qy,qz)\n\n";
    for (const HeatFluxPoint& p : t.hfl_points)
      out << i10(p.elem_id) << i10(p.gp).substr(6) << e12_5(p.flux[0])
          << e12_5(p.flux[1]) << e12_5(p.flux[2]) << "\n";
  }
}

}  // namespace cxpp::io
