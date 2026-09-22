// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#include "VtkDebug.hpp"

#include <fstream>

namespace ninja {

void DebugVtk::write(const std::string& path) const {
    std::ofstream os(path);
    os << "# vtk DataFile Version 3.0\nninja layers debug\nASCII\nDATASET POLYDATA\n";
    os << "POINTS " << pts.size() << " double\n";
    for (const Vec3& p : pts) os << p.x << " " << p.y << " " << p.z << "\n";

    if (!polys.empty()) {
        std::size_t sz = 0;
        for (const auto& poly : polys) sz += poly.size() + 1;
        os << "POLYGONS " << polys.size() << " " << sz << "\n";
        for (const auto& poly : polys) {
            os << poly.size();
            for (int i : poly) os << " " << i;
            os << "\n";
        }
        if (!cellScalars.empty()) {
            os << "CELL_DATA " << polys.size() << "\n";
            for (const auto& [name, vals] : cellScalars) {
                os << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
                for (double v : vals) os << v << "\n";
            }
        }
    } else if (!lines.empty()) {
        os << "LINES " << lines.size() << " " << lines.size() * 3 << "\n";
        for (const auto& [a, b] : lines) os << "2 " << a << " " << b << "\n";
        if (!cellScalars.empty()) {
            os << "CELL_DATA " << lines.size() << "\n";
            for (const auto& [name, vals] : cellScalars) {
                os << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
                for (double v : vals) os << v << "\n";
            }
        }
    } else if (!vertexPts.empty()) {
        os << "VERTICES " << vertexPts.size() << " " << vertexPts.size() * 2 << "\n";
        for (int i : vertexPts) os << "1 " << i << "\n";
        if (!cellScalars.empty()) {
            os << "CELL_DATA " << vertexPts.size() << "\n";
            for (const auto& [name, vals] : cellScalars) {
                os << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
                for (double v : vals) os << v << "\n";
            }
        }
    }

    if (!pointScalars.empty() || !pointVectors.empty()) {
        os << "POINT_DATA " << pts.size() << "\n";
        for (const auto& [name, vals] : pointScalars) {
            os << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
            for (double v : vals) os << v << "\n";
        }
        for (const auto& [name, vals] : pointVectors) {
            os << "VECTORS " << name << " double\n";
            for (const Vec3& v : vals) os << v.x << " " << v.y << " " << v.z << "\n";
        }
    }
}

} // namespace ninja
