// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#include "FoamWriter.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace ninja {

namespace {

namespace fs = std::filesystem;

void writeFoamHeader(std::ostream& out, const std::string& cls, const std::string& object) {
    out << "FoamFile\n"
           "{\n"
           "    version     2.0;\n"
           "    format      ascii;\n"
           "    class       "
        << cls
        << ";\n"
           "    object      "
        << object
        << ";\n"
           "}\n"
           "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
}

// %.15g round-trips integer-valued and simple decimal coordinates
// without printing spurious precision noise.
std::string fmtDouble(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.15g", v);
    return buf;
}

void writePointsFile(const GeneratedMesh& mesh, const fs::path& dir) {
    std::ofstream out(dir / "points");
    if (!out) {
        throw std::runtime_error("FoamWriter: cannot open points file for writing");
    }
    writeFoamHeader(out, "vectorField", "points");
    out << mesh.points.size() << "\n(\n";
    for (const Vec3& p : mesh.points) {
        out << "(" << fmtDouble(p.x) << " " << fmtDouble(p.y) << " " << fmtDouble(p.z) << ")\n";
    }
    out << ")\n";
}

void writeFacesFile(const GeneratedMesh& mesh, const fs::path& dir) {
    std::ofstream out(dir / "faces");
    if (!out) {
        throw std::runtime_error("FoamWriter: cannot open faces file for writing");
    }
    writeFoamHeader(out, "faceList", "faces");
    out << mesh.nFaces() << "\n(\n";
    for (int f = 0; f < mesh.nFaces(); ++f) {
        const IntSpan pts = mesh.faces.pointsOf(f);
        out << pts.size() << "(";
        for (int i = 0; i < pts.size(); ++i) {
            if (i != 0) {
                out << " ";
            }
            out << pts[i];
        }
        out << ")\n";
    }
    out << ")\n";
}

void writeOwnerFile(const GeneratedMesh& mesh, const fs::path& dir, int nCells) {
    std::ofstream out(dir / "owner");
    if (!out) {
        throw std::runtime_error("FoamWriter: cannot open owner file for writing");
    }
    writeFoamHeader(out, "labelList", "owner");
    out << "// nPoints:" << mesh.points.size() << " nCells:" << nCells
        << " nFaces:" << mesh.nFaces() << " nInternalFaces:" << mesh.nInternalFaces << "\n\n";
    out << mesh.nFaces() << "\n(\n";
    for (int f = 0; f < mesh.nFaces(); ++f) {
        out << mesh.faces.owner[static_cast<std::size_t>(f)] << "\n";
    }
    out << ")\n";
}

void writeNeighbourFile(const GeneratedMesh& mesh, const fs::path& dir, int nCells) {
    std::ofstream out(dir / "neighbour");
    if (!out) {
        throw std::runtime_error("FoamWriter: cannot open neighbour file for writing");
    }
    writeFoamHeader(out, "labelList", "neighbour");
    out << "// nPoints:" << mesh.points.size() << " nCells:" << nCells
        << " nFaces:" << mesh.nFaces() << " nInternalFaces:" << mesh.nInternalFaces << "\n\n";
    out << mesh.nInternalFaces << "\n(\n";
    for (int i = 0; i < mesh.nInternalFaces; ++i) {
        out << mesh.faces.neighbour[static_cast<std::size_t>(i)] << "\n";
    }
    out << ")\n";
}

void writeBoundaryFile(const GeneratedMesh& mesh, const fs::path& dir) {
    std::ofstream out(dir / "boundary");
    if (!out) {
        throw std::runtime_error("FoamWriter: cannot open boundary file for writing");
    }
    writeFoamHeader(out, "polyBoundaryMesh", "boundary");
    out << mesh.patches.size() << "\n(\n";
    for (const PatchInfo& p : mesh.patches) {
        out << "    " << p.name << "\n"
            << "    {\n"
            << "        type " << p.type << ";\n"
            << "        nFaces " << p.nFaces << ";\n"
            << "        startFace " << p.startFace << ";\n"
            << "    }\n";
    }
    out << ")\n";
}

void writeLabelListFile(const std::vector<int>& values, const fs::path& dir, const std::string& object) {
    std::ofstream out(dir / object);
    if (!out) {
        throw std::runtime_error("FoamWriter: cannot open " + object + " file for writing");
    }
    writeFoamHeader(out, "labelList", object);
    out << values.size() << "\n(\n";
    for (int v : values) {
        out << v << "\n";
    }
    out << ")\n";
}

// A volScalarField file under 0/: same FoamFile header as the polyMesh
// files plus the location entry, uniform-free internalField, and a
// wildcard zeroGradient boundaryField (accepted by both OpenFOAM and
// ParaView's reader for every patch type we emit).
void writeVolScalarField(const std::vector<double>& values, const fs::path& dir,
                         const std::string& object) {
    std::ofstream out(dir / object);
    if (!out) {
        throw std::runtime_error("FoamWriter: cannot open " + object + " file for writing");
    }
    out << "FoamFile\n"
           "{\n"
           "    version     2.0;\n"
           "    format      ascii;\n"
           "    class       volScalarField;\n"
           "    location    \"0\";\n"
           "    object      "
        << object
        << ";\n"
           "}\n"
           "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
    out << "dimensions      [0 0 0 0 0 0 0];\n\n";
    out << "internalField   nonuniform List<scalar>\n" << values.size() << "\n(\n";
    for (double v : values) {
        out << fmtDouble(v) << "\n";
    }
    out << ")\n;\n\n";
    out << "boundaryField\n{\n    \".*\"\n    {\n        type            zeroGradient;\n    }\n}\n";
}

// Per-cell volumes by the divergence theorem: each face loop is fanned
// into triangles about its centroid; each tetrahedron (origin, centroid,
// p0, p1) contributes +vol to the owner and -vol to the neighbour (face
// normals point owner->neighbour / outward on the boundary).
std::vector<double> computeCellVolumes(const GeneratedMesh& mesh) {
    std::vector<double> vol(static_cast<std::size_t>(mesh.nCells()), 0.0);
    for (int f = 0; f < mesh.nFaces(); ++f) {
        const IntSpan pts = mesh.faces.pointsOf(f);
        const int n = pts.size();
        Vec3 c{0.0, 0.0, 0.0};
        for (int i = 0; i < n; ++i) {
            const Vec3& p = mesh.points[static_cast<std::size_t>(pts[i])];
            c.x += p.x;
            c.y += p.y;
            c.z += p.z;
        }
        c.x /= n;
        c.y /= n;
        c.z /= n;
        double contrib = 0.0;
        for (int i = 0; i < n; ++i) {
            const Vec3& a = mesh.points[static_cast<std::size_t>(pts[i])];
            const Vec3& b = mesh.points[static_cast<std::size_t>(pts[(i + 1) % n])];
            // scalar triple product c . (a x b) / 6
            contrib += (c.x * (a.y * b.z - a.z * b.y) + c.y * (a.z * b.x - a.x * b.z) +
                        c.z * (a.x * b.y - a.y * b.x)) /
                       6.0;
        }
        vol[static_cast<std::size_t>(mesh.faces.owner[static_cast<std::size_t>(f)])] += contrib;
        if (f < mesh.nInternalFaces) {
            vol[static_cast<std::size_t>(mesh.faces.neighbour[static_cast<std::size_t>(f)])] -= contrib;
        }
    }
    return vol;
}

} // namespace

void writeDiagnosticFields(const GeneratedMesh& mesh, const std::vector<int>& cellLevel,
                           const std::vector<int>& cellLayerIndex, const std::string& caseDir) {
    const fs::path dir = fs::path(caseDir) / "0";
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec) {
        throw std::runtime_error("FoamWriter: cannot create directory '" + dir.string() +
                                  "': " + ec.message());
    }
    const std::size_t nCells = static_cast<std::size_t>(mesh.nCells());

    std::vector<double> refinementLevel(nCells, 0.0);
    for (std::size_t c = 0; c < nCells && c < cellLevel.size(); ++c) {
        refinementLevel[c] = cellLevel[c];
    }
    std::vector<double> isLayer(nCells, 0.0);
    std::vector<double> layerNumber(nCells, -1.0);
    for (std::size_t c = 0; c < nCells && c < cellLayerIndex.size(); ++c) {
        if (cellLayerIndex[c] >= 0) {
            isLayer[c] = 1.0;
            layerNumber[c] = cellLayerIndex[c];
        }
    }
    std::vector<double> cellSize = computeCellVolumes(mesh);
    for (double& v : cellSize) {
        v = std::cbrt(v > 0.0 ? v : 0.0);
    }

    writeVolScalarField(refinementLevel, dir, "refinementLevel");
    writeVolScalarField(isLayer, dir, "isLayer");
    writeVolScalarField(layerNumber, dir, "layerNumber");
    writeVolScalarField(cellSize, dir, "cellSize");
}

void writeLevelFields(const std::vector<int>& cellLevel, const std::vector<int>& pointLevel,
                       const std::string& caseDir) {
    const fs::path dir = fs::path(caseDir) / "constant" / "polyMesh";
    writeLabelListFile(cellLevel, dir, "cellLevel");
    writeLabelListFile(pointLevel, dir, "pointLevel");
}

void writeCellIntField(const std::vector<int>& values, const std::string& caseDir, const std::string& fieldName) {
    const fs::path dir = fs::path(caseDir) / "constant" / "polyMesh";
    writeLabelListFile(values, dir, fieldName);
}

void writePolyMesh(const GeneratedMesh& mesh, const std::string& caseDir) {
    const fs::path dir = fs::path(caseDir) / "constant" / "polyMesh";
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec) {
        throw std::runtime_error("FoamWriter: cannot create directory '" + dir.string() +
                                  "': " + ec.message());
    }

    // Empty marker file so the case opens directly in ParaView.
    std::ofstream foam(fs::path(caseDir) / "case.foam");
    if (!foam) {
        throw std::runtime_error("FoamWriter: cannot create case.foam marker file");
    }

    const int nCells = mesh.nCells();
    writePointsFile(mesh, dir);
    writeFacesFile(mesh, dir);
    writeOwnerFile(mesh, dir, nCells);
    writeNeighbourFile(mesh, dir, nCells);
    writeBoundaryFile(mesh, dir);
}

} // namespace ninja
