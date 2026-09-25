// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#include <mpi.h>
#include <sys/resource.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <sstream>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "BaseMesh.hpp"
#include "CutData.hpp"
#include "Cutter.hpp"
#include "Dict.hpp"
#include "FoamWriter.hpp"
#include "Geometry.hpp"
#include "Layers.hpp"
#include "Par.hpp"
#include "Refine.hpp"
#include "Stl.hpp"
#include "VtkDebug.hpp"

namespace {

constexpr const char* kVersion = "0.4.0";

// Opt-in wall-clock stage timing, stderr-only, enabled by
// NINJA_STAGE_TIMES=1 (never printed in normal runs, so no gate/log
// output is affected).
bool stageTimesOn() {
    static const bool on = std::getenv("NINJA_STAGE_TIMES") != nullptr;
    return on;
}

// Opt-in per-stage memory attribution, stderr-only, enabled by
// NINJA_MEM_LOG=1 -- prints this process's current VmRSS and peak VmHWM
// (from /proc/self/status) at every stage-timer boundary. Zero output
// (and one cached getenv) when unset; never touches stdout, so no
// gate/log output is affected.
bool memLogOn() {
    static const bool on = std::getenv("NINJA_MEM_LOG") != nullptr;
    return on;
}
void memLog(const char* name) {
    if (!memLogOn()) {
        return;
    }
    long rssKb = -1, hwmKb = -1;
    std::ifstream st("/proc/self/status");
    std::string line;
    while (std::getline(st, line)) {
        if (line.rfind("VmRSS:", 0) == 0) {
            rssKb = std::atol(line.c_str() + 6);
        } else if (line.rfind("VmHWM:", 0) == 0) {
            hwmKb = std::atol(line.c_str() + 6);
        }
    }
    std::cerr << "mem " << name << " rank " << ninja::parRank() << ": VmRSS " << rssKb
              << " kB, VmHWM " << hwmKb << " kB\n";
}

struct StageTimer {
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    void lap(const char* name) {
        if (stageTimesOn()) {
            const auto t1 = std::chrono::steady_clock::now();
            std::cerr << "stage " << name << " rank " << ninja::parRank() << ": "
                      << std::chrono::duration<double>(t1 - t0).count() << " s\n";
        }
        memLog(name);
        t0 = std::chrono::steady_clock::now();
    }
};

// Rank-local run context, owned by main() after MPI_Init — no globals,
// no singletons.
struct RunContext {
    int rank = 0;
    int size = 1;
};

std::vector<std::string> requiredSideKeys() {
    return {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
}

ninja::MeshConfig buildMeshConfig(const ninja::DictEntry& root) {
    ninja::MeshConfig cfg;

    const ninja::DictEntry& domain = ninja::requireEntry(root, "domain");
    auto parseVec3 = [&](const std::string& key) {
        std::vector<std::string> comps = ninja::requireList(domain, key);
        if (comps.size() != 3) {
            throw std::runtime_error("Dict error: key 'domain." + key + "' must have 3 components");
        }
        ninja::Vec3 v;
        v.x = std::stod(comps[0]);
        v.y = std::stod(comps[1]);
        v.z = std::stod(comps[2]);
        return v;
    };
    cfg.min = parseVec3("min");
    cfg.max = parseVec3("max");

    std::vector<std::string> nComps = ninja::requireList(domain, "n");
    if (nComps.size() != 3) {
        throw std::runtime_error("Dict error: key 'domain.n' must have 3 components");
    }
    cfg.nx = std::stoi(nComps[0]);
    cfg.ny = std::stoi(nComps[1]);
    cfg.nz = std::stoi(nComps[2]);

    const ninja::DictEntry& patches = ninja::requireEntry(root, "patches");
    for (const std::string& sideKey : requiredSideKeys()) {
        const ninja::DictEntry& side = ninja::requireEntry(patches, sideKey);
        ninja::PatchSpec spec;
        spec.sideKey = sideKey;
        spec.type = ninja::requireScalar(side, "type");
        spec.name = ninja::requireScalar(side, "name");
        cfg.patches.push_back(spec);
    }

    // Warn on unknown top-level keys.
    for (const auto& [key, value] : root.dict) {
        if (key != "domain" && key != "patches" && key != "FoamFile" && key != "geometry" &&
            key != "locationInMesh" && key != "refinement" && key != "layers" && key != "layersDebug" &&
            key != "refinementGeometry") {
            std::cerr << "warning: unknown top-level key '" << key << "' ignored\n";
        }
    }

    return cfg;
}

ninja::Vec3 parseVec3Entry(const ninja::DictEntry& d, const std::string& key) {
    std::vector<std::string> comps = ninja::requireList(d, key);
    if (comps.size() != 3) {
        throw std::runtime_error("Dict error: key '" + key + "' must have 3 components");
    }
    ninja::Vec3 v;
    v.x = std::stod(comps[0]);
    v.y = std::stod(comps[1]);
    v.z = std::stod(comps[2]);
    return v;
}

// refinement { zoneN { type box|sphere|surface; ...; level N; } }
// Absent `refinement` => no refinement regions.
// `type surface` takes `stl <rawFile>; distance <d>;` -- requires
// `geometry` to be present and `rawFile` declared there (clear fatal
// error otherwise, checked here since both dicts are available at this
// call site).
std::vector<ninja::RefineRegion> buildRefineRegions(const ninja::DictEntry& root, const ninja::GeometryConfig& geom) {
    std::vector<ninja::RefineRegion> regions;
    auto it = root.dict.find("refinement");
    if (it == root.dict.end()) {
        return regions;
    }
    const ninja::DictEntry& refDict = it->second;
    if (refDict.kind != ninja::DictEntry::Kind::Dict) {
        throw std::runtime_error("Dict error: 'refinement' must be a sub-dict");
    }
    for (const auto& [name, entry] : refDict.dict) {
        if (entry.kind != ninja::DictEntry::Kind::Dict) {
            throw std::runtime_error("Dict error: refinement region '" + name + "' must be a sub-dict");
        }
        ninja::RefineRegion r;
        r.type = ninja::requireScalar(entry, "type");
        r.level = std::stoi(ninja::requireScalar(entry, "level"));
        if (r.type == "box") {
            r.min = parseVec3Entry(entry, "min");
            r.max = parseVec3Entry(entry, "max");
        } else if (r.type == "sphere") {
            r.centre = parseVec3Entry(entry, "centre");
            r.radius = std::stod(ninja::requireScalar(entry, "radius"));
        } else if (r.type == "surface") {
            r.stlKey = ninja::requireScalar(entry, "stl");
            r.distance = std::stod(ninja::requireScalar(entry, "distance"));
            bool found = false;
            for (const std::vector<ninja::StlEntry>* list : {&geom.stls, &geom.sizing}) {
                for (const ninja::StlEntry& se : *list) {
                    if (se.rawFile == r.stlKey) found = true;
                }
            }
            if (!found) {
                throw std::runtime_error("Dict error: refinement region '" + name + "' references stl '" +
                                          r.stlKey + "' which is declared in neither 'geometry' nor"
                                          " 'refinementGeometry'");
            }
        } else {
            throw std::runtime_error("Dict error: refinement region '" + name + "' has unknown type '" +
                                      r.type + "'");
        }
        regions.push_back(r);
    }
    return regions;
}

// geometry { <file1.stl> { name ...; } <file2.stl> { name ...; } ... }
// and locationInMesh (x y z). Absent/empty `geometry` => no cut.
// Any number of STL entries is accepted, each its own named wall
// patch, in dict declaration order (per
// `DictEntry::order` — `geomDict.dict` itself is alphabetical and MUST
// NOT be iterated for patch ordering).
ninja::GeometryConfig buildGeometryConfig(const ninja::DictEntry& root, const std::string& caseDir) {
    ninja::GeometryConfig gc;
    // refinementGeometry { <file.stl> { } ... }: surfaces that only size
    // the grid (`type surface` refinement rules). A file may not be both
    // a wall and a sizing surface -- one role per STL keeps "which
    // surfaces bound the fluid" answerable from the geometry block alone.
    if (auto sz = root.dict.find("refinementGeometry"); sz != root.dict.end()) {
        if (sz->second.kind != ninja::DictEntry::Kind::Dict) {
            throw std::runtime_error("Dict error: 'refinementGeometry' must be a sub-dict");
        }
        for (const std::string& stlFile : sz->second.order) {
            if (sz->second.dict.at(stlFile).kind != ninja::DictEntry::Kind::Dict) {
                throw std::runtime_error("Dict error: refinementGeometry entry '" + stlFile +
                                          "' must be a sub-dict (e.g. '" + stlFile + " { }')");
            }
            ninja::StlEntry se;
            se.stlPath = (std::filesystem::path(caseDir) / stlFile).string();
            se.rawFile = stlFile;
            gc.sizing.push_back(std::move(se));
        }
    }
    auto it = root.dict.find("geometry");
    if (it == root.dict.end()) {
        return gc;
    }
    const ninja::DictEntry& geomDict = it->second;
    if (geomDict.kind != ninja::DictEntry::Kind::Dict) {
        throw std::runtime_error("Dict error: 'geometry' must be a sub-dict");
    }
    if (geomDict.dict.empty()) {
        return gc;
    }

    for (const std::string& stlFile : geomDict.order) {
        const ninja::DictEntry& stlEntry = geomDict.dict.at(stlFile);
        if (stlEntry.kind != ninja::DictEntry::Kind::Dict) {
            throw std::runtime_error("Dict error: geometry entry '" + stlFile + "' must be a sub-dict");
        }
        for (const ninja::StlEntry& sz : gc.sizing) {
            if (sz.rawFile == stlFile) {
                throw std::runtime_error("Dict error: '" + stlFile +
                                          "' is declared in both 'geometry' and 'refinementGeometry'");
            }
        }
        ninja::StlEntry se;
        se.stlPath = (std::filesystem::path(caseDir) / stlFile).string();
        se.patchName = ninja::requireScalar(stlEntry, "name");
        se.rawFile = stlFile;
        gc.stls.push_back(std::move(se));
    }
    gc.present = true;

    std::vector<std::string> loc = ninja::requireList(root, "locationInMesh");
    if (loc.size() != 3) {
        throw std::runtime_error("Dict error: 'locationInMesh' must have 3 components");
    }
    gc.locationInMesh.x = std::stod(loc[0]);
    gc.locationInMesh.y = std::stod(loc[1]);
    gc.locationInMesh.z = std::stod(loc[2]);
    return gc;
}

// --- `layers` dict block -------------------------------------------------
//
// layers
// {
//     sphere.stl { }                                   // every default
//     wall.stl   { nLayers 6; expansionRatio 1.3; finalLayerThickness 0.4; }
//     absorbVolFrac 0.3;   // optional, default 0.3
// }
//
// LAYER THICKNESS IS ALWAYS RELATIVE TO THE LOCAL CELL, never an
// absolute length. `finalLayerThickness` is the OUTERMOST (thickest)
// layer as a fraction of the undistorted cell it grows into -- snappy's
// `relativeSizes true` convention, and the form every case here was
// originally written in before being hand-converted to absolute lengths.
// The total is derived:
//     total = finalLayerThickness * h * sum_{i<nLayers} expansionRatio^-i
// with h = dx_base / 2^level and `level` from refineLevelForStl below.
//
// WHY ABSOLUTE LENGTHS WERE REMOVED. Layer quality depends on t/h, and h
// varies with refinement, so an absolute t makes the dict author
// hand-compute the conversion once per surface -- bm_layers_gate's header
// literally carried that table -- and silently misstates the request
// wherever the grid is not the size the author assumed. Making the
// relative form the ONLY form deletes the conversion step, and with it
// the class of mistake where a passage is refined specifically to be
// resolved and a thickness sized elsewhere closes it anyway.
//
// Defaults applied per key when absent: nLayers 4, expansionRatio 1.5,
// finalLayerThickness 0.3 -- so a bare `some.stl { }` is a valid and
// sensible layer request. 0.3 is the more conservative of the two values
// the shipped cases were built at (0.3 and 0.4): both are far from the
// aspect-ratio floor, so the discriminator is that the offset cut closes
// passages narrower than 2t, which 0.3 keeps at ~1.4 cells against
// ~1.9 for 0.4.
//
// `totalThickness` and `firstLayerThickness` are hard errors whose
// message quotes the equivalent fraction for that entry.
struct LayersConfig {
    bool present = false;
    // total thickness, keyed by raw `geometry` dict filename (matches
    // StlEntry::rawFile / Triangle::solidId ordinal via geometry order).
    std::unordered_map<std::string, double> thicknessByRawFile;
    // Layer-by-layer march: the march
    // advances one graded layer step at a time. Defaults (absent keys):
    // 1 layer, ratio 1.0. For combo 3 the stored ratio is the DERIVED
    // one.
    std::unordered_map<std::string, int> nLayersByRawFile;
    std::unordered_map<std::string, double> ratioByRawFile;
    // Per-STL smoothing radius fraction, from the dict's
    // `smoothRadius` key. Absent keys keep LayerStlSpec's production
    // default (1.0); 0 disables smoothing for that STL.
    std::unordered_map<std::string, double> smoothRadiusByRawFile;
    // `localThickness true` (default false): the total thickness PER UNIT
    // CELL SIZE (finalLayerThickness * sum ratio^-i); the offset cut and the
    // march take the thickness from ThicknessField instead of the constant
    // above.
    std::unordered_map<std::string, double> localPerHByRawFile;
    double absorbVolFrac = 0.3;
    // Test-only forced-drop scaffold, parsed from a top-level
    // `layersDebug { forceDropSphere (x y z r); }` block. Ignored when
    // absent; never used outside tests.
    std::vector<ninja::ForceDropSphere> debugForceDropSpheres;
};

// Layer-thickness defaults, applied per `layers` entry when the key is
// absent. See the block comment above for why 0.3.
constexpr double kDefaultFinalLayerFrac = 0.3;
constexpr int kDefaultNLayers = 4;
constexpr double kDefaultExpansionRatio = 1.5;

// Total thickness in units of the FINAL (outermost, thickest) layer,
// sum_{i=0}^{n-1} r^-i. With t_final = t0*r^(n-1) and the total the
// geometric sum t0*(r^n - 1)/(r - 1), total / t_final is exactly this.
double finalToTotalFactor(double expansionRatio, int nLayers) {
    if (nLayers < 1) {
        throw std::runtime_error("Dict error: 'layers' entry 'nLayers' must be >= 1");
    }
    if (std::fabs(expansionRatio - 1.0) < 1e-12) {
        return static_cast<double>(nLayers);
    }
    double s = 0.0;
    for (int i = 0; i < nLayers; ++i) {
        s += std::pow(expansionRatio, -i);
    }
    return s;
}

// The refinement level that applies to an STL: the deepest `type
// surface` region naming it, else 0 (base).
//
// Box and sphere regions are deliberately NOT considered. They are not
// tied to a surface, and a box that merely overlaps a wall says nothing
// about the cell size AT that wall -- only a surface rule does.
//
// One level per STL: the level the whole STL's thickness is sized for,
// exact only while the STL sits at one level (bm_layers_wfp's
// `wfp_fw.stl` spans levels 2, 3 and 4 and takes level 2 here). With
// `localThickness true` it is only the REFERENCE the thickness is written
// against -- every face rescales it to its own cell.
int refineLevelForStl(const std::vector<ninja::RefineRegion>& regions, const std::string& rawFile) {
    int lvl = 0;
    for (const ninja::RefineRegion& r : regions) {
        if (r.type == "surface" && r.stlKey == rawFile) {
            lvl = std::max(lvl, r.level);
        }
    }
    return lvl;
}

LayersConfig buildLayersConfig(const ninja::DictEntry& root, const ninja::GeometryConfig& geom,
                                const ninja::MeshConfig& cfg,
                                const std::vector<ninja::RefineRegion>& regions) {
    LayersConfig lc;
    auto it = root.dict.find("layers");
    if (it == root.dict.end()) {
        return lc;
    }
    const ninja::DictEntry& ld = it->second;
    if (ld.kind != ninja::DictEntry::Kind::Dict) {
        throw std::runtime_error("Dict error: 'layers' must be a sub-dict");
    }
    lc.present = true;

    if (ld.dict.find("absorbVolFrac") != ld.dict.end()) {
        lc.absorbVolFrac = std::stod(ninja::requireScalar(ld, "absorbVolFrac"));
        // A volume FRACTION, so it means nothing outside [0, 1]. 0 is
        // legal and disables absorption: no volFrac is ever below it.
        //
        // A NEGATIVE value is the case this really guards. cutMesh reads
        // `volFracThreshold >= 0.0 ? volFracThreshold : kSmallVolFrac`,
        // so it treats any negative number as "not supplied" and
        // silently falls back to the PLAIN-cut threshold of 1e-3 -- 300x
        // smaller than the 0.3 default, changing which cells are
        // absorbed rather than dropped, with nothing printed. Refusing
        // here is the difference between a dict error and a mesh that is
        // quietly not the one that was asked for.
        if (!(lc.absorbVolFrac >= 0.0 && lc.absorbVolFrac <= 1.0)) {
            std::ostringstream m;
            m << "Dict error: 'layers.absorbVolFrac' must be in [0, 1] -- it is the fraction of its"
                 " UNCUT volume a cut cell must keep to avoid being absorbed into a neighbour. Got "
              << lc.absorbVolFrac << ".";
            throw std::runtime_error(m.str());
        }
    }

    for (const std::string& key : ld.order) {
        if (key == "absorbVolFrac") {
            continue;
        }
        const ninja::DictEntry& entry = ld.dict.at(key);
        if (entry.kind != ninja::DictEntry::Kind::Dict) {
            throw std::runtime_error("Dict error: 'layers' entry '" + key + "' must be a sub-dict");
        }
        bool declared = false;
        for (const ninja::StlEntry& se : geom.stls) {
            if (se.rawFile == key) {
                declared = true;
                break;
            }
        }
        if (!declared) {
            for (const ninja::StlEntry& se : geom.sizing) {
                if (se.rawFile == key) {
                    throw std::runtime_error("Dict error: 'layers' entry '" + key +
                                              "' is a 'refinementGeometry' surface; layers grow only on"
                                              " 'geometry' walls");
                }
            }
            throw std::runtime_error("Dict error: 'layers' entry '" + key +
                                      "' references an stl not declared in 'geometry'");
        }

        // Thickness is RELATIVE: the outermost layer as a fraction of
        // the undistorted local cell. h uses the x-direction base size
        // and this STL's refinement level, matching the march's own
        // h_local (Layers.cpp: dx0 / 2^cellLevel), so the guard
        // thresholds there are measured against the same length.
        const int stlLevel = refineLevelForStl(regions, key);
        const double dx0 = (cfg.max.x - cfg.min.x) / static_cast<double>(std::max(1, cfg.nx));
        const double hLocal = std::ldexp(dx0, -stlLevel);

        int nLayers = kDefaultNLayers;
        if (entry.dict.find("nLayers") != entry.dict.end()) {
            nLayers = std::stoi(ninja::requireScalar(entry, "nLayers"));
            if (nLayers < 1) {
                throw std::runtime_error("Dict error: 'layers." + key + ".nLayers' must be >= 1");
            }
        }
        double expansionRatio = kDefaultExpansionRatio;
        if (entry.dict.find("expansionRatio") != entry.dict.end()) {
            expansionRatio = std::stod(ninja::requireScalar(entry, "expansionRatio"));
            // A non-positive ratio is a dict error, caught here rather
            // than silently coerced deep in the march (Layers.cpp).
            if (!(expansionRatio > 0.0)) {
                throw std::runtime_error("Dict error: 'layers." + key + ".expansionRatio' must be > 0");
            }
        }
        const double factor = finalToTotalFactor(expansionRatio, nLayers);

        // The removed ABSOLUTE keys. Refusing them is not enough: the
        // message quotes the fraction that reproduces the very length
        // the dict asked for, so migrating an old dict is mechanical
        // and cannot be got wrong by re-deriving it by hand -- which is
        // exactly the error this whole change exists to remove.
        auto refuseAbsolute = [&](const char* oldKey, double equivFinal) {
            std::ostringstream m;
            m << "Dict error: 'layers." << key << "." << oldKey
              << "' is no longer accepted -- layer thickness is now RELATIVE to the local cell. Write"
              << " 'finalLayerThickness " << (equivFinal / hLocal) << ";' instead (local cell h = " << hLocal
              << " m = base dx " << dx0 << " / 2^" << stlLevel << ", this STL's 'refinement' level).";
            throw std::runtime_error(m.str());
        };
        if (entry.dict.find("totalThickness") != entry.dict.end()) {
            refuseAbsolute("totalThickness",
                           std::stod(ninja::requireScalar(entry, "totalThickness")) / factor);
        }
        if (entry.dict.find("firstLayerThickness") != entry.dict.end()) {
            refuseAbsolute("firstLayerThickness",
                           std::stod(ninja::requireScalar(entry, "firstLayerThickness")) *
                               std::pow(expansionRatio, nLayers - 1));
        }
        // Accepted only as the identity it now always carries, so a
        // snappy dict copied verbatim either reads correctly or says why.
        if (entry.dict.find("relativeSizes") != entry.dict.end()) {
            if (ninja::requireScalar(entry, "relativeSizes") != "true") {
                throw std::runtime_error("Dict error: 'layers." + key +
                                          ".relativeSizes' must be 'true' -- absolute layer sizes are not"
                                          " supported; thickness is always a fraction of the local cell");
            }
        }

        double finalFrac = kDefaultFinalLayerFrac;
        if (entry.dict.find("finalLayerThickness") != entry.dict.end()) {
            finalFrac = std::stod(ninja::requireScalar(entry, "finalLayerThickness"));
            if (!(finalFrac > 0.0)) {
                throw std::runtime_error("Dict error: 'layers." + key +
                                          ".finalLayerThickness' must be > 0 (a fraction of the local cell)");
            }
        }

        lc.thicknessByRawFile[key] = finalFrac * hLocal * factor;
        // `localThickness true`: the thickness follows the LOCAL cell, face by
        // face (snappy's relativeSizes), instead of the one level above --
        // for a wall whose parts are refined to different levels.
        bool local = false;
        if (entry.dict.find("localThickness") != entry.dict.end()) {
            const std::string v = ninja::requireScalar(entry, "localThickness");
            if (v != "true" && v != "false") {
                throw std::runtime_error("Dict error: 'layers." + key + ".localThickness' must be true or false");
            }
            local = v == "true";
        }
        if (local) lc.localPerHByRawFile[key] = finalFrac * factor;
        lc.nLayersByRawFile[key] = nLayers;
        lc.ratioByRawFile[key] = expansionRatio;

        // `smoothRadius`, e.g. `layers { fp.stl { ...; smoothRadius 1.0; } }`.
        // Optional; absent keeps the production default (1.0), 0 disables.
        if (entry.dict.find("smoothRadius") != entry.dict.end()) {
            const double mf = std::stod(ninja::requireScalar(entry, "smoothRadius"));
            if (mf < 0.0) {
                throw std::runtime_error("Dict error: 'layers." + key + ".smoothRadius' must be >= 0");
            }
            lc.smoothRadiusByRawFile[key] = mf;
        }
    }
    return lc;
}

// layersDebug { forceDropSphere (x y z r); } -- forced-drop test
// scaffold. Absent => empty (no forced drops); test-only, never used
// outside tests.
void buildLayersDebugConfig(const ninja::DictEntry& root, LayersConfig& lc) {
    auto it = root.dict.find("layersDebug");
    if (it == root.dict.end()) {
        return;
    }
    const ninja::DictEntry& dbg = it->second;
    if (dbg.kind != ninja::DictEntry::Kind::Dict) {
        throw std::runtime_error("Dict error: 'layersDebug' must be a sub-dict");
    }
    auto sIt = dbg.dict.find("forceDropSphere");
    if (sIt != dbg.dict.end()) {
        std::vector<std::string> comps = ninja::requireList(dbg, "forceDropSphere");
        if (comps.size() != 4) {
            throw std::runtime_error("Dict error: 'layersDebug.forceDropSphere' must have 4 components (x y z r)");
        }
        ninja::ForceDropSphere sph;
        sph.centre.x = std::stod(comps[0]);
        sph.centre.y = std::stod(comps[1]);
        sph.centre.z = std::stod(comps[2]);
        sph.radius = std::stod(comps[3]);
        lc.debugForceDropSpheres.push_back(sph);
    }
}

// Reads every STL in `gc.stls`, stamping each triangle's `solidId` with
// its declaration ordinal. Classification/intercepts run on the
// combined triangle soup — the union of solids; STLs are assumed
// disjoint and non-nested.
std::vector<ninja::Triangle> readCombinedStls(const ninja::GeometryConfig& gc) {
    std::vector<ninja::Triangle> combined;
    for (std::size_t s = 0; s < gc.stls.size(); ++s) {
        std::vector<ninja::Triangle> tris = ninja::readStl(gc.stls[s].stlPath);
        for (ninja::Triangle& t : tris) {
            t.solidId = static_cast<int>(s);
        }
        combined.insert(combined.end(), tris.begin(), tris.end());
    }
    return combined;
}

// Keyed by raw geometry-dict filename (`type surface` regions
// reference an STL that way, not by patch name). Each STL is
// read individually (not combined) since distance-mode regions query a
// single named surface, not the union soup.
std::unordered_map<std::string, std::vector<ninja::Triangle>> readStlsByRawFile(const ninja::GeometryConfig& gc) {
    std::unordered_map<std::string, std::vector<ninja::Triangle>> out;
    for (const ninja::StlEntry& se : gc.stls) {
        out[se.rawFile] = ninja::readStl(se.stlPath);
    }
    for (const ninja::StlEntry& se : gc.sizing) {
        out[se.rawFile] = ninja::readStl(se.stlPath);
    }
    return out;
}

// One triangle list per STL, in `geom.stls` declaration order (== the
// `Triangle::solidId` ordinal readCombinedStls stamps). The offset
// classification/root-finder needs a PER-STL TriangleAabbBins
// (`anyTriangleWithin(bins_s, p, t_s)`), separate from the
// combined-soup bins parity still uses.
std::vector<std::vector<ninja::Triangle>> readPerStlTriangles(const ninja::GeometryConfig& gc) {
    std::vector<std::vector<ninja::Triangle>> out;
    out.reserve(gc.stls.size());
    for (const ninja::StlEntry& se : gc.stls) {
        out.push_back(ninja::readStl(se.stlPath));
    }
    return out;
}

// Per-STL total thickness vector aligned with `geom.stls` order (0.0 for
// an STL with no `layers` entry -- t = 0 => that STL cuts at the true
// wall).
std::vector<double> stlThicknessOf(const ninja::GeometryConfig& gc, const LayersConfig& layers) {
    std::vector<double> t(gc.stls.size(), 0.0);
    for (std::size_t s = 0; s < gc.stls.size(); ++s) {
        auto it = layers.thicknessByRawFile.find(gc.stls[s].rawFile);
        if (it != layers.thicknessByRawFile.end()) {
            t[s] = it->second;
        }
    }
    return t;
}

// Per-STL thickness PER UNIT CELL SIZE for `localThickness` STLs, 0 for
// the rest (aligned with `geom.stls`).
std::vector<double> stlLocalPerHOf(const ninja::GeometryConfig& gc, const LayersConfig& layers) {
    std::vector<double> v(gc.stls.size(), 0.0);
    for (std::size_t s = 0; s < gc.stls.size(); ++s) {
        auto it = layers.localPerHByRawFile.find(gc.stls[s].rawFile);
        if (it != layers.localPerHByRawFile.end()) v[s] = it->second;
    }
    return v;
}

// --- `localThickness`: the gradient-limited local thickness field --------
//
// t(p) = min over levels L of ( perH * h_L + kThicknessSlope * dist(p, Z_L) ),
// Z_L = the region the dict refines to level >= L, h_L = dx0 / 2^L.
//
// WHY NOT SIMPLY perH * h(cell at p). The offset cut removes everything
// within t(p) of the wall, and a thickness that jumps 2x at every level
// boundary makes that boundary a CLIFF: MEASURED on a hydrofoil, a 1 mm
// level-7 front face had corners 1.1 to 3.1 mm off the wall, the march
// sheared its prisms to 2-20% of the requested height, and every stack along
// the 11->10->9 transitions behind the leading edge was dropped
// (collapsed / meanHeight / bottomFolded + their one-ring neighbours). The
// slope limit turns each cliff into a ramp no steeper than kThicknessSlope,
// the same job snappy's thickness smoothing does.
//
// dist(p, Z_L) comes from the refinement RULES (box / sphere distance,
// max(0, dist(p, S) - band) for a surface band), not from the realized
// octree, so t is a pure function of position -- identical on every rank,
// and one field for the cut (per vertex) and the march (per face).
constexpr double kThicknessSlope = 0.25;

class ThicknessField {
public:
    ThicknessField(const ninja::MeshConfig& cfg, const std::vector<ninja::RefineRegion>& regions,
                   const ninja::GeometryConfig& geom)
        : dx0_((cfg.max.x - cfg.min.x) / static_cast<double>(std::max(1, cfg.nx))), regions_(regions) {
        std::vector<std::string> keys;
        for (const ninja::RefineRegion& r : regions_) {
            maxLevel_ = std::max(maxLevel_, r.level);
            int k = -1;
            if (r.type == "surface") {
                auto it = std::find(keys.begin(), keys.end(), r.stlKey);
                k = static_cast<int>(it - keys.begin());
                if (it == keys.end()) keys.push_back(r.stlKey);
            }
            regionStl_.push_back(k);
        }
        tris_.resize(keys.size());   // sized once: the bins point into it
        for (std::size_t k = 0; k < keys.size(); ++k) {
            for (const std::vector<ninja::StlEntry>* list : {&geom.stls, &geom.sizing}) {
                for (const ninja::StlEntry& se : *list) {
                    if (se.rawFile == keys[k]) tris_[k] = ninja::readStl(se.stlPath);
                }
            }
        }
        for (const std::vector<ninja::Triangle>& t : tris_) bins_.push_back(ninja::buildTriangleAabbBins(t));
    }
    ThicknessField(const ThicknessField&) = delete;
    ThicknessField& operator=(const ThicknessField&) = delete;

    // Total thickness at p for an STL whose thickness is perH x the local cell.
    double at(double perH, const ninja::Vec3& p) const {
        std::vector<double> surfDist(tris_.size(), -1.0);
        std::vector<double> distToLevel(static_cast<std::size_t>(maxLevel_) + 1,
                                        std::numeric_limits<double>::max());
        for (std::size_t r = 0; r < regions_.size(); ++r) {
            const ninja::RefineRegion& reg = regions_[r];
            double d = 0.0;
            if (reg.type == "box") {
                const double dx = std::max({reg.min.x - p.x, 0.0, p.x - reg.max.x});
                const double dy = std::max({reg.min.y - p.y, 0.0, p.y - reg.max.y});
                const double dz = std::max({reg.min.z - p.z, 0.0, p.z - reg.max.z});
                d = std::sqrt(dx * dx + dy * dy + dz * dz);
            } else if (reg.type == "sphere") {
                d = std::max(0.0, ninja::norm(p - reg.centre) - reg.radius);
            } else {
                double& sd = surfDist[static_cast<std::size_t>(regionStl_[r])];
                if (sd < 0.0) sd = std::sqrt(ninja::closestPointOnSoup(bins_[static_cast<std::size_t>(regionStl_[r])], p).distSq);
                d = std::max(0.0, sd - reg.distance);
            }
            double& dl = distToLevel[static_cast<std::size_t>(reg.level)];
            dl = std::min(dl, d);
        }
        double t = perH * dx0_; // level 0 covers everything
        double reach = std::numeric_limits<double>::max();
        for (int L = maxLevel_; L >= 1; --L) {
            reach = std::min(reach, distToLevel[static_cast<std::size_t>(L)]); // region refined to >= L
            if (reach < std::numeric_limits<double>::max()) {
                t = std::min(t, perH * std::ldexp(dx0_, -L) + kThicknessSlope * reach);
            }
        }
        return t;
    }

private:
    double dx0_;
    int maxLevel_ = 0;
    std::vector<ninja::RefineRegion> regions_;
    std::vector<int> regionStl_; // surface rules: index into tris_/bins_, else -1
    std::vector<std::vector<ninja::Triangle>> tris_;
    std::vector<ninja::TriangleAabbBins> bins_;
};

std::vector<std::string> patchNamesOf(const ninja::GeometryConfig& gc); // forward decl, defined below

// --- The MPI gather wire records ---------------------------------------
// ONE packing site, small explicit structs. Points/edges are keyed by
// their GLOBAL fine-grid lattice flat (Refine.hpp's fineLatticeFlat) --
// exact integer keys, no tolerance matching anywhere; values are pure
// functions of position, so duplicates arriving from several ranks are
// bit-identical (asserted at reconstruction: assert, don't exchange).
struct PointRec {
    long long flat = 0;
    int solid = 0;
    int pad = 0;
};
struct EdgeRec {
    long long flatA = 0;
    long long flatB = 0;
    int solidId = 0;
    int pad = 0;
    double x = 0.0, y = 0.0, z = 0.0;
};

std::vector<std::string> patchNamesOf(const ninja::GeometryConfig& gc) {
    std::vector<std::string> names;
    names.reserve(gc.stls.size());
    for (const ninja::StlEntry& se : gc.stls) {
        names.push_back(se.patchName);
    }
    return names;
}

// Runs the march + prism emission over every layered STL
// (thickness > 0) in `geom`, re-homing their wall patches
// and appending prism cells/points to `mesh`/`cellLevel`/`pointLevel`
// in place. A no-op (mesh unchanged) when no STL has a `layers` entry,
// so it is safe to call unconditionally whenever `layers.present`.
// Shared by runPipeline and the `--cut-stats` CLI path so both report
// residual/frozenPoints/droppedFaces consistently.
void runLayersStage(ninja::GeneratedMesh& mesh, std::vector<int>& cellLevel, std::vector<int>& pointLevel,
                     const ninja::MeshConfig& cfg, const ninja::GeometryConfig& geom, const LayersConfig& layers,
                     ninja::CutStats& stats, std::vector<char>* pointFrozenOut = nullptr,
                     std::vector<char>* pointDroppedOut = nullptr,
                     std::vector<int>* cellStackIdOut = nullptr,
                     std::vector<int>* cellLayerIndexOut = nullptr,
                     const ThicknessField* field = nullptr) {
    std::vector<std::vector<ninja::Triangle>> perStlTris = readPerStlTriangles(geom);
    std::vector<double> thickness = stlThicknessOf(geom, layers);
    std::vector<ninja::TriangleAabbBins> perStlBins;
    perStlBins.reserve(perStlTris.size());
    for (const std::vector<ninja::Triangle>& tris : perStlTris) {
        perStlBins.push_back(ninja::buildTriangleAabbBins(tris));
    }
    std::vector<ninja::LayerStlSpec> specs;
    for (std::size_t s = 0; s < geom.stls.size(); ++s) {
        if (thickness[s] > 0.0) {
            ninja::LayerStlSpec spec;
            spec.stlIndex = static_cast<int>(s);
            spec.wallPatchName = geom.stls[s].patchName;
            spec.thickness = thickness[s];
            const std::string& raw = geom.stls[s].rawFile;
            auto nlIt = layers.nLayersByRawFile.find(raw);
            if (nlIt != layers.nLayersByRawFile.end()) spec.nLayers = nlIt->second;
            auto rIt = layers.ratioByRawFile.find(raw);
            if (rIt != layers.ratioByRawFile.end()) spec.expansionRatio = rIt->second;
            auto mfIt = layers.smoothRadiusByRawFile.find(raw);
            if (mfIt != layers.smoothRadiusByRawFile.end()) spec.smoothRadius = mfIt->second;
            auto ltIt = layers.localPerHByRawFile.find(raw);
            if (ltIt != layers.localPerHByRawFile.end()) {
                if (field == nullptr) throw std::logic_error("localThickness without a thickness field");
                const double perH = ltIt->second;
                spec.localThickness = [field, perH](const ninja::Vec3& p) { return field->at(perH, p); };
            }
            specs.push_back(spec);
        }
    }
    if (specs.empty()) {
        return;
    }
    ninja::LayersResult lr = ninja::applyLayers(mesh, cellLevel, pointLevel, cfg, specs, perStlTris, perStlBins,
                                                 layers.debugForceDropSpheres, geom.locationInMesh);
    mesh = std::move(lr.mesh);
    cellLevel = std::move(lr.cellLevel);
    pointLevel = std::move(lr.pointLevel);
    if (pointFrozenOut) *pointFrozenOut = std::move(lr.pointFrozen);
    if (pointDroppedOut) *pointDroppedOut = std::move(lr.pointDropped);
    if (cellStackIdOut) *cellStackIdOut = std::move(lr.cellStackId);
    if (cellLayerIndexOut) *cellLayerIndexOut = std::move(lr.cellLayerIndex);
    stats.residualMean = lr.stats.residualMean;
    stats.residualMax = lr.stats.residualMax;
    stats.frozenPoints = lr.stats.frozenPoints;
    stats.droppedFaces = lr.stats.droppedFaces;
    stats.qualityDroppedFaces = lr.stats.qualityDroppedFaces;
    stats.droppedByReason = lr.stats.droppedByReason;
    stats.surfaceClampedSteps = lr.stats.surfaceClampedSteps;
    stats.buriedLandingFaces = lr.stats.buriedLandingFaces;
    stats.buriedLandingArea = lr.stats.buriedLandingArea;
    stats.buriedLandingMaxDepth = lr.stats.buriedLandingMaxDepth;
    stats.buriedLandingMaxDepthOverH = lr.stats.buriedLandingMaxDepthOverH;
    stats.perStepDropped = lr.stats.perStepDropped;
    stats.perStepFrozen = lr.stats.perStepFrozen;
    stats.perStepPrismCells = lr.stats.perStepPrismCells;
    stats.infeasibleWallArea = lr.stats.infeasibleWallArea;
    stats.smoothGradFallback = lr.stats.smoothGradFallback;
    stats.smoothMovedFaces = lr.stats.smoothMovedFaces;
    stats.smoothMovedArea = lr.stats.smoothMovedArea;
    stats.smoothResidWeightedSum = lr.stats.smoothResidWeightedSum;
    stats.smoothResidAreaSum = lr.stats.smoothResidAreaSum;
    stats.smoothResidMaxOverH = lr.stats.smoothResidMaxOverH;
    stats.smoothSealedRegionCount = static_cast<int>(lr.stats.smoothSealedRegions.size());
    for (const ninja::LayerStats::SealedRegion& region : lr.stats.smoothSealedRegions) {
        stats.smoothSealedRegionArea += region.area;
        // Keep the PER-REGION record, not just the sum. `count` and
        // `area` alone cannot answer the question the sealed-region instrument
        // exists to answer -- WHERE the smoothing changed inside/outside, and
        // whether the same passages seal as the radius grows.
        stats.smoothSealedRegionDetail.push_back({region.count, region.area, region.centroid});
    }
    std::sort(stats.smoothSealedRegionDetail.begin(), stats.smoothSealedRegionDetail.end(),
              [](const ninja::SealedRegionRecord& a, const ninja::SealedRegionRecord& b) { return a.area > b.area; });
    stats.smoothLandingBisectFallback = lr.stats.smoothLandingBisectFallback;
    stats.smoothLandingDegenerateGradFallback = lr.stats.smoothLandingDegenerateGradFallback;
    stats.smoothLandingNoBracketFallback = lr.stats.smoothLandingNoBracketFallback;
    stats.smoothLandingEvalSum = lr.stats.smoothLandingEvalSum;
    stats.smoothLandingEvalMax = lr.stats.smoothLandingEvalMax;
    stats.smoothLandingCandidateCount = lr.stats.smoothLandingCandidateCount;
    stats.stackThicknessBelowHalfT = lr.stats.stackThicknessBelowHalfT;
    stats.stackThicknessBelowQuarterT = lr.stats.stackThicknessBelowQuarterT;
    stats.stackThicknessCount = static_cast<long>(lr.stats.stackThickness.size());
    if (!lr.stats.stackThickness.empty()) {
        std::vector<double> sorted = lr.stats.stackThickness;
        std::sort(sorted.begin(), sorted.end());
        stats.stackThicknessMin = sorted.front();
        double sum = 0.0;
        for (double v : sorted) sum += v;
        stats.stackThicknessMean = sum / static_cast<double>(sorted.size());
        const std::size_t p1idx = static_cast<std::size_t>(0.01 * static_cast<double>(sorted.size() - 1));
        stats.stackThicknessP1 = sorted[p1idx];
    }

    // `stats.wallArea` (as reported by cutMesh) is the PRE-MARCH offset
    // wall's area; the march moves the wall (bottom prism faces +
    // reverted/seam faces), so re-derive it from the FINAL mesh's STL
    // wall patches (the same Newell-normal/2 area convention cutMesh
    // itself uses). Restricted to patch ordinals >= the domain patch
    // count -- a DOMAIN patch may also carry `type wall` (e.g. this
    // suite's "bottom"/"top" box walls), which must NOT be counted as
    // STL wall area.
    double wallArea = 0.0;
    // Domain patches are MERGED by name (several sides may share one,
    // see mergePatchSpecs), so the count of side specs is NOT the count
    // of domain patches on the mesh -- using it here skipped the first
    // STL wall patches and under-reported wallArea.
    const int nDomainPatchesForArea = static_cast<int>(ninja::mergePatchSpecs(cfg.patches).patches.size());
    for (int pIdx = nDomainPatchesForArea; pIdx < static_cast<int>(mesh.patches.size()); ++pIdx) {
        const ninja::PatchInfo& p = mesh.patches[static_cast<std::size_t>(pIdx)];
        for (int i = p.startFace; i < p.startFace + p.nFaces; ++i) {
            const ninja::IntSpan fp = mesh.faces.pointsOf(i);
            ninja::Vec3 nw{0, 0, 0};
            const int n = fp.size();
            for (int k = 0; k < n; ++k) {
                const ninja::Vec3& pc = mesh.points[static_cast<std::size_t>(fp[k])];
                const ninja::Vec3& pn = mesh.points[static_cast<std::size_t>(fp[(k + 1) % n])];
                nw.x += (pc.y - pn.y) * (pc.z + pn.z);
                nw.y += (pc.z - pn.z) * (pc.x + pn.x);
                nw.z += (pc.x - pn.x) * (pc.y + pn.y);
            }
            wallArea += 0.5 * std::sqrt(nw.x * nw.x + nw.y * nw.y + nw.z * nw.z);
        }
    }
    stats.wallArea = wallArea;
}

// Wall-distance fidelity export. Write-only post-pass, run at the very
// end of the pipeline (after writePolyMesh) for EVERY case with a
// `geometry` block, layered or not: for each STL wall patch, the distance
// from every node of that patch's FINAL mesh surface to the ORIGINAL
// (un-offset) STL, via the existing `closestPointOnSoup` machinery
// against a fresh per-STL bins build (independent of whatever offset/
// layers bins the cut stage used -- this is deliberately the TRUE STL,
// not the offset solid). Output: `constant/wallDist_<patch>.vtk` (shared
// DebugVtk legacy-VTK writer, POLYGONS + a `distToStl` POINT_DATA
// scalar); always prints a per-STL mean/max summary line, matching the
// design doc's "always-on stdout summary" requirement regardless of
// whether NINJA_*_DEBUG_DIR is set. Never touches `mesh` -- reads only.
//
// Large distances can appear in easy regions when the offset cut's
// "thin gap seals" rule fires on near-zero-thickness double-sided
// panels (an open shell modelled as two coincident skins).
// `pointFrozen`/`pointDropped` (optional, size-0 when `layers` is
// absent or the caller doesn't have them -- check against
// `mesh.points.size()` before indexing) make that attribution visible
// directly in the file: a node whose `dropped` flag is set sits on a
// face that never completed marching (kept at the offset/intermediate
// position -- exactly the "infeasibleWallArea" mechanism); `frozen`
// marks a node the non-inversion clamp fought to near-zero movement.
void writeWallDistanceExport(const ninja::GeneratedMesh& mesh, const ninja::GeometryConfig& geom,
                              const std::string& caseDir, const std::vector<char>& pointFrozen,
                              const std::vector<char>& pointDropped) {
    if (!geom.present) return;
    const bool haveFrozen = pointFrozen.size() == mesh.points.size();
    const bool haveDropped = pointDropped.size() == mesh.points.size();
    for (const ninja::StlEntry& se : geom.stls) {
        // Locate this STL's own wall patch by name (cutMesh always emits
        // one patch per `wallPatchNames` entry, even if empty -- see
        // Cutter.hpp).
        const ninja::PatchInfo* patch = nullptr;
        for (const ninja::PatchInfo& p : mesh.patches) {
            if (p.name == se.patchName) {
                patch = &p;
                break;
            }
        }
        if (patch == nullptr || patch->nFaces == 0) {
            std::cout << "wallDist " << se.patchName << ": mean=0 max=0 (no faces)\n";
            continue;
        }

        const std::vector<ninja::Triangle> tris = ninja::readStl(se.stlPath);
        const ninja::TriangleAabbBins bins = ninja::buildTriangleAabbBins(tris);

        // Local point remap: only the points this patch's faces actually
        // reference, in first-seen order -- keeps the .vtk small and
        // keeps CellData/PointData indices aligned with DebugVtk's
        // convention (poly indices are LOCAL to `d.pts`).
        std::unordered_map<int, int> globalToLocal;
        std::vector<int> localToGlobal;
        ninja::DebugVtk d;
        for (int fi = patch->startFace; fi < patch->startFace + patch->nFaces; ++fi) {
            const ninja::IntSpan fp = mesh.faces.pointsOf(fi);
            std::vector<int> loop;
            loop.reserve(static_cast<std::size_t>(fp.size()));
            for (int gp : fp) {
                auto it = globalToLocal.find(gp);
                int local;
                if (it == globalToLocal.end()) {
                    local = static_cast<int>(localToGlobal.size());
                    globalToLocal.emplace(gp, local);
                    localToGlobal.push_back(gp);
                } else {
                    local = it->second;
                }
                loop.push_back(local);
            }
            d.polys.push_back(std::move(loop));
        }
        d.pts.reserve(localToGlobal.size());
        std::vector<double> distToStl, frozenArr, droppedArr;
        distToStl.reserve(localToGlobal.size());
        double sumDist = 0.0, maxDist = 0.0;
        for (int gp : localToGlobal) {
            const ninja::Vec3& p = mesh.points[static_cast<std::size_t>(gp)];
            d.pts.push_back(p);
            const ninja::ClosestHit hit = ninja::closestPointOnSoup(bins, p);
            const double dist = std::sqrt(hit.distSq);
            distToStl.push_back(dist);
            sumDist += dist;
            maxDist = std::max(maxDist, dist);
            frozenArr.push_back(haveFrozen && pointFrozen[static_cast<std::size_t>(gp)] ? 1.0 : 0.0);
            droppedArr.push_back(haveDropped && pointDropped[static_cast<std::size_t>(gp)] ? 1.0 : 0.0);
        }
        d.pointScalars = {{"distToStl", distToStl}};
        if (haveFrozen) d.pointScalars.emplace_back("frozen", frozenArr);
        if (haveDropped) d.pointScalars.emplace_back("dropped", droppedArr);
        d.write(caseDir + "/constant/wallDist_" + se.patchName + ".vtk");

        const double meanDist = localToGlobal.empty() ? 0.0 : sumDist / static_cast<double>(localToGlobal.size());
        std::cout << "wallDist " << se.patchName << ": mean=" << meanDist << " max=" << maxDist
                  << " (" << localToGlobal.size() << " nodes)\n";
    }
}

// Per-march-step counters, printed as a space
// separated list on one line ("<name> = a b c d"); an empty vector
// (no `layers` block) prints nothing at all.
void printPerStep(const char* name, const std::vector<int>& v) {
    if (v.empty()) {
        return;
    }
    std::cout << name << " =";
    for (int x : v) std::cout << " " << x;
    std::cout << "\n";
}

// Core-mesh integrity report (NINJA_CORE_INTEGRITY). Zero is the
// contract at every one of these stages; a non-zero count names the
// stage that broke it. See Cutter.hpp's countNonClosedCells.
static void coreIntegrity(const ninja::GeneratedMesh& m, const char* stage) {
    if (!std::getenv("NINJA_CORE_INTEGRITY")) return;
    const int bad = ninja::countNonClosedCells(m);
    std::cout << "coreIntegrity " << stage << ": nCells=" << m.nCells() << " nonClosedCells=" << bad
              << (bad ? "   <-- CONTRACT VIOLATED" : "") << std::endl;
}

// NINJA_VERBOSE=1 enables the long niche stat tail (per-step
// drop/freeze vectors, smoothing internals, stack-thickness block).
// Default output keeps only the lines a user acts on plus everything
// the test/benchmark harnesses grep.
bool verboseStats() {
    const char* v = std::getenv("NINJA_VERBOSE");
    return v != nullptr && std::string(v) == "1";
}

void printCutStats(const ninja::CutStats& s) {
    std::cout << "cut-stats:\n"
              << "  solidVertices = " << s.solidVertices << "\n"
              << "  fluidVertices = " << s.fluidVertices << "\n"
              << "  cutEdges      = " << s.cutEdges << "\n"
              << "  keptCells     = " << s.keptCells << "\n"
              << "  cutCells      = " << s.cutCells << "\n"
              << "  removedCells  = " << s.removedCells << "\n"
              << "  totalVolume   = " << s.totalVolume << "\n"
              << "  wallArea      = " << s.wallArea << "\n"
              << "  multiRootEdges = " << s.multiRootEdges << "\n"
              << "  quantizedIntercepts = " << s.quantizedIntercepts << "\n"
              << "  onVertexIntercepts = " << s.onVertexIntercepts << "\n"
              << "  residualMean  = " << s.residualMean << "\n"
              << "  residualMax   = " << s.residualMax << "\n"
              << "  frozenPoints  = " << s.frozenPoints << "\n"
              << "  droppedFaces  = " << s.droppedFaces << "\n"
              << "  infeasibleWallArea = " << s.infeasibleWallArea << "\n";
    printPerStep("  perStepDropped", s.perStepDropped);
    printPerStep("  perStepFrozen", s.perStepFrozen);
    printPerStep("  perStepPrismCells", s.perStepPrismCells);
    std::cout << "  mergedSliverCells    = " << s.mergedSliverCells << "\n"
              << "  mergeFallbackRemovals = " << s.mergeFallbackRemovals << "\n"
              << "  mergeRefusedTopology = " << s.mergeRefusedTopology << "\n"
              << "  mergeRefusedConvexity = " << s.mergeRefusedConvexity << "\n"
              << "  mergeRefusedGeometry = " << s.mergeRefusedGeometry << "\n"
              << "  mergeReselectedSkew = " << s.mergeReselectedSkew << "\n"
              << "  mergeRefusedSkew = " << s.mergeRefusedSkew << "\n"
              << "  mergeKeptOverSkew = " << s.mergeKeptOverSkew << "\n"
              << "  mergeRefusedWallArea = " << s.mergeRefusedWallArea << "\n"
              << "  disconnectedCellsDropped = " << s.disconnectedCellsDropped << "\n"
              << "  discardedComponents  = " << s.discardedComponents << "\n"
              << "  discardedVolume      = " << s.discardedVolume << "\n";
    std::cout << "  smoothGradFallback = " << s.smoothGradFallback << "\n"
              << "  smoothMovedFaces = " << s.smoothMovedFaces << "\n"
              << "  smoothMovedArea = " << s.smoothMovedArea << "\n"
              << "  smoothResidMean(/h) = "
              << (s.smoothResidAreaSum > 0.0 ? s.smoothResidWeightedSum / s.smoothResidAreaSum : 0.0) << "\n"
              << "  smoothResidMaxOverH = " << s.smoothResidMaxOverH << "\n"
              << "  smoothSealedRegionCount = " << s.smoothSealedRegionCount << "\n"
              << "  smoothSealedRegionArea = " << s.smoothSealedRegionArea << "\n";
    std::cout << "  smoothLandingBisectFallback = " << s.smoothLandingBisectFallback << "\n"
              << "  smoothLandingDegenerateGradFallback = " << s.smoothLandingDegenerateGradFallback << "\n"
              << "  smoothLandingNoBracketFallback = " << s.smoothLandingNoBracketFallback << "\n"
              << "  smoothLandingFallbackRate = "
              << (s.smoothLandingCandidateCount + s.smoothLandingBisectFallback > 0
                      ? static_cast<double>(s.smoothLandingBisectFallback) /
                            static_cast<double>(s.smoothLandingCandidateCount + s.smoothLandingBisectFallback)
                      : 0.0)
              << "  (candidateCount+bisectFallback = "
              << (s.smoothLandingCandidateCount + s.smoothLandingBisectFallback)
              << " total landing attempts, summed over every quality-gate pass)\n"
              << "  smoothLandingEvalMean = "
              << (s.smoothLandingCandidateCount > 0
                      ? static_cast<double>(s.smoothLandingEvalSum) / static_cast<double>(s.smoothLandingCandidateCount)
                      : 0.0)
              << "\n"
              << "  smoothLandingEvalMax = " << s.smoothLandingEvalMax << "\n"
              << "  smoothLandingCandidateCount = " << s.smoothLandingCandidateCount << "\n"
              << "  stackThicknessCount = " << s.stackThicknessCount << "\n"
              << "  stackThicknessMin = " << s.stackThicknessMin << "\n"
              << "  stackThicknessMean = " << s.stackThicknessMean << "\n"
              << "  stackThicknessP1 = " << s.stackThicknessP1 << "\n"
              << "  stackThicknessBelowHalfT = " << s.stackThicknessBelowHalfT << "\n"
              << "  stackThicknessBelowQuarterT = " << s.stackThicknessBelowQuarterT << "\n";
}

void printRefineStats(const ninja::RefineStats& s) {
    std::cout << "refine-stats:\n"
              << "  maxAdjacentLevelDiff = " << s.maxAdjacentLevelDiff << "\n"
              << "  levelCounts =";
    for (std::size_t l = 0; l < s.countByLevel.size(); ++l) {
        std::cout << " " << l << ":" << s.countByLevel[l];
    }
    std::cout << "\n";
}

void printConfig(const ninja::MeshConfig& cfg) {
    std::cout << "domain.min = (" << cfg.min.x << " " << cfg.min.y << " " << cfg.min.z << ")\n";
    std::cout << "domain.max = (" << cfg.max.x << " " << cfg.max.y << " " << cfg.max.z << ")\n";
    std::cout << "domain.n   = (" << cfg.nx << " " << cfg.ny << " " << cfg.nz << ")\n";
    std::cout << "patches:\n";
    for (const ninja::PatchSpec& p : cfg.patches) {
        std::cout << "  " << p.sideKey << ": type=" << p.type << " name=" << p.name << "\n";
    }
}

void printUsage() {
    std::cerr << "usage: ninjaMesher [<caseDir>] | --help | --version | --parse-only <caseDir> | "
                 "--cut-stats <caseDir> | --refine-stats <caseDir>\n";
}

// `--help` is a request, not an error, so it answers on stdout with exit 0 --
// unlike printUsage(), which is the diagnostic printed alongside a bad
// invocation.
int printHelp() {
    std::cout
        << "usage: ninjaMesher [<caseDir>]\n"
           "\n"
           "Meshes an OpenFOAM-layout case directory containing system/ninjaMeshDict\n"
           "and the STL files it references; writes <caseDir>/constant/polyMesh.\n"
           "With no argument, meshes the current directory.\n"
           "\n"
           "options:\n"
           "  --help, -h                this message\n"
           "  --version                 version and copyright\n"
           "  --parse-only <caseDir>    validate the dict and print the config, no meshing\n"
           "  --cut-stats <caseDir>     cut-stage statistics only\n"
           "  --refine-stats <caseDir>  refinement statistics only\n"
           "\n"
           "parallel:\n"
           "  mpirun --bind-to none -np <N> ninjaMesher <caseDir>\n"
           "  (--bind-to none matters: the default core binding pins all OpenMP\n"
           "   threads of a rank to one core and forfeits the threaded speedup)\n"
           "\n"
           "diagnostics (environment, all default off):\n"
           "  NINJA_VERBOSE=1         full statistics tail\n"
           "  NINJA_STAGE_TIMES=1     per-stage wall-clock timing\n"
           "  NINJA_MEM_LOG=1         memory usage logging\n"
           "  NINJA_PRELAYERS_DUMP=<dir>  also write the pre-layer mesh\n"
           "\n"
           "See README.md for the ninjaMeshDict reference and DESIGN.md for how the\n"
           "mesher works.\n";
    return 0;
}

int printVersion(const RunContext& ctx) {
    std::cout << "ninjaMesher version " << kVersion << " (rank " << ctx.rank << "/" << ctx.size
              << ")\n"
              << "Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/\n"
              << "Apache License 2.0\n";
    return 0;
}

// Full pipeline: base mesh -> refinement (if any) -> STL cut (if any).
// Always produces cellLevel/pointLevel alongside the mesh (all-zero
// when refinement is absent: emit always).
struct PipelineResult {
    ninja::GeneratedMesh mesh;
    std::vector<int> cellLevel;
    std::vector<int> pointLevel;
    ninja::CutStats cutStats;
    // Per-point march diagnostics (empty when `layers`
    // is absent, or aligned 1:1 with `mesh.points` otherwise -- see
    // Layers.hpp's LayersResult).
    std::vector<char> pointFrozen;
    std::vector<char> pointDropped;
    // Per-cell stack identity for the layer prisms (-1 elsewhere).
    std::vector<int> cellStackId;
    std::vector<int> cellLayerIndex;
};

PipelineResult runPipeline(const ninja::MeshConfig& cfg, const std::vector<ninja::RefineRegion>& regions,
                            const ninja::GeometryConfig& geom, const LayersConfig& layers,
                            const std::string& caseDir) {
    PipelineResult out;
    // Built only when an STL asks for `localThickness`; read by the offset
    // cut and the march alike.
    std::unique_ptr<ThicknessField> thicknessField;
    if (layers.present && !layers.localPerHByRawFile.empty()) {
        thicknessField = std::make_unique<ThicknessField>(cfg, regions, geom);
    }
    StageTimer timer;
    const bool root = ninja::parRank() == 0;

    // --- Stage 1 (all ranks, block-local): rank-local octree
    // refinement over this rank's base-cell block + halo, with the 2:1
    // grading halo exchange iterated to an Allreduce fixpoint, plus the
    // block-local cut-entity (point/edge) enumeration. R=1 and no
    // marking work when `regions` is empty -- one code path for all N
    // and for both the refined and unrefined pipelines.
    std::unordered_map<std::string, std::vector<ninja::Triangle>> stlTriangles = readStlsByRawFile(geom);
    ninja::LocalRefine lr = ninja::refineLocal(cfg, regions, stlTriangles, geom.present);
    // Peak-RSS lifecycle rule: free consumed inputs at the
    // last-consumer boundary instead of end of pipeline.
    stlTriangles = {}; // only refineLocal reads the per-STL triangle map
    timer.lap("refine-local(parallel)");

    // --- Stage 2 (all ranks, block-local): classification + edge
    // intercepts on this rank's own points/edges only -- the serial
    // CutData functions applied to the block-local entity
    // lists (values are pure functions of (position, replicated STL),
    // identical to what the serial full-mesh sweep computes there).
    std::vector<PointRec> pointRecs;
    std::vector<char> onSurface;
    std::vector<EdgeRec> edgeRecs;
    ninja::OffsetCutStats offsetStats;
    if (geom.present) {
        std::vector<ninja::Triangle> tris = readCombinedStls(geom);
        std::vector<bool> vertexSolid;
        std::unordered_map<ninja::EdgeKey, ninja::Intercept, ninja::EdgeKeyHash> intercepts;
        if (!layers.present) {
            vertexSolid = ninja::classifyVertices(lr.points, tris, geom.locationInMesh);
            intercepts = ninja::computeEdgeIntercepts(lr.points, tris, vertexSolid, lr.edges);
        } else {
            std::vector<std::vector<ninja::Triangle>> perStlTris = readPerStlTriangles(geom);
            std::vector<double> thickness = stlThicknessOf(geom, layers);
            // The half-grid cut. PER-VERTEX band = h_v / 4, with h_v the size
            // of the FINEST leaf touching the vertex
            // (LocalRefine::pointFinestLevel -- a pure function of global
            // position, identical on every rank). h/4 is the Voronoi radius of
            // a half-grid node, i.e. the exact width at which "nearest node"
            // changes answer, so the snap displaces the surface by at most h/4
            // in cells of ANY level. (MEASURED, bm_layers_wfp: one global band
            // sized from the finest level in the case sent almost every
            // crossing on a coarser wall to the edge MIDPOINT -- displacement
            // up to h/2 against layer stacks ~0.65 h thick -- and the
            // collapsed-prism guard then dropped those stacks.)
            // A vertex shared by two levels takes the FINE side's band: one
            // verdict per vertex (watertightness), and the fine edges ending at
            // it keep a true nearest-node rule; the coarse edges there lose a
            // little accuracy, never validity.
            const double hx = (cfg.max.x - cfg.min.x) / cfg.nx;
            const double hy = (cfg.max.y - cfg.min.y) / cfg.ny;
            const double hz = (cfg.max.z - cfg.min.z) / cfg.nz;
            const double h0 = std::min(hx, std::min(hy, hz));
            std::vector<double> bandPerPoint(lr.points.size());
            for (std::size_t i = 0; i < bandPerPoint.size(); ++i) {
                bandPerPoint[i] = std::ldexp(h0, -lr.pointFinestLevel[i]) / 4.0;
            }
            // `localThickness` STLs: the local field at every vertex that can
            // be within reach of the wall (farther than the level-0
            // thickness plus the band, no radius can reach it).
            std::vector<std::vector<double>> pointThickness(geom.stls.size());
            const std::vector<double> perH = stlLocalPerHOf(geom, layers);
            const double dx0 = (cfg.max.x - cfg.min.x) / static_cast<double>(std::max(1, cfg.nx));
            for (std::size_t s = 0; s < geom.stls.size(); ++s) {
                if (perH[s] <= 0.0) continue;
                const ninja::TriangleAabbBins wallBins = ninja::buildTriangleAabbBins(perStlTris[s]);
                std::vector<double>& pt = pointThickness[s];
                pt.assign(lr.points.size(), perH[s] * dx0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
                for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(lr.points.size()); ++i) {
                    const std::size_t k = static_cast<std::size_t>(i);
                    const double reachSq = std::pow(perH[s] * dx0 + bandPerPoint[k], 2);
                    if (ninja::closestPointOnSoup(wallBins, lr.points[k]).distSq > reachSq) continue;
                    pt[k] = thicknessField->at(perH[s], lr.points[k]);
                }
            }
            vertexSolid = ninja::classifyVerticesOffsetHalfGrid(lr.points, tris, perStlTris, geom.locationInMesh,
                                                                thickness, bandPerPoint, onSurface, pointThickness);
            intercepts = ninja::computeEdgeInterceptsOffsetHalfGrid(lr.points, perStlTris, thickness, pointThickness,
                                                                    vertexSolid, onSurface, lr.edges, offsetStats);
        }
        pointRecs.reserve(lr.points.size());
        for (std::size_t i = 0; i < lr.points.size(); ++i) {
            pointRecs.push_back(PointRec{lr.pointFlat[i], vertexSolid[i] ? 1 : 0,
                                          (i < onSurface.size() && onSurface[i]) ? 1 : 0});
        }
        edgeRecs.reserve(intercepts.size());
        for (const ninja::EdgeKey& e : lr.edges) {
            auto it = intercepts.find(e);
            if (it == intercepts.end()) continue;
            edgeRecs.push_back(EdgeRec{lr.pointFlat[static_cast<std::size_t>(e.first)],
                                       lr.pointFlat[static_cast<std::size_t>(e.second)], it->second.solidId, 0,
                                       it->second.point.x, it->second.point.y, it->second.point.z});
        }
        timer.lap("cutdata-local(parallel)");
    }
    // The block-local points/edges/flat keys are fully consumed by
    // Stage 2 (the records now carry everything downstream);
    // ownedLeaves is consumed by the gather just below.
    lr.points = {};
    lr.pointFlat = {};
    lr.edges = {};

    // --- Stage 3: gather to rank 0 (leaves + CutData records). Offset
    // stats counters partition exactly with the uniquely-owned edges,
    // so their SUM equals the np=1 value.
    std::vector<ninja::LeafRec> allLeaves = ninja::parGatherRecordsToRoot(lr.ownedLeaves);
    lr.ownedLeaves = {}; // gathered; only lr.R is read below
    std::vector<PointRec> allPts = ninja::parGatherRecordsToRoot(pointRecs);
    pointRecs = {};
    std::vector<EdgeRec> allEdges = ninja::parGatherRecordsToRoot(edgeRecs);
    edgeRecs = {};
    if (geom.present && layers.present) {
        int sums[3] = {offsetStats.multiRootEdges, offsetStats.quantizedIntercepts,
                       offsetStats.onVertexIntercepts};
        ninja::parAllreduceSumInts(sums, 3);
        offsetStats.multiRootEdges = sums[0];
        offsetStats.quantizedIntercepts = sums[1];
        offsetStats.onVertexIntercepts = sums[2];
    }
    timer.lap("gather");
    if (!root) {
        return out; // non-root ranks are done; only rank 0 holds the full mesh (at the gather, as decided)
    }

    // --- Stage 4 (rank 0): serial assembly. `regions` empty takes the
    // direct generateBaseMesh call (the
    // R=1 lattice flat IS the base-mesh point index, no map needed);
    // refined runs the unchanged serial assembly over the gathered
    // global leaf set.
    ninja::GeneratedMesh preCut;
    std::vector<int> preCellLevel;
    std::vector<int> prePointLevel;
    std::unordered_map<long long, int> flatToId;
    const bool directFlat = regions.empty();
    if (directFlat) {
        preCut = ninja::generateBaseMesh(cfg);
    } else {
        ninja::AssembledRefine ar = ninja::assembleRefined(cfg, allLeaves, lr.R);
        preCut = std::move(ar.mesh);
        preCellLevel = std::move(ar.cellLevel);
        prePointLevel = std::move(ar.pointLevel);
        flatToId.reserve(ar.pointFlat.size());
        for (std::size_t id = 0; id < ar.pointFlat.size(); ++id) {
            flatToId[ar.pointFlat[id]] = static_cast<int>(id);
        }
    }
    allLeaves = {}; // consumed by assembleRefined
    timer.lap("assemble(serial)");
    // Progress line (flushed): a long run must show signs of life.
    std::cout << "refine: " << preCut.nCells() << " cells (R=" << lr.R << ")" << std::endl;

    if (!geom.present) {
        out.mesh = std::move(preCut);
        if (directFlat) {
            out.cellLevel.assign(static_cast<std::size_t>(out.mesh.nCells()), 0);
            out.pointLevel.assign(out.mesh.points.size(), 0);
        } else {
            out.cellLevel = std::move(preCellLevel);
            out.pointLevel = std::move(prePointLevel);
        }
        return out;
    }

    // --- Stage 5 (rank 0): reconstitute the serial CutData from the
    // gathered records (exact integer-key mapping; coverage and
    // duplicate-consistency asserted), then run the UNCHANGED serial
    // tail: cutMesh -> [levels] -> mergeSlivers -> applyLayers ->
    // dropDisconnectedCells.
    auto idOfFlat = [&](long long f) {
        if (directFlat) {
            return static_cast<int>(f);
        }
        auto it = flatToId.find(f);
        if (it == flatToId.end()) {
            throw std::runtime_error("MPI gather error: point key not present in assembled mesh");
        }
        return it->second;
    };
    ninja::CutData cd;
    cd.vertexSolid.assign(preCut.points.size(), false);
    // Carry the ON mask through the gather alongside the binary verdict
    // (PointRec's spare field). Only populated in half-grid mode (the layered
    // cut); an empty vector means "binary path", which Cutter.cpp checks for.
    const bool halfGridOn = layers.present;
    if (halfGridOn) cd.vertexOnSurface.assign(preCut.points.size(), 0);
    std::vector<char> seen(preCut.points.size(), 0);
    for (const PointRec& r : allPts) {
        const int id = idOfFlat(r.flat);
        const bool solid = r.solid != 0;
        if (seen[static_cast<std::size_t>(id)] != 0 &&
            cd.vertexSolid[static_cast<std::size_t>(id)] != solid) {
            throw std::runtime_error("MPI gather error: inconsistent duplicate classification");
        }
        cd.vertexSolid[static_cast<std::size_t>(id)] = solid;
        if (halfGridOn) cd.vertexOnSurface[static_cast<std::size_t>(id)] = static_cast<char>(r.pad != 0 ? 1 : 0);
        seen[static_cast<std::size_t>(id)] = 1;
    }
    for (std::size_t i = 0; i < seen.size(); ++i) {
        if (seen[i] == 0) {
            throw std::runtime_error("MPI gather error: mesh point not covered by any rank");
        }
    }
    for (const EdgeRec& r : allEdges) {
        cd.edgeIntercept[ninja::makeEdgeKey(idOfFlat(r.flatA), idOfFlat(r.flatB))] =
            ninja::Intercept{ninja::Vec3{r.x, r.y, r.z}, r.solidId};
    }
    // The gathered records and the flat->id map end their life
    // here -- CutData now carries the reconstituted truth.
    allPts = {};
    allEdges = {};
    flatToId = {};
    seen = {};
    timer.lap("reconstitute(serial)");

    std::vector<ninja::Triangle> tris = readCombinedStls(geom);
    std::vector<int> originCell;
    std::vector<bool> sliverFlags;
    std::vector<int>* originPtr = directFlat ? nullptr : &originCell;
    if (!layers.present) {
        out.mesh = ninja::cutMesh(preCut, cd, patchNamesOf(geom), tris, out.cutStats, originPtr);
    } else {
        out.mesh = ninja::cutMesh(preCut, cd, patchNamesOf(geom), tris, out.cutStats, originPtr,
                                   layers.absorbVolFrac, /*keepSliversForMerge=*/true, &sliverFlags);
        out.cutStats.multiRootEdges = offsetStats.multiRootEdges;
        out.cutStats.quantizedIntercepts = offsetStats.quantizedIntercepts;
        out.cutStats.onVertexIntercepts = offsetStats.onVertexIntercepts;
    }
    tris = {}; // consumed by cutMesh
    timer.lap("cutmesh(serial)");
    coreIntegrity(out.mesh, "afterCut");
    std::cout << "cut: kept " << out.cutStats.keptCells << " cells, cut " << out.cutStats.cutCells << ", removed "
              << out.cutStats.removedCells << std::endl;
    if (out.cutStats.invalidCoreCells > 0) {
        std::cout << "cut: invalidCoreCells = " << out.cutStats.invalidCoreCells
                  << " (failed the closed-polyhedron post-condition)" << std::endl;
    }
    if (out.cutStats.invalidGeometryCells > 0) {
        std::cout << "cut: invalidGeometryCells = " << out.cutStats.invalidGeometryCells
                  << " (closed, but failed the checkMesh-mimicking face pyramid / tet tests)" << std::endl;
    }
    if (out.cutStats.skewBoundaryCells > 0) {
        std::cout << "cut: skewBoundaryCells = " << out.cutStats.skewBoundaryCells
                  << " (valid, but a face shipping as boundary breached checkMesh's skewness)" << std::endl;
    }
    if (directFlat) {
        out.cellLevel.assign(static_cast<std::size_t>(out.mesh.nCells()), 0);
        out.pointLevel.assign(out.mesh.points.size(), 0);
    } else {
        ninja::PostCutLevels postLevels =
            ninja::propagateLevelsThroughCut(preCut, preCellLevel, prePointLevel, cd, out.mesh, &originCell);
        out.cellLevel = std::move(postLevels.cellLevel);
        out.pointLevel = std::move(postLevels.pointLevel);
    }
    // The pre-cut mesh, CutData and level fields have no
    // consumer past this point -- free them BEFORE the layers march /
    // drop stages instead of at end of pipeline (cutMesh + applyLayers
    // each hold full input + full output; the pre-cut copy dying here
    // removes one of those two peaks).
    preCut = ninja::GeneratedMesh{};
    preCellLevel = {};
    prePointLevel = {};
    cd = ninja::CutData{};
    originCell = {};
    if (layers.present) {
        // Merge sliver cut cells into a kept neighbour BEFORE
        // the march (offset-cut path only).
        ninja::mergeSlivers(out.mesh, out.cellLevel, out.pointLevel, sliverFlags, out.cutStats,
                             geom.locationInMesh);
        std::cout << "slivers: merged " << out.cutStats.mergedSliverCells << ", refused "
                  << (out.cutStats.mergeRefusedTopology + out.cutStats.mergeRefusedConvexity +
                      out.cutStats.mergeRefusedGeometry + out.cutStats.mergeRefusedSkew)
                  << std::endl;
        coreIntegrity(out.mesh, "afterMergeSlivers");
        ninja::auditCellGeometry(out.mesh, "afterMergeSlivers");
        // Planarize the layer top. Runs AFTER mergeSlivers so it sees the
        // final layer-top face set, and BEFORE the pre-layers dump so the
        // artifact shows what the march will actually be handed. See
        // Cutter.hpp for why this is watertightness-safe.
        std::vector<std::string> layeredPatches;
        {
            const std::vector<double> th = stlThicknessOf(geom, layers);
            for (std::size_t si = 0; si < geom.stls.size(); ++si) {
                if (th[si] > 0.0) layeredPatches.push_back(geom.stls[si].patchName);
            }
        }
        // Re-parent coplanar flaps. Runs after mergeSlivers (so it
        // sees the final face set) and BEFORE planarize, because a flap makes
        // its wall face spuriously non-planar -- repairing first is what lets
        // planarize see the genuinely warped faces only. Always on: it is a
        // topology repair with an exact closure identity, not a tolerance.
        // See Cutter.hpp.
        {
            const int nFlaps = ninja::repairCoplanarFlaps(out.mesh, layeredPatches);
            out.cutStats.repairedCoplanarFlaps = nFlaps;
            coreIntegrity(out.mesh, "afterRepairFlaps");
            ninja::auditCellGeometry(out.mesh, "afterRepairFlaps");
        }
        if (!layeredPatches.empty()) {
            // Planarize tolerance: 0.1 * t0, a tenth of the thinnest wall-layer
            // thickness across the layered STLs -- the quantity the layer top
            // has to resolve. The warp distribution is bimodal (numerical noise
            // near zero, genuine warp far above); any value in the gap gives an
            // identical mesh, and 0.1*t0 sits comfortably inside it.
            // Rules based on dx or h_local land past the gap edge and are wrong.
            double minT0 = std::numeric_limits<double>::max();
            for (const auto& kv : layers.thicknessByRawFile) {
                const double t = kv.second;
                if (t <= 0.0) continue;
                const auto nIt = layers.nLayersByRawFile.find(kv.first);
                const auto rIt = layers.ratioByRawFile.find(kv.first);
                const int n = nIt != layers.nLayersByRawFile.end() ? nIt->second : 1;
                const double r = rIt != layers.ratioByRawFile.end() ? rIt->second : 1.0;
                const double t0 = std::fabs(r - 1.0) < 1e-12
                                      ? t / static_cast<double>(n)
                                      : t * (r - 1.0) / (std::pow(r, n) - 1.0);
                minT0 = std::min(minT0, t0);
            }
            const double warpTol = 0.1 * minT0;
            int tris = 0, nonConvex = 0, apexFallback = 0;
            const int nSplit = ninja::planarizeBoundaryFaces(out.mesh, layeredPatches, warpTol, &tris, &nonConvex,
                                                            &apexFallback);
            if (root) {
                std::cout << "planarizeSplitFaces = " << nSplit << "\n"
                          << "planarizeTrianglesAdded = " << tris << "\n"
                          << "planarizeNonConvexSkipped = " << nonConvex << "\n"
                          << "planarizeApexPyramidFallback = " << apexFallback << "\n";
            }
            coreIntegrity(out.mesh, "afterPlanarize");
            ninja::auditCellGeometry(out.mesh, "afterPlanarize");
        }
        // Face-pyramid tripwire on the exact mesh the march is handed
        // (post-merge, post-planarize), so a regression shows up here rather
        // than as condemned layer stacks three stages later. Diagnostic only --
        // the repair lives in cutMesh's coplanar-flap pass; see Cutter.hpp for
        // why flipping faces here is the wrong treatment.
        out.cutStats.badCutFacePyramids = ninja::countBadFacePyramids(out.mesh);
        if (const char* dbg = std::getenv("NINJA_PRELAYERS_DUMP")) {
            ninja::writePolyMesh(out.mesh, dbg);
            // Tag with keptFinal so this pre-layers/pre-drop
            // dump (both fluid components still present -- dropDisconnectedCells
            // hasn't run yet) is one ParaView threshold away from showing
            // only the side dropDisconnectedCells will keep. Redundant
            // flood fill on the CURRENT mesh state, debug-only cost.
            std::vector<int> keptFinal(static_cast<std::size_t>(out.mesh.nCells()), 0);
            const std::vector<char> reached = ninja::reachableCellsFromLocation(out.mesh, geom.locationInMesh);
            for (std::size_t c = 0; c < reached.size(); ++c) keptFinal[c] = reached[c] ? 1 : 0;
            ninja::writeCellIntField(keptFinal, dbg, "keptFinal");
        }
        runLayersStage(out.mesh, out.cellLevel, out.pointLevel, cfg, geom, layers, out.cutStats, &out.pointFrozen,
                       &out.pointDropped, &out.cellStackId, &out.cellLayerIndex, thicknessField.get());
    }
    ninja::auditCellGeometry(out.mesh, "afterLayers");
    ninja::dropDisconnectedCells(out.mesh, out.cellLevel, out.pointLevel, geom.locationInMesh, out.cutStats,
                                 {&out.cellStackId, &out.cellLayerIndex});
    timer.lap("tail(serial)");
    return out;
}

// stderr, one write() per LINE. std::cerr is unbuffered, so every `<<` is its
// own write(); mpirun forwards stdout and stderr through separate pipes and
// can slot a whole stdout line between two fragments. MEASURED, win_desc_slot
// under `ctest -j 12`: "warning: layers: the offset cut (d = t) left
// 22.5524offsetCutDisconnectedComponents = 3" -- a gate line that no longer
// starts its line, so an anchored grep reports it missing. The lock is for the
// OpenMP regions that warn.
class LineAtomicStderr : public std::streambuf {
  public:
    ~LineAtomicStderr() override { flushLine(); }

  protected:
    int_type overflow(int_type ch) override {
        if (ch == traits_type::eof()) return traits_type::not_eof(ch);
        std::lock_guard<std::mutex> lock(mutex_);
        line_.push_back(traits_type::to_char_type(ch));
        if (line_.back() == '\n') flushLine();
        return ch;
    }
    std::streamsize xsputn(const char* s, std::streamsize n) override {
        std::lock_guard<std::mutex> lock(mutex_);
        line_.append(s, static_cast<std::size_t>(n));
        const std::size_t nl = line_.rfind('\n');
        if (nl != std::string::npos) {
            std::string tail = line_.substr(nl + 1);
            line_.resize(nl + 1);
            flushLine();
            line_ = std::move(tail);
        }
        return n;
    }

  private:
    void flushLine() {
        std::size_t done = 0;
        while (done < line_.size()) {
            const ssize_t w = ::write(STDERR_FILENO, line_.data() + done, line_.size() - done);
            if (w <= 0) break;
            done += static_cast<std::size_t>(w);
        }
        line_.clear();
    }
    std::mutex mutex_;
    std::string line_;
};

} // namespace

int main(int argc, char** argv) {
    std::cerr.rdbuf(new LineAtomicStderr); // never freed: std::cerr outlives main()
    MPI_Init(&argc, &argv);

    RunContext ctx;
    MPI_Comm_rank(MPI_COMM_WORLD, &ctx.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ctx.size);

    int exitCode = 0;
    try {
        std::vector<std::string> args(argv + 1, argv + argc);

        if (args.empty()) {
            args.push_back(".");
        }
        if (args[0] == "--help" || args[0] == "-h") {
            exitCode = ctx.rank == 0 ? printHelp() : 0;
        } else if (args[0] == "--version") {
            exitCode = printVersion(ctx);
        } else if (args[0] == "--parse-only") {
            if (args.size() < 2) {
                std::cerr << "error: --parse-only requires a caseDir argument\n";
                printUsage();
                exitCode = 1;
            } else {
                const std::string dictPath = args[1] + "/system/ninjaMeshDict";
                ninja::DictEntry root = ninja::parseDictFile(dictPath);
                ninja::MeshConfig cfg = buildMeshConfig(root);
                printConfig(cfg);
                exitCode = 0;
            }
        } else if (args[0] == "--cut-stats") {
            if (args.size() < 2) {
                std::cerr << "error: --cut-stats requires a caseDir argument\n";
                printUsage();
                exitCode = 1;
            } else {
                const std::string caseDir = args[1];
                const std::string dictPath = caseDir + "/system/ninjaMeshDict";
                ninja::DictEntry root = ninja::parseDictFile(dictPath);
                ninja::MeshConfig cfg = buildMeshConfig(root);
                ninja::GeometryConfig geom = buildGeometryConfig(root, caseDir);
                if (!geom.present) {
                    std::cerr << "error: --cut-stats requires a case with a 'geometry' block\n";
                    exitCode = 1;
                } else {
                    // No regions: this path has always ignored the
                    // `refinement` block, so every cell IS at base
                    // level and level 0 is the honest h for thickness.
                    LayersConfig layers = buildLayersConfig(root, geom, cfg, {});
                    buildLayersDebugConfig(root, layers);
                    // Same collective pipeline as a normal
                    // run, minus refinement (this CLI path has always
                    // ignored the `refinement` block) and minus the
                    // writer; rank 0 prints.
                    PipelineResult res = runPipeline(cfg, {}, geom, layers, caseDir);
                    if (ctx.rank == 0) {
                        printCutStats(res.cutStats);
                    }
                    exitCode = 0;
                }
            }
        } else if (args[0] == "--refine-stats") {
            if (args.size() < 2) {
                std::cerr << "error: --refine-stats requires a caseDir argument\n";
                printUsage();
                exitCode = 1;
            } else {
                const std::string caseDir = args[1];
                const std::string dictPath = caseDir + "/system/ninjaMeshDict";
                ninja::DictEntry root = ninja::parseDictFile(dictPath);
                ninja::MeshConfig cfg = buildMeshConfig(root);
                ninja::GeometryConfig geom = buildGeometryConfig(root, caseDir);
                std::vector<ninja::RefineRegion> regions = buildRefineRegions(root, geom);
                LayersConfig layers = buildLayersConfig(root, geom, cfg, regions);
                buildLayersDebugConfig(root, layers);
                PipelineResult res = runPipeline(cfg, regions, geom, layers, caseDir);
                if (ctx.rank == 0) {
                    ninja::RefineStats stats = ninja::computeRefineStats(res.mesh, res.cellLevel);
                    printRefineStats(stats);
                }
                exitCode = 0;
            }
        } else if (args[0].size() > 2 && args[0][0] == '-' && args[0][1] == '-') {
            std::cerr << "error: unknown option '" << args[0] << "'\n";
            printUsage();
            exitCode = 1;
        } else {
            const std::string caseDir = args[0];
            if (ctx.rank == 0) {
                std::cout << "ninjaMesher " << kVersion << " -- " << caseDir << "\n";
            }
            const std::string dictPath = caseDir + "/system/ninjaMeshDict";
            ninja::DictEntry root = ninja::parseDictFile(dictPath);
            ninja::MeshConfig cfg = buildMeshConfig(root);
            ninja::GeometryConfig geom = buildGeometryConfig(root, caseDir);
            std::vector<ninja::RefineRegion> regions = buildRefineRegions(root, geom);
            LayersConfig layers = buildLayersConfig(root, geom, cfg, regions);
            buildLayersDebugConfig(root, layers);

            PipelineResult res = runPipeline(cfg, regions, geom, layers, caseDir);

            if (ctx.rank == 0) {
            StageTimer wtimer;
            std::cout << "write: " << caseDir << "/constant/polyMesh" << std::endl;
            ninja::writePolyMesh(res.mesh, caseDir);
            // Opt-in (NINJA_LAYER_STACK_FIELDS): write the per-cell
            // stack identity alongside the mesh so ONE stack can be isolated in
            // ParaView (threshold on stackId) and its ACTUAL topology inspected
            // cell by cell. `stackId` is the originating wall seed face --
            // everything one seed extrudes over the whole march shares it --
            // and `layerIndex` is the march step that created the cell; both
            // are -1 for non-prism cells. Aggregate face/cell counts cannot
            // answer "what is this particular face", which is what this is for.
            if (std::getenv("NINJA_LAYER_STACK_FIELDS") && !res.cellStackId.empty()) {
                ninja::writeCellIntField(res.cellStackId, caseDir, "stackId");
                ninja::writeCellIntField(res.cellLayerIndex, caseDir, "layerIndex");
            }
            ninja::writeLevelFields(res.cellLevel, res.pointLevel, caseDir);
            // Opt-in (NINJA_DIAG_FIELDS): visualization-only volScalarFields
            // under 0/ (refinementLevel, isLayer, layerNumber, cellSize).
            if (std::getenv("NINJA_DIAG_FIELDS")) {
                ninja::writeDiagnosticFields(res.mesh, res.cellLevel, res.cellLayerIndex, caseDir);
            }
            wtimer.lap("write(serial)");
            std::cout << "ninjaMesher: wrote " << res.mesh.points.size() << " points, "
                      << res.mesh.nCells() << " cells, " << res.mesh.nFaces() << " faces ("
                      << res.mesh.nInternalFaces << " internal) to " << caseDir << "/constant/polyMesh\n";
            // Write-only post-pass, after the polyMesh
            // write above -- never touches `res.mesh`, so polyMesh bytes
            // are unaffected (see writeWallDistanceExport's own doc).
            writeWallDistanceExport(res.mesh, geom, caseDir, res.pointFrozen, res.pointDropped);
            wtimer.lap("wallDistExport");
            if (geom.present) {
                std::cout << "wallArea = " << res.cutStats.wallArea << "\n";
                std::cout << "disconnectedCellsDropped = " << res.cutStats.disconnectedCellsDropped << "\n";
                std::cout << "discardedComponents = " << res.cutStats.discardedComponents << "\n";
                std::cout << "discardedVolume = " << res.cutStats.discardedVolume << "\n";
            }
            if (layers.present) {
                const bool verbose = verboseStats();
                if (verbose) {
                    std::cout << "multiRootEdges = " << res.cutStats.multiRootEdges << "\n";
                    std::cout << "quantizedIntercepts = " << res.cutStats.quantizedIntercepts << "\n";
                    std::cout << "onVertexIntercepts = " << res.cutStats.onVertexIntercepts << "\n";
                }
                std::cout << "repairedCoplanarFlaps = " << res.cutStats.repairedCoplanarFlaps << "\n";
                std::cout << "badCutFacePyramids = " << res.cutStats.badCutFacePyramids << "\n";
                std::cout << "residualMean = " << res.cutStats.residualMean << "\n";
                std::cout << "residualMax = " << res.cutStats.residualMax << "\n";
                std::cout << "frozenPoints = " << res.cutStats.frozenPoints << "\n";
                std::cout << "droppedFaces = " << res.cutStats.droppedFaces << "\n";
                // The same faces by the first thing that refused them, largest first.
                {
                    std::vector<std::pair<int, int>> byReason;
                    for (std::size_t r = 0; r < res.cutStats.droppedByReason.size(); ++r) {
                        if (res.cutStats.droppedByReason[r] > 0) {
                            byReason.emplace_back(res.cutStats.droppedByReason[r], static_cast<int>(r));
                        }
                    }
                    std::sort(byReason.begin(), byReason.end(), [](const auto& a, const auto& b) {
                        return a.first != b.first ? a.first > b.first : a.second < b.second;
                    });
                    std::cout << "droppedByReason =";
                    for (const auto& [n, r] : byReason) std::cout << " " << ninja::layerDropReasonName(r) << " " << n;
                    std::cout << "\n";
                }
                std::cout << "qualityDroppedFaces = " << res.cutStats.qualityDroppedFaces << "\n";
                std::cout << "surfaceClampedSteps = " << res.cutStats.surfaceClampedSteps << "\n";
                std::cout << "buriedLandingFaces = " << res.cutStats.buriedLandingFaces << "\n";
                std::cout << "buriedLandingArea = " << res.cutStats.buriedLandingArea << "\n";
                std::cout << "buriedLandingMaxDepth = " << res.cutStats.buriedLandingMaxDepth << " ("
                          << res.cutStats.buriedLandingMaxDepthOverH << " h)\n";
                if (res.cutStats.surfaceClampedSteps > 0) {
                    std::cerr << "warning: layers: " << res.cutStats.surfaceClampedSteps
                              << " march step(s) would have carried a front point through a surface into fluid"
                              << " on the far side, and were landed on that surface instead. A handful is a plate"
                              << " thinner than a layer step. Thousands on one wall is a surface the signed distance"
                              << " field cannot see -- two coincident faces (one CAD body lying exactly on another)"
                              << " or an open shell: layers stay incomplete there until the STL is cleaned\n";
                }
                std::cout << "infeasibleWallArea = " << res.cutStats.infeasibleWallArea << "\n";
                if (verbose) {
                    printPerStep("perStepDropped", res.cutStats.perStepDropped);
                    printPerStep("perStepFrozen", res.cutStats.perStepFrozen);
                }
                printPerStep("perStepPrismCells", res.cutStats.perStepPrismCells);
                std::cout << "mergedSliverCells = " << res.cutStats.mergedSliverCells << "\n";
                std::cout << "mergeFallbackRemovals = " << res.cutStats.mergeFallbackRemovals << "\n";
                std::cout << "invalidSliverRemovals = " << res.cutStats.invalidSliverRemovals << "\n";
                std::cout << "mergeRefusedTopology = " << res.cutStats.mergeRefusedTopology << "\n";
                std::cout << "mergeRefusedConvexity = " << res.cutStats.mergeRefusedConvexity << "\n";
                std::cout << "mergeRefusedGeometry = " << res.cutStats.mergeRefusedGeometry << "\n";
                std::cout << "mergeReselectedSkew = " << res.cutStats.mergeReselectedSkew << "\n";
                std::cout << "mergeRefusedSkew = " << res.cutStats.mergeRefusedSkew << "\n";
                std::cout << "mergeKeptOverSkew = " << res.cutStats.mergeKeptOverSkew << "\n";
                std::cout << "mergeRefusedWallArea = " << res.cutStats.mergeRefusedWallArea << "\n";
                // Smoothed-field stats: all zero when smoothRadius is 0.
                if (verbose) {
                    std::cout << "smoothGradFallback = " << res.cutStats.smoothGradFallback << "\n";
                }
                std::cout << "smoothMovedFaces = " << res.cutStats.smoothMovedFaces << "\n";
                std::cout << "smoothMovedArea = " << res.cutStats.smoothMovedArea << "\n";
                if (verbose) {
                    std::cout << "smoothResidMean(/h) = "
                              << (res.cutStats.smoothResidAreaSum > 0.0
                                      ? res.cutStats.smoothResidWeightedSum / res.cutStats.smoothResidAreaSum
                                      : 0.0)
                              << "\n";
                    std::cout << "smoothResidMaxOverH = " << res.cutStats.smoothResidMaxOverH << "\n";
                }
                std::cout << "smoothSealedRegionCount = " << res.cutStats.smoothSealedRegionCount << "\n";
                std::cout << "smoothSealedRegionArea = " << res.cutStats.smoothSealedRegionArea << "\n";
                for (std::size_t ri = 0; ri < res.cutStats.smoothSealedRegionDetail.size(); ++ri) {
                    const ninja::SealedRegionRecord& r = res.cutStats.smoothSealedRegionDetail[ri];
                    std::cout << "smoothSealedRegion[" << ri << "] faces = " << r.count
                              << " area = " << r.area << " centroid = (" << r.centroid.x << " " << r.centroid.y
                              << " " << r.centroid.z << ")\n";
                }
                if (verbose) {
                std::cout << "smoothLandingBisectFallback = " << res.cutStats.smoothLandingBisectFallback << "\n";
                std::cout << "smoothLandingDegenerateGradFallback = "
                          << res.cutStats.smoothLandingDegenerateGradFallback << "\n";
                std::cout << "smoothLandingNoBracketFallback = " << res.cutStats.smoothLandingNoBracketFallback
                          << "\n";
                std::cout << "smoothLandingFallbackRate = "
                          << (res.cutStats.smoothLandingCandidateCount + res.cutStats.smoothLandingBisectFallback > 0
                                  ? static_cast<double>(res.cutStats.smoothLandingBisectFallback) /
                                        static_cast<double>(res.cutStats.smoothLandingCandidateCount +
                                                             res.cutStats.smoothLandingBisectFallback)
                                  : 0.0)
                          << "\n";
                std::cout << "smoothLandingEvalMean = "
                          << (res.cutStats.smoothLandingCandidateCount > 0
                                  ? static_cast<double>(res.cutStats.smoothLandingEvalSum) /
                                        static_cast<double>(res.cutStats.smoothLandingCandidateCount)
                                  : 0.0)
                          << "\n";
                std::cout << "smoothLandingEvalMax = " << res.cutStats.smoothLandingEvalMax << "\n";
                std::cout << "smoothLandingCandidateCount = " << res.cutStats.smoothLandingCandidateCount << "\n";
                std::cout << "stackThicknessCount = " << res.cutStats.stackThicknessCount << "\n";
                std::cout << "stackThicknessMin = " << res.cutStats.stackThicknessMin << "\n";
                std::cout << "stackThicknessMean = " << res.cutStats.stackThicknessMean << "\n";
                std::cout << "stackThicknessP1 = " << res.cutStats.stackThicknessP1 << "\n";
                std::cout << "stackThicknessBelowHalfT = " << res.cutStats.stackThicknessBelowHalfT << "\n";
                std::cout << "stackThicknessBelowQuarterT = " << res.cutStats.stackThicknessBelowQuarterT << "\n";
                } // verbose
            }
            } // ctx.rank == 0
            exitCode = 0;
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        exitCode = 1;
        if (ctx.size > 1) {
            // A rank that errors out must not leave the other
            // ranks blocked in a collective -- abort the whole job.
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    if (std::getenv("NINJA_RSS") != nullptr) {
        // Per-rank peak RSS (stderr, opt-in) -- the
        // memory-scaling measurement; /usr/bin/time -v around mpirun
        // only reports the max across processes (i.e. rank 0's
        // gather-time ceiling), which hides the per-rank reduction.
        struct rusage ru {};
        getrusage(RUSAGE_SELF, &ru);
        std::cerr << "rank " << ctx.rank << " maxrss " << ru.ru_maxrss << " KB\n";
    }
    MPI_Finalize();
    return exitCode;
}
