// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#include "Stl.hpp"

#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ninja {

namespace {

std::vector<Triangle> readBinaryStl(const std::vector<char>& buf) {
    std::vector<Triangle> tris;
    uint32_t count = 0;
    std::memcpy(&count, buf.data() + 80, 4);
    tris.reserve(count);
    std::size_t off = 84;
    auto readVec = [&](std::size_t o) {
        float f[3];
        std::memcpy(f, buf.data() + o, 12);
        return Vec3{static_cast<double>(f[0]), static_cast<double>(f[1]), static_cast<double>(f[2])};
    };
    for (uint32_t t = 0; t < count; ++t) {
        // normal (12 bytes, ignored), v0,v1,v2 (12 bytes each), attr (2 bytes)
        Triangle tri;
        tri.v0 = readVec(off + 12);
        tri.v1 = readVec(off + 24);
        tri.v2 = readVec(off + 36);
        tris.push_back(tri);
        off += 50;
    }
    return tris;
}

std::vector<Triangle> readAsciiStl(const std::string& text) {
    std::vector<Triangle> tris;
    std::istringstream in(text);
    std::string tok;
    std::vector<Vec3> pending;
    while (in >> tok) {
        if (tok == "vertex") {
            double x, y, z;
            in >> x >> y >> z;
            pending.push_back(Vec3{x, y, z});
            if (pending.size() == 3) {
                tris.push_back(Triangle{pending[0], pending[1], pending[2]});
                pending.clear();
            }
        }
    }
    return tris;
}

} // namespace

std::vector<Triangle> readStl(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Stl error: cannot open file '" + path + "'");
    }
    std::vector<char> buf((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    if (buf.size() >= 84) {
        uint32_t count = 0;
        std::memcpy(&count, buf.data() + 80, 4);
        const std::size_t expected = 84 + static_cast<std::size_t>(count) * 50;
        if (expected == buf.size()) {
            return readBinaryStl(buf);
        }
    }
    return readAsciiStl(std::string(buf.begin(), buf.end()));
}

} // namespace ninja
