// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <map>
#include <string>
#include <vector>

namespace ninja {

// A minimal recursive-descent parser for the subset of OpenFOAM dict
// grammar used by ninjaMeshDict: nested braces, "key value;",
// "key (list);", "//" line comments and "/* */" block comments. The
// optional FoamFile header block (if present) is parsed like any other
// sub-dict and simply ignored by callers.
//
// A DictEntry is either a scalar string value, a parenthesised list of
// tokens (stored as strings, split on whitespace), or a nested
// sub-dictionary. Exactly one of these is populated.
struct DictEntry {
    enum class Kind { Scalar, List, Dict };
    Kind kind = Kind::Scalar;

    std::string scalar;
    std::vector<std::string> list;
    std::map<std::string, DictEntry> dict;
    // Sub-entry keys in source declaration order (Kind::Dict only). The
    // `dict` map above is alphabetical (std::map), so callers that need
    // declaration order (e.g. multi-STL patch ordering) must use
    // this instead of iterating `dict` directly.
    std::vector<std::string> order;
};

// Parses OpenFOAM-dict-formatted text into a top-level DictEntry of
// Kind::Dict. Throws std::runtime_error with a descriptive message on
// malformed input.
DictEntry parseDict(const std::string& text);

// Reads the file at path and parses it. Throws std::runtime_error on
// I/O failure or parse failure.
DictEntry parseDictFile(const std::string& path);

// Convenience accessors on a Kind::Dict entry. Throw std::runtime_error
// (with the missing key name in the message) if the key is absent or
// has the wrong kind.
const DictEntry& requireEntry(const DictEntry& d, const std::string& key);
std::string requireScalar(const DictEntry& d, const std::string& key);
std::vector<std::string> requireList(const DictEntry& d, const std::string& key);

} // namespace ninja
