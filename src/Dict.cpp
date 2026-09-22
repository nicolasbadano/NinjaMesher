// NinjaMesher -- hex-dominant cut-cell mesher for CFD
// Copyright 2026 Nicolás Diego Badano -- https://hydronumerical.com/nicolasbadano/
// SPDX-License-Identifier: Apache-2.0

#include "Dict.hpp"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ninja {

namespace {

// Strips "//" line comments and "/* */" block comments, leaving
// whitespace layout otherwise intact so error messages could later be
// extended with line numbers if needed.
std::string stripComments(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    for (std::size_t i = 0; i < text.size();) {
        if (text[i] == '/' && i + 1 < text.size() && text[i + 1] == '/') {
            while (i < text.size() && text[i] != '\n') {
                ++i;
            }
        } else if (text[i] == '/' && i + 1 < text.size() && text[i + 1] == '*') {
            i += 2;
            while (i + 1 < text.size() && !(text[i] == '*' && text[i + 1] == '/')) {
                ++i;
            }
            i += 2;
        } else {
            out.push_back(text[i]);
            ++i;
        }
    }
    return out;
}

class Tokenizer {
public:
    explicit Tokenizer(const std::string& text) : text_(text) {}

    // Returns the next token: one of '{', '}', '(', ')', ';', or a bare
    // word. Returns an empty string at end of input.
    std::string next() {
        skipWhitespace();
        if (pos_ >= text_.size()) {
            return {};
        }
        const char c = text_[pos_];
        if (c == '{' || c == '}' || c == '(' || c == ')' || c == ';') {
            ++pos_;
            return std::string(1, c);
        }
        std::size_t start = pos_;
        while (pos_ < text_.size() && !std::isspace(static_cast<unsigned char>(text_[pos_])) &&
               text_[pos_] != '{' && text_[pos_] != '}' && text_[pos_] != '(' &&
               text_[pos_] != ')' && text_[pos_] != ';') {
            ++pos_;
        }
        return text_.substr(start, pos_ - start);
    }

    // Reads raw tokens up to (and consuming) the matching ')' for a
    // list body; the caller has already consumed the opening '('.
    std::vector<std::string> readListBody() {
        std::vector<std::string> items;
        for (;;) {
            skipWhitespace();
            if (pos_ >= text_.size()) {
                throw std::runtime_error("Dict parse error: unterminated list (missing ')')");
            }
            if (text_[pos_] == ')') {
                ++pos_;
                break;
            }
            std::size_t start = pos_;
            while (pos_ < text_.size() && !std::isspace(static_cast<unsigned char>(text_[pos_])) &&
                   text_[pos_] != ')') {
                ++pos_;
            }
            items.push_back(text_.substr(start, pos_ - start));
        }
        return items;
    }

private:
    void skipWhitespace() {
        while (pos_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[pos_]))) {
            ++pos_;
        }
    }

    const std::string& text_;
    std::size_t pos_ = 0;
};

// Parses entries until a '}' (consumed) or end of input (top level).
std::map<std::string, DictEntry> parseEntries(Tokenizer& tok, bool topLevel, std::vector<std::string>& order) {
    std::map<std::string, DictEntry> entries;
    for (;;) {
        std::string key = tok.next();
        if (key.empty()) {
            if (!topLevel) {
                throw std::runtime_error("Dict parse error: unterminated dict (missing '}')");
            }
            break;
        }
        if (key == "}") {
            if (topLevel) {
                throw std::runtime_error("Dict parse error: unexpected '}' at top level");
            }
            break;
        }

        std::string next = tok.next();
        if (next.empty()) {
            throw std::runtime_error("Dict parse error: expected value after key '" + key + "'");
        }

        DictEntry entry;
        if (next == "{") {
            entry.kind = DictEntry::Kind::Dict;
            entry.dict = parseEntries(tok, false, entry.order);
        } else if (next == "(") {
            entry.kind = DictEntry::Kind::List;
            entry.list = tok.readListBody();
            std::string semi = tok.next();
            if (semi != ";") {
                throw std::runtime_error("Dict parse error: expected ';' after list value for key '" +
                                          key + "'");
            }
        } else {
            // Scalar value: collect tokens until ';' (values are
            // normally a single word, but be lenient and join with
            // spaces in case of e.g. "type patch").
            std::string scalar = next;
            std::string t = tok.next();
            while (t != ";") {
                if (t.empty()) {
                    throw std::runtime_error("Dict parse error: expected ';' after value for key '" +
                                              key + "'");
                }
                scalar += " " + t;
                t = tok.next();
            }
            entry.kind = DictEntry::Kind::Scalar;
            entry.scalar = scalar;
        }

        // A REPEATED KEY IS AN ERROR, not an override. `entries[key] =
        // ...` silently discarded the earlier entry while
        // `order.push_back` recorded the key twice, so the two halves
        // of a DictEntry disagreed and the loss was invisible.
        //
        // MEASURED (Tests/cases/dict_dupkey): two refinement regions
        // both named `dup`, level 2 declared first and level 1 second,
        // produced 7 672 cells -- bit-identical to the level-1 region
        // ALONE, against the 60 648 the same two regions give under
        // distinct names, where they combine by max. Nothing warned. The
        // exposure is worse outside `refinement`: buildGeometryConfig
        // iterates `order` for patch ordering (the map itself is
        // alphabetical and must not be), so a repeated geometry key
        // would emit a patch twice from one STL entry.
        //
        // No dict in the tree has a duplicate sub-dict key, so this
        // refuses only inputs that were already being silently
        // mis-read.
        if (entries.find(key) != entries.end()) {
            throw std::runtime_error("Dict parse error: duplicate key '" + key +
                                      "' in the same block -- each key may appear only once (an earlier"
                                      " entry would otherwise be silently discarded)");
        }
        entries[key] = std::move(entry);
        order.push_back(key);
    }
    return entries;
}

} // namespace

DictEntry parseDict(const std::string& text) {
    const std::string cleaned = stripComments(text);
    Tokenizer tok(cleaned);
    DictEntry root;
    root.kind = DictEntry::Kind::Dict;
    root.dict = parseEntries(tok, /*topLevel=*/true, root.order);
    return root;
}

DictEntry parseDictFile(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Dict parse error: cannot open file '" + path + "'");
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    return parseDict(ss.str());
}

const DictEntry& requireEntry(const DictEntry& d, const std::string& key) {
    auto it = d.dict.find(key);
    if (it == d.dict.end()) {
        throw std::runtime_error("Dict error: missing required key '" + key + "'");
    }
    return it->second;
}

std::string requireScalar(const DictEntry& d, const std::string& key) {
    const DictEntry& e = requireEntry(d, key);
    if (e.kind != DictEntry::Kind::Scalar) {
        throw std::runtime_error("Dict error: key '" + key + "' is not a scalar value");
    }
    return e.scalar;
}

std::vector<std::string> requireList(const DictEntry& d, const std::string& key) {
    const DictEntry& e = requireEntry(d, key);
    if (e.kind != DictEntry::Kind::List) {
        throw std::runtime_error("Dict error: key '" + key + "' is not a list value");
    }
    return e.list;
}

} // namespace ninja
