#ifndef GAUDI_WINDYCHIEN_CATALOG_HPP
#define GAUDI_WINDYCHIEN_CATALOG_HPP

#include "gaudi/paths.hpp"
#include "gaudi/windychien/braid.hpp"

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace gaudi {
namespace windychien {

inline std::vector<braid> load_braid_catalog(const std::string &path_name =
                                                "assets/knots/braids.json") {
  const std::filesystem::path path = gaudi::resolve_path(path_name);
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("windychien: cannot open braid catalog '" +
                             path.string() + "'");
  }

  nlohmann::json root;
  in >> root;
  if (!root.contains("knots") || !root["knots"].is_array()) {
    throw std::runtime_error("windychien: catalog missing 'knots' array");
  }

  std::vector<braid> out;
  out.reserve(root["knots"].size());
  for (const auto &item : root["knots"]) {
    braid b;
    b.name = item.at("name").get<std::string>();
    b.strands = item.at("strands").get<int>();
    b.word = item.at("word").get<std::vector<int>>();
    validate_braid(b);
    out.push_back(std::move(b));
  }
  return out;
}

inline braid find_braid(const std::vector<braid> &catalog,
                        const std::string &name) {
  for (const braid &b : catalog) {
    if (b.name == name) {
      return b;
    }
  }
  throw std::runtime_error("windychien: braid '" + name + "' not in catalog");
}

inline braid load_braid(const std::string &name,
                        const std::string &path_name =
                            "assets/knots/braids.json") {
  return find_braid(load_braid_catalog(path_name), name);
}

} // namespace windychien
} // namespace gaudi

#endif
