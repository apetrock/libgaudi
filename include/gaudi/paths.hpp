#ifndef GAUDI_PATHS_HPP
#define GAUDI_PATHS_HPP

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace gaudi {
namespace detail {

inline std::filesystem::path strip_assets_prefix(
    const std::filesystem::path &path) {
  auto it = path.begin();
  if (it == path.end() || it->string() != "assets")
    return path;

  std::filesystem::path stripped;
  for (++it; it != path.end(); ++it)
    stripped /= *it;
  return stripped;
}

inline bool fs_exists(const std::filesystem::path &p) {
  std::error_code ec;
  return std::filesystem::exists(p, ec);
}

inline void push_unique(std::vector<std::filesystem::path> &candidates,
                        const std::filesystem::path &candidate) {
  if (candidate.empty())
    return;
  const auto normalized = candidate.lexically_normal();
  if (std::find(candidates.begin(), candidates.end(), normalized) ==
      candidates.end())
    candidates.push_back(normalized);
}

} // namespace detail

inline std::filesystem::path repo_root() {
#ifdef GAUDI_REPO_ROOT
  return std::filesystem::path(GAUDI_REPO_ROOT);
#else
  return {};
#endif
}

/// Resolve a repo-relative path (e.g. "assets/bunny.obj") by searching
/// several candidate locations.  Throws std::runtime_error with a
/// diagnostic listing on failure.
inline std::filesystem::path resolve_path(const std::string &name) {
  const std::filesystem::path requested(name);
  const std::filesystem::path asset_relative =
      detail::strip_assets_prefix(requested);

  std::vector<std::filesystem::path> candidates;

#ifdef GAUDI_REPO_ROOT
  detail::push_unique(candidates,
                      std::filesystem::path(GAUDI_REPO_ROOT) / requested);
#endif

  if (const char *repo_env = std::getenv("GAUDI_REPO_ROOT"))
    detail::push_unique(candidates,
                        std::filesystem::path(repo_env) / requested);

  if (const char *asset_env = std::getenv("GAUDI_ASSET_ROOT"))
    detail::push_unique(candidates,
                        std::filesystem::path(asset_env) / asset_relative);

  detail::push_unique(candidates, requested);
  detail::push_unique(candidates,
                      std::filesystem::current_path() / requested);

  for (const auto &c : candidates) {
    if (detail::fs_exists(c))
      return c;
  }

  std::ostringstream msg;
  msg << "gaudi::resolve_path failed for '" << name << "'. Searched:";
  for (const auto &c : candidates)
    msg << "\n  - " << c.string();
#ifdef GAUDI_REPO_ROOT
  msg << "\n  GAUDI_REPO_ROOT (compile) = \"" << GAUDI_REPO_ROOT << "\"";
#else
  msg << "\n  GAUDI_REPO_ROOT (compile) = (not defined)";
#endif
  if (const char *repo_env = std::getenv("GAUDI_REPO_ROOT"))
    msg << "\n  GAUDI_REPO_ROOT (env) = \"" << repo_env << "\"";
  else
    msg << "\n  GAUDI_REPO_ROOT (env) = (not set)";
  if (const char *asset_env = std::getenv("GAUDI_ASSET_ROOT"))
    msg << "\n  GAUDI_ASSET_ROOT (env) = \"" << asset_env << "\"";
  else
    msg << "\n  GAUDI_ASSET_ROOT (env) = (not set)";
  throw std::runtime_error(msg.str());
}

/// Non-throwing variant: returns fallback on resolution failure.
inline std::filesystem::path resolve_path_or(
    const std::string &name, const std::filesystem::path &fallback) {
  try {
    return resolve_path(name);
  } catch (...) {
    return fallback;
  }
}

} // namespace gaudi

#endif // GAUDI_PATHS_HPP
