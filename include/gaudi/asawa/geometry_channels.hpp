#ifndef GAUDI_ASAWA_GEOMETRY_CHANNELS_HPP
#define GAUDI_ASAWA_GEOMETRY_CHANNELS_HPP

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "gaudi/common.h"

namespace gaudi {
namespace asawa {

class datum;

/// Attribute pool adapter: owns datums; topology stays on shell/rod.
/// Datum handles are views (slot_id + data()). Sugar for known roles; string
/// map for app-specific names. Safety is API convention only.
class GeometryChannels {
public:
  using datum_ptr = std::shared_ptr<datum>;

  index_t insert(datum_ptr d) {
    pool_.push_back(d);
    return static_cast<index_t>(pool_.size() - 1);
  }

  index_t insert(const std::string &name, datum_ptr d) {
    const index_t id = insert(d);
    by_name_[name] = id;
    return id;
  }

  void bind_name(const std::string &name, index_t slot) {
    by_name_[name] = slot;
  }

  void set_velocity_slot(index_t slot) { velocity_slot_ = slot; }
  index_t velocity_slot() const { return velocity_slot_; }

  datum_ptr &get(index_t i) { return pool_[static_cast<size_t>(i)]; }
  const datum_ptr &get(index_t i) const {
    return pool_[static_cast<size_t>(i)];
  }

  /// App-specific / named channel lookup (nullptr if missing).
  datum_ptr get(const std::string &name) const { return get_named(name); }

  datum_ptr get_named(const std::string &name) const {
    auto it = by_name_.find(name);
    if (it == by_name_.end())
      return nullptr;
    return pool_[static_cast<size_t>(it->second)];
  }

  index_t named_slot(const std::string &name) const {
    auto it = by_name_.find(name);
    if (it == by_name_.end())
      return -1;
    return it->second;
  }

  std::vector<datum_ptr> &pool() { return pool_; }
  const std::vector<datum_ptr> &pool() const { return pool_; }

  bool empty() const { return pool_.empty(); }
  size_t size() const { return pool_.size(); }

private:
  std::vector<datum_ptr> pool_;
  std::unordered_map<std::string, index_t> by_name_;
  index_t velocity_slot_ = -1;
};

} // namespace asawa
} // namespace gaudi

#endif
