#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "gaudi/common.h"

namespace gaudi {
namespace duchamp {

struct mesh_snapshot {
  std::vector<vec3> positions;
  std::vector<uint32_t> indices;
  std::vector<vec3> colors;
};

struct rod_snapshot {
  std::vector<vec3> positions;
  std::vector<vec3> colors;
};

// Trait-style interface most Duchamp demos can satisfy with thin wrappers.
class demo_trait {
public:
  using ptr = std::shared_ptr<demo_trait>;
  virtual ~demo_trait() = default;

  virtual void step(int frame) = 0;
  virtual void reset() = 0;
  virtual int frame_count() const = 0;
  virtual double time_seconds() const = 0;
  virtual std::string name() const = 0;

  virtual std::optional<mesh_snapshot> shell_mesh() const { return std::nullopt; }
  virtual std::optional<rod_snapshot> rod_polyline() const { return std::nullopt; }
};

template <typename DemoT>
class demo_adapter : public demo_trait {
public:
  using demo_ptr = std::shared_ptr<DemoT>;

  static ptr create(demo_ptr demo, std::string demo_name) {
    return std::make_shared<demo_adapter<DemoT>>(std::move(demo), std::move(demo_name));
  }

  demo_adapter(demo_ptr demo, std::string demo_name)
      : _demo(std::move(demo)), _name(std::move(demo_name)) {}

  void step(int frame) override {
    _frame = frame;
    _demo->step(frame);
  }

  void reset() override {
    _frame = 0;
    _demo = DemoT::create();
  }

  int frame_count() const override { return _frame; }
  double time_seconds() const override { return _frame * (1.0 / 60.0); }
  std::string name() const override { return _name; }

protected:
  demo_ptr demo() const { return _demo; }

private:
  demo_ptr _demo;
  std::string _name;
  int _frame = 0;
};

} // namespace duchamp
} // namespace gaudi
