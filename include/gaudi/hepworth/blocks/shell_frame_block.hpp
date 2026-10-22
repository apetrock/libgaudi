#ifndef __GAUDI_HEPWORTH_SHELL_FRAME_BLOCK__
#define __GAUDI_HEPWORTH_SHELL_FRAME_BLOCK__

#include <cmath>
#include <queue>
#include <vector>

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/blocks/block_wrapper.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// One quaternion per shell vertex. Columns are red, green, blue. The datum
/// is VERTEX, so subdivision slerps and pack permutes. Validity is a separate
/// real: a zero quaternion is Identity, so norm alone cannot mark a new slot.
class shell_frame_block : public dof_block_wrapper {
public:
  using ptr = std::shared_ptr<shell_frame_block>;

  explicit shell_frame_block(asawa::shell::shell::ptr mesh)
      : mesh(std::move(mesh)) {
    _datum_id = asawa::init_vert_datum<quat>(*this->mesh, quat::Identity());
    _valid_id = asawa::init_vert_datum<real>(*this->mesh, real(0));
  }

  asawa::shell::shell::ptr mesh;
  std::vector<vec4> J;
  std::vector<quat> o;

  const std::vector<char> &valid() const { return _valid; }

  std::vector<quat> &frames() {
    return static_pointer_cast<asawa::datum_t<quat>>(mesh->get_datum(_datum_id))
        ->data();
  }

  std::vector<real> &valid_weights() {
    return asawa::get_real_data(*mesh, _valid_id);
  }

  void prepare(solver_context &ctx) override {
    (void)ctx;
    std::vector<quat> &u = frames();
    const int n = static_cast<int>(mesh->vert_count());
    if (static_cast<int>(u.size()) < n)
      u.resize(static_cast<size_t>(n), quat::Identity());
    else if (static_cast<int>(u.size()) > n)
      u.resize(static_cast<size_t>(n));

    const std::vector<vec3> &x = asawa::const_get_vec_data(*mesh, 0);
    ensure_frames(x);

    o.assign(static_cast<size_t>(n), quat(0, 0, 0, 0));
    const std::vector<real> area = asawa::shell::vertex_areas(*mesh, x);
    J.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
      const real a =
          i < static_cast<int>(area.size()) ? area[static_cast<size_t>(i)] : 0.0;
      J[static_cast<size_t>(i)] = vec4(a, a, a, a);
    }
  }

  void apply_external_input(const std::vector<vec3> &input, real scale = 1.0) override {
    (void)input;
    (void)scale;
  }

  void flush_external_inputs() override {}

  void emit_dof_blocks(std::vector<sim_block::ptr> &out) override {
    _quat = quat_block::create(J, frames(), o);
    out.push_back(_quat);
  }

  quat_block::ptr dof_quat() const { return _quat; }

private:
  static bool live(const asawa::shell::shell &M, int i) {
    if (i < 0 || i >= static_cast<int>(M.vert_count()))
      return false;
    return static_cast<int>(M.vbegin(asawa::shell::vert_id(i))) >= 0;
  }

  bool fit_longest_edge(int i, const std::vector<vec3> &x, quat &q) const {
    if (!live(*mesh, i) || i >= static_cast<int>(x.size()))
      return false;
    const asawa::shell::VertId vi = asawa::shell::vert_id(i);
    const std::vector<asawa::shell::VertId> ring = mesh->get_one_ring(vi);
    if (ring.size() < 3)
      return false;
    if (asawa::shell::vert_area(*mesh, vi, x) < 1e-8)
      return false;
    const vec3 N = asawa::shell::vert_normal(*mesh, vi, x);
    if (!N.allFinite() || N.squaredNorm() < 1e-12)
      return false;
    real best = 0;
    vec3 T = vec3::Zero();
    const vec3 xi = x[static_cast<size_t>(i)];
    for (asawa::shell::VertId vj : ring) {
      const int j = static_cast<int>(vj);
      if (j < 0 || j >= static_cast<int>(x.size()))
        continue;
      vec3 e = x[static_cast<size_t>(j)] - xi;
      e -= e.dot(N) * N;
      const real l2 = e.squaredNorm();
      if (l2 > best) {
        best = l2;
        T = e;
      }
    }
    if (!(best > 1e-16) || !T.allFinite())
      return false;
    T.normalize();
    vec3 B = N.cross(T);
    if (B.squaredNorm() < 1e-24)
      return false;
    B.normalize();
    mat3 R;
    R.col(0) = T;
    R.col(1) = B;
    R.col(2) = N.normalized();
    q = quat(R);
    q.normalize();
    const real n2 = q.squaredNorm();
    return std::isfinite(n2) && n2 > 0.25;
  }

  void ensure_frames(const std::vector<vec3> &x) {
    std::vector<quat> &u = frames();
    std::vector<real> &vw = valid_weights();
    const int n = static_cast<int>(u.size());
    if (static_cast<int>(vw.size()) < n)
      vw.resize(static_cast<size_t>(n), real(0));
    else if (static_cast<int>(vw.size()) > n)
      vw.resize(static_cast<size_t>(n));
    _valid.assign(static_cast<size_t>(n), 0);
    std::vector<char> fitted(static_cast<size_t>(n), 0);

    for (int i = 0; i < n; ++i) {
      if (!live(*mesh, i)) {
        vw[static_cast<size_t>(i)] = 0.0;
        continue;
      }
      quat &qi = u[static_cast<size_t>(i)];
      const real n2 = qi.squaredNorm();
      if (vw[static_cast<size_t>(i)] >= 0.5 && n2 >= 0.25 && std::isfinite(n2)) {
        qi.normalize();
        fitted[static_cast<size_t>(i)] = 1;
        continue;
      }
      quat q = quat::Identity();
      if (fit_longest_edge(i, x, q)) {
        qi = q;
        fitted[static_cast<size_t>(i)] = 1;
      }
    }

    std::queue<int> grow;
    for (int i = 0; i < n; ++i)
      if (fitted[static_cast<size_t>(i)])
        grow.push(i);
    while (!grow.empty()) {
      const int i = grow.front();
      grow.pop();
      for (asawa::shell::VertId vj :
           mesh->get_one_ring(asawa::shell::vert_id(i))) {
        const int j = static_cast<int>(vj);
        if (j < 0 || j >= n || fitted[static_cast<size_t>(j)] || !live(*mesh, j))
          continue;
        u[static_cast<size_t>(j)] = u[static_cast<size_t>(i)];
        fitted[static_cast<size_t>(j)] = 1;
        grow.push(j);
      }
    }

    std::vector<char> seen(static_cast<size_t>(n), 0);
    for (int seed = 0; seed < n; ++seed) {
      if (!fitted[static_cast<size_t>(seed)] || seen[static_cast<size_t>(seed)])
        continue;
      std::queue<int> align;
      seen[static_cast<size_t>(seed)] = 1;
      align.push(seed);
      while (!align.empty()) {
        const int i = align.front();
        align.pop();
        for (asawa::shell::VertId vj :
             mesh->get_one_ring(asawa::shell::vert_id(i))) {
          const int j = static_cast<int>(vj);
          if (j < 0 || j >= n || !fitted[static_cast<size_t>(j)] ||
              seen[static_cast<size_t>(j)])
            continue;
          seen[static_cast<size_t>(j)] = 1;
          if (u[static_cast<size_t>(i)].coeffs().dot(
                  u[static_cast<size_t>(j)].coeffs()) < 0.0)
            u[static_cast<size_t>(j)].coeffs() *= real(-1);
          align.push(j);
        }
      }
    }

    for (int i = 0; i < n; ++i) {
      quat &qi = u[static_cast<size_t>(i)];
      if (!fitted[static_cast<size_t>(i)] || !live(*mesh, i)) {
        const real n2 = qi.squaredNorm();
        if (!(n2 >= 0.25) || !std::isfinite(n2))
          qi = quat::Identity();
        vw[static_cast<size_t>(i)] = 0.0;
        continue;
      }
      qi.normalize();
      vw[static_cast<size_t>(i)] = 1.0;
      _valid[static_cast<size_t>(i)] = 1;
    }
  }

  index_t _datum_id = -1;
  index_t _valid_id = -1;
  std::vector<char> _valid;
  quat_block::ptr _quat;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
