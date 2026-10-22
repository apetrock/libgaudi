#ifndef GAUDI_ASAWA_ROD_DATUMS_HPP
#define GAUDI_ASAWA_ROD_DATUMS_HPP

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/vec_addendum.h"

#include <memory>

namespace gaudi {
namespace asawa {
namespace rod {

namespace datum_id {
constexpr index_t x = 0;
constexpr index_t v = 1;
constexpr index_t u = 2;
constexpr index_t o = 3;
constexpr index_t M = 4;
constexpr index_t J = 5;
constexpr index_t l0 = 6;
constexpr index_t t = 7;
constexpr index_t rest_l = 8;
} // namespace datum_id

/// Positions: Catmull–Rom on split when a 4-stencil exists.
struct vec3_catmull_datum : public vec3_datum {
  using ptr = std::shared_ptr<vec3_catmull_datum>;
  static ptr create(prim_type type, const std::vector<vec3> &data) {
    return std::make_shared<vec3_catmull_datum>(type, data);
  }
  using vec3_datum::vec3_datum;

  void rod_split(rod & /*R*/, index_t c_prev, index_t c0, index_t c1,
                 index_t c_next, index_t c_new) override {
    ensure_index(c_new);
    if (c_prev < 0 || c_next < 0) {
      __data[c_new] = va::mix(real(0.5), __data[c0], __data[c1]);
    } else {
      __data[c_new] =
          va::catmull_rom(__data[c_prev], __data[c0], __data[c1],
                          __data[c_next], real(0.5));
    }
  }

  // Collapse: inherit midpoint mix from vec3_datum / datum_t.
};

/// Material / angular frames: slerp.
struct quat_slerp_datum : public quat_datum {
  using ptr = std::shared_ptr<quat_slerp_datum>;
  static ptr create(prim_type type, const std::vector<quat> &data) {
    return std::make_shared<quat_slerp_datum>(type, data);
  }
  using quat_datum::quat_datum;

  void calc(const rod & /*R*/, const index_t &i, const index_t & /*prim_id*/,
            const real & /*C*/, const std::vector<index_t> &vals) override {
    if (vals.empty())
      return;
    ensure_index(i);
    __data[i] = __data[vals[0]];
  }

  void rod_split(rod & /*R*/, index_t /*c_prev*/, index_t c0, index_t c1,
                 index_t /*c_next*/, index_t c_new) override {
    ensure_index(c_new);
    __data[c_new] = va::slerp(__data[c0], __data[c1], 0.5f);
  }

  void rod_collapse(rod & /*R*/, index_t c_prev, index_t /*c0*/, index_t c1,
                    index_t /*c_next*/) override {
    if (c_prev < 0)
      return;
    ensure_index(c1);
    __data[c1] = va::slerp(__data[c_prev], __data[c1], 0.5f);
  }
};

/// Rest edge lengths: split proportionally to live geometry (interp_l).
struct real_interp_l_datum : public scalar_datum {
  using ptr = std::shared_ptr<real_interp_l_datum>;
  static ptr create(prim_type type, const std::vector<real> &data) {
    return std::make_shared<real_interp_l_datum>(type, data);
  }
  using scalar_datum::scalar_datum;

  void rod_split(rod &R, index_t /*c_prev*/, index_t c0, index_t c1,
                 index_t /*c_next*/, index_t c_new) override {
    ensure_index(c_new);
    const std::vector<vec3> &px = R.x();
    const vec3 &x0 = px[c0];
    const vec3 &x1 = px[c1];
    const vec3 &xn = px[c_new];
    const real l10 = (x1 - x0).norm();
    if (!(l10 > real(0.0))) {
      __data[c_new] = real(0.0);
      return;
    }
    const real ln0 = (xn - x0).norm();
    const real l1n = (x1 - xn).norm();
    const real l0 = __data[c0];
    __data[c0] = ln0 / l10 * l0;
    __data[c_new] = l1n / l10 * l0;
  }

  void rod_collapse(rod & /*R*/, index_t c_prev, index_t /*c0*/, index_t c1,
                    index_t /*c_next*/) override {
    // Legacy interp_l(c_prev, c1, c1): xn == c1 ⇒ l1n = 0 ⇒ survivor rest = 0.
    // (step() usually rebuilds l0 via update_lengths afterward.)
    if (c_prev < 0)
      return;
    ensure_index(c1);
    __data[c1] = real(0.0);
  }
};

/// Constraint rest length. Split shares the parent rest. Collapse adds the
/// deleted edge's rest onto the surviving previous corner. update_lengths
/// does not write this channel.
struct real_rest_length_datum : public real_interp_l_datum {
  using ptr = std::shared_ptr<real_rest_length_datum>;
  static ptr create(prim_type type, const std::vector<real> &data) {
    return std::make_shared<real_rest_length_datum>(type, data);
  }
  using real_interp_l_datum::real_interp_l_datum;

  void rod_collapse(rod & /*R*/, index_t c_prev, index_t c0, index_t /*c1*/,
                    index_t /*c_next*/) override {
    if (c_prev < 0 || c0 < 0)
      return;
    ensure_index(c_prev);
    ensure_index(c0);
    __data[c_prev] += __data[c0];
  }
};

inline void install_default_datums(rod &R) {
  if (!R.get_data().empty())
    return;
  const size_t n = R.corner_count();
  const prim_type pt = prim_type::CORNER;

  R.insert_datum(vec3_catmull_datum::create(pt, std::vector<vec3>(n, vec3::Zero())));
  R.insert_datum(vec3_datum::create(pt, std::vector<vec3>(n, vec3::Zero())));
  R.insert_datum(quat_slerp_datum::create(
      pt, std::vector<quat>(n, quat::Identity())));
  R.insert_datum(quat_slerp_datum::create(
      pt, std::vector<quat>(n, quat::Identity())));
  R.insert_datum(vec3_datum::create(pt, std::vector<vec3>(n, vec3::Zero())));
  R.insert_datum(vec4_datum::create(pt, std::vector<vec4>(n, vec4::Zero())));
  R.insert_datum(
      real_interp_l_datum::create(pt, std::vector<real>(n, real(0.0))));
  R.insert_datum(scalar_datum::create(pt, std::vector<real>(n, real(0.0))));
  R.insert_datum(
      real_rest_length_datum::create(pt, std::vector<real>(n, real(0.0))));

  auto &ch = R.channels();
  ch.bind_name("x", datum_id::x);
  ch.bind_name("v", datum_id::v);
  ch.bind_name("u", datum_id::u);
  ch.bind_name("o", datum_id::o);
  ch.bind_name("M", datum_id::M);
  ch.bind_name("J", datum_id::J);
  ch.bind_name("l0", datum_id::l0);
  ch.bind_name("t", datum_id::t);
  ch.bind_name("rest_l", datum_id::rest_l);
  ch.set_velocity_slot(datum_id::v);
}

inline std::vector<vec3> &rod_vec3(rod &R, index_t h) {
  return static_pointer_cast<vec3_datum>(R.get_datum(h))->data();
}
inline const std::vector<vec3> &rod_vec3(const rod &R, index_t h) {
  return static_pointer_cast<vec3_datum>(R.const_get_datum(h))->data();
}
inline std::vector<quat> &rod_quat(rod &R, index_t h) {
  return static_pointer_cast<quat_datum>(R.get_datum(h))->data();
}
inline const std::vector<quat> &rod_quat(const rod &R, index_t h) {
  return static_pointer_cast<quat_datum>(R.const_get_datum(h))->data();
}
inline std::vector<real> &rod_real(rod &R, index_t h) {
  return static_pointer_cast<scalar_datum>(R.get_datum(h))->data();
}
inline const std::vector<real> &rod_real(const rod &R, index_t h) {
  return static_pointer_cast<scalar_datum>(R.const_get_datum(h))->data();
}
inline std::vector<vec4> &rod_vec4(rod &R, index_t h) {
  return static_pointer_cast<vec4_datum>(R.get_datum(h))->data();
}
inline const std::vector<vec4> &rod_vec4(const rod &R, index_t h) {
  return static_pointer_cast<vec4_datum>(R.const_get_datum(h))->data();
}

} // namespace rod
} // namespace asawa
} // namespace gaudi

#endif // GAUDI_ASAWA_ROD_DATUMS_HPP
