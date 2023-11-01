
#ifndef __CCD_MODULE__
#define __CCD_MODULE__

#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/common.h"
#include "module_base_shell.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace gaudi {
namespace duchamp {

/////////////////////////////////////////////////////////////

class continuous_collision_detection : public module_base_shell {
public:
  DEFINE_CREATE_FUNC(continuous_collision_detection)
  continuous_collision_detection(asawa::shell::shell::ptr M,
                                 asawa::shell::dynamic::ptr D, real eps)
      : module_base_shell(M), __dynamic(D), _eps(eps) {

    std::vector<vec3> &x = asawa::get_vec_data(*_M, 0);

    _i = gaudi::asawa::init_vert_datum<vec3>(*_M, vec3::Zero());
  };

  virtual ~continuous_collision_detection(){};

  void init_step(real h) {
    std::vector<vec3> &f = asawa::get_vec_data(*_M, _i);
    f = std::vector<vec3>(f.size(), vec3::Zero());
  }

  void add_shell_force(const std::vector<vec3> &fs, const real &h = 1.0) {
    std::vector<vec3> &f = asawa::get_vec_data(*_M, _i);

    for (int i = 0; i < fs.size(); i++) {
      f[i] += h * fs[i];
    }
  }

  void clip_outliers(std::vector<vec3> &v, real clip) {
    // takes in velocity vector, calc mean/std of norm, then iterates through
    // vector and clips outliers based on mean/std

    std::array<real, 2> mean_std = va::mean_std<vec3, real>(v);
    real mean = mean_std[0];
    real stddev = mean_std[1];
    std::cout << "clip mean/std: " << mean << " " << stddev << std::endl;
    //

    if (stddev < 0.5 * mean)
      return;

    std::vector<vec3> v_clip(v.size(), vec3::Zero());
    for (int i = 0; i < v.size(); i++) {
      real v_norm = v[i].norm();
      if (v_norm > mean + clip * stddev) {
        v[i] = (mean + clip * stddev) * v[i].normalized();
      }
    }
  }

  virtual real calc_dt(const std::vector<vec3> &v0, real h0) {
    const std::vector<vec3> &x = asawa::get_vec_data(*_M, 0);
    std::vector<vec3> v0f = asawa::shell::vert_to_face<vec3>(*_M, x, v0);
    std::vector<vec3> v1 = calder::mls_avg(*_M, v0f, x, 2.0 * _eps, 2.0);
    std::vector<real> h_vec(v0.size(), 0.0);
    real hmin = std::numeric_limits<real>::infinity();
    for (int i = 0; i < v0.size(); i++) {
      vec3 v0i = v0[i];
      vec3 v1i = v1[i];
      // project v0i onto v1i

      if (v0i.norm() < 1e-5)
        h_vec[i] = 1.0;

      vec3 vproj = (v1i).dot(v0i) / v0i.dot(v0i) * v0i;
      real vproj_norm = vproj.norm();
      real v0i_norm = v0i.norm();
      real vratio = vproj_norm / v0i_norm;
      // todo: why is vratio nan?
      if (isnan(vratio))
        continue;
      vratio = std::min(vratio, 1.0);
      real hi = vratio * h0;

      h_vec[i] = vratio;
      hmin = std::min(hmin, hi);
    }

    std::array<real, 2> mean_std = va::mean_std<real, real>(h_vec);
    real h_mean = mean_std[0];
    real h_stddev = mean_std[1];

    real h = (h_mean - 0.5 * h_stddev) * h0;
    if (h_stddev > h_mean)
      h = h0;

    return h;
  }

  void step(real h) {
    real h_step = 0.0;

    real h_step_max = h;

    std::vector<vec3> &f = asawa::get_vec_data(*_M, _i);
    std::vector<vec3> &x_s = asawa::get_vec_data(*_M, 0);

    std::array<real, 2> meanstd = va::mean_std<vec3, real>(f);
    real fnorm = meanstd[0];
    real fstddev = meanstd[1];
    std::cout << "f mean/std: " << fnorm << " " << fstddev << std::endl;

    if (fnorm < 1.0e-3 && fstddev < 1.0e-3)
      return;

    while (h_step <= h_step_max) {
      std::cout << "h_step: " << h_step << std::endl;

      std::cout << "calc_dt..." << std::endl;
      // clip_outliers(f, 2.0);
      real hi = this->calc_dt(f, h);
      std::cout << "hi: " << hi << std::endl;
      hi = std::min(hi, h);
      hi = std::max(hi, kmin * h);
      // hi = std::min(h_step_max - h_step, hi);
      std::cout << "hi min: " << hi << std::endl;
      for (int i = 0; i < f.size(); i++) {
        x_s[i] += hi * f[i];
      }

      h_step += hi;
      std::cout << "before step..." << std::endl;
      for (int i = 0; i < f.size(); i++) {
        if (f[i].hasNaN())
          std::cout << f[i].transpose() << " ";
      }
      for (int i = 0; i < 2; i++)
        __dynamic->step();

      std::cout << "after step..." << std::endl;
      for (int i = 0; i < f.size(); i++) {
        if (f[i].hasNaN())
          std::cout << f[i].transpose() << " ";
      }
    }
  }
  void set_kmin(real kmin) { this->kmin = kmin; }

  real kmin = 1e-1;
  asawa::shell::dynamic::ptr __dynamic;
  index_t _i;
  real _eps, _h;
};

} // namespace duchamp
} // namespace gaudi

#endif