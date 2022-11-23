#ifndef __GROWTH_EXPERIMENT__
#define __GROWTH_EXPERIMENT__

#include "experiment.hpp"

namespace duchamp {

template <typename SPACE>
class action_policy
    : public asawa::vertex_policy_t<SPACE, typename SPACE::coordinate_type> {
public:
  M2_TYPEDEFS;
  action_policy(typename SPACE::vertex_index id)
      : vertex_policy_t<SPACE, coordinate_type>(id) {
    std::cout << " id: " << int(this->_id) << std::endl;
  }
  virtual void calc(int i, edge_ptr &e, asawa::op_type op) {

    face_vertex_ptr v0 = e->v1();
    face_vertex_ptr v1 = v0->prev();
    face_vertex_ptr v2 = e->v2();
    face_vertex_ptr v3 = v2->prev();

    coordinate_type t0 = v0->vertex()->template get<coordinate_type>(this->_id);
    coordinate_type t1 = v1->vertex()->template get<coordinate_type>(this->_id);
    coordinate_type t2 = v2->vertex()->template get<coordinate_type>(this->_id);
    coordinate_type t3 = v3->vertex()->template get<coordinate_type>(this->_id);

    if (op == asawa::op_type::split) {
      this->_vals[i] =
          coordinate_type::Zero(); // 0.0 * (t0 + t2) - 0.0 * (t1 + t3);
      //_vals[i] = 0.0;
    } else {
      coordinate_type mx = t0.norm() > t2.norm() ? t0 : t2;
      // this->_vals[i] = t0; // might be zero might not
      this->_vals[i] = mx; // might be zero might not
    }

    //_vals[i] = -0.25 * (t0 + t2);
    return;
  }
};

template <typename T> class growth_space {
public:
  typedef T double_type;
  typedef T real;
  typedef std::complex<T> complex;
  // always use homogeneous coordinates, provides decent error checking
  typedef Eigen::Matrix<T, 3, 1> coordinate_type;

  typedef line<T, coordinate_type> line_type;
  typedef triangle<T, coordinate_type> triangle_type;
  typedef swept_point<T, coordinate_type> swept_point_type;
  typedef swept_triangle<T, coordinate_type> swept_triangle_type;
  typedef bounding_box<T, coordinate_type> box_type;
  typedef Eigen::Matrix<T, 2, 2> mat2;
  typedef Eigen::Matrix<T, 3, 3> mat3;
  typedef Eigen::Matrix<T, 4, 4> mat4;
  typedef Eigen::Matrix<T, 4, 3> mat43;

  typedef unsigned short ushort;
  typedef unsigned int uint;
  typedef unsigned long ulong;

  typedef double double_t;
  typedef float float_t;

  typedef Eigen::Quaternion<T> quat;
  typedef Eigen::Matrix<T, 2, 1> vec2;
  typedef Eigen::Matrix<T, 3, 1> vec3;
  typedef Eigen::Matrix<T, 4, 1> vec4;

  typedef Eigen::Matrix<uint, 2, 1> uint2;
  typedef Eigen::Matrix<uint, 2, 1> uint4;

  typedef Eigen::Matrix<T, 2, 1> int2;
  typedef Eigen::Matrix<T, 4, 1> int3;
  typedef Eigen::Matrix<T, 4, 1> int4;

  enum class face_vertex_index { BARY = 0, MAXINDEX = 1 };

  enum class edge_index { PHI0 = 0, MAXINDEX = 1 };

  enum class vertex_index {
    COORDINATE = 0,
    COLOR = 1,
    SMOOTH = 2,
    ACT = 3,

    MAXINDEX = 4
  };

  enum class face_index { NORMAL = 0, CENTER = 1, AREA = 2, MAXINDEX = 3 };

  static storage_type get_type(face_vertex_index idx) {
    switch (idx) {
    case face_vertex_index::BARY:
      return storage_type::REAL;
    default:
      return storage_type::SIZE;
    };
  }
  static storage_type get_type(edge_index idx) {
    switch (idx) {
    case edge_index::PHI0:
      return storage_type::REAL;
    default:
      return storage_type::SIZE;
    }
  }

  static storage_type get_type(face_index idx) {
    switch (idx) {
    case face_index::NORMAL:
      return storage_type::VEC3;
    case face_index::CENTER:
      return storage_type::VEC3;
    case face_index::AREA:
      return storage_type::REAL;
    default:
      return storage_type::SIZE;
    }
  }

  static storage_type get_type(vertex_index idx) {
    switch (idx) {
    case vertex_index::COORDINATE:
      return storage_type::VEC3;
    case vertex_index::COLOR:
      return storage_type::VEC3;
    case vertex_index::SMOOTH:
      return storage_type::REAL;
    case vertex_index::ACT:
      return storage_type::VEC3;
    default:
      return storage_type::SIZE;
    }
  }
};

using namespace asawa;
template <typename SPACE> class growth_target_base {
public:
  M2_TYPEDEFS;

  typedef std::shared_ptr<growth_target_base<SPACE>> ptr;
  static ptr create() { return std::make_shared<growth_target_base<SPACE>>(); }

  virtual void init(const real &scale) {}

  virtual void pre() {}

  virtual size_t size() { return 0; }

  virtual real weight(const coordinate_type &X, const coordinate_type &N) {
    return 0.0;
  }

  virtual void accumWeights(const real &area, const coordinate_type &x,
                            const coordinate_type &N) {}

  virtual coordinate_type gradient(const real &area, const coordinate_type &x,
                                   const coordinate_type &N, const real &C) {
    return coordinate_type::Zero();
  }

  virtual void render() {}
};

template <typename SPACE>
class point_target : public growth_target_base<SPACE> {
public:
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;

  typedef std::shared_ptr<point_target<SPACE>> ptr;
  static ptr create() { return std::make_shared<point_target<SPACE>>(); }

  virtual size_t size() { return _targets.size(); }
  virtual void init(const real &scale) {
    std::mt19937_64 rng;
    std::uniform_real_distribution<real> unif(-1.0, 1.0);

    for (int i = 0; i < 128; i++) {
      double c0 = 2.0 * unif(rng);
      double c1 = 2.0 * unif(rng);
      double c2 = 2.0 * unif(rng);
      coordinate_type c(c0, c1, c2);
      coordinate_type trans(2.0, 0.0, 0.0);

      _targets.push_back(c + trans);
    }
    _scale = scale;
    _tdist = this->avgTargDist();
  }

  real avgTargDist() {
    real accum = 0.0;
    real N = 0.0;
    for (int i = 0; i < _targets.size(); i++) {
      for (int j = 0; j < _targets.size(); j++) {
        if (i == j)
          continue;
        coordinate_type ti = _targets[i];
        coordinate_type tj = _targets[j];
        real dX = (tj - ti).norm();
        accum += dX;
        N += 1.0;
      }
    }
    return (accum / N);
  }

  real lweight(coordinate_type X, coordinate_type N) {
    double dX = X.norm();
    double eps = 0.1 * _scale;
    double pdx = pow(dX, 0.5);

    double pdx2 = pow(dX, 4.0);
    double NdX = N.dot(X / dX);

    return pow(NdX, 32) / (pdx + eps);
  }

  virtual real weight(const coordinate_type &X, const coordinate_type &N) {
    std::mt19937_64 rng;
    std::normal_distribution<> norm(8.0, 1.0);

    double dX = X.norm();
    double eps = 1e-6 * _scale;
    double pdx = pow(dX, norm(rng));
    return 1.0 / (pdx + eps);
  }

  virtual void pre() { _weights = std::vector<real>(_targets.size(), 0.0); }

  virtual void accumWeights(const real &area, const coordinate_type &x,
                            const coordinate_type &N) {

    std::mt19937_64 rng;
    std::uniform_real_distribution<real> unif01(0.0, 1.0);
    _totalWeight = 0.0;
    for (int j = 0; j < _targets.size(); j++) {
      coordinate_type X = _targets[j] - x;
      real w = area * weight(X, N);
      _weights[j] += w;
      _totalWeight += w;
    }
  }

  coordinate_type localGradient(const int &j, const real &area,
                                const coordinate_type &x,
                                const coordinate_type &N, const real &tdist) {

    coordinate_type g = coordinate_type::Zero();
    coordinate_type X = _targets[j] - x;
    double dX = X.norm();
    if (dX < 0.25 * tdist && N.dot(X) > 0) {
      coordinate_type u = N.cross(X.normalized());
      coordinate_type Nu = Eigen::AngleAxisd(M_PI / 12.0, u) * N;
      coordinate_type gi = area * lweight(X, N) * Nu;
      g += gi;
    }

    return g;
  }

  virtual coordinate_type gradient(const real &area, const coordinate_type &x,
                                   const coordinate_type &N, const real &C) {

    coordinate_type g = coordinate_type::Zero();
    for (int j = 0; j < _targets.size(); j++) {

      coordinate_type X = _targets[j] - x;

      double dX = X.norm();
      coordinate_type Xn = X / dX;

      double p = 3.0;
      double k = 2.0;

      double lc = _scale;
      double eps = _scale;

      double el = (_weights[j] + 1e0 * eps);
      double w = area * weight(X, N) * (1e-2 / sqrt(el) + 1.0 / el);

      double NdX = N.dot(X);

      double Xlc = dX - lc;
      coordinate_type gXlc = 2.0 * Xlc * Xn;
      coordinate_type gNdX = (-NdX * Xn + N) / (dX + _scale);

      coordinate_type gNdXp =
          p * pow(NdX, p - 1) * gNdX; // alighnment with point gradient

      // coordinate_type NxX = X.cross(N);    // perp to grad
      coordinate_type gNxX = -Xn.cross(N); // perp to grad
      coordinate_type gNxXdX = gNxX;
      g += 1.0 * w * X;

      // g += 1.0 * w * (1.0 * X + 0.01 * gXlc + 0.01 * gNdX);
      g += C * localGradient(j, area, x, N, _tdist);
    }
    return g;
  }

  virtual void render() {
    for (int j = 0; j < _targets.size(); j++) {
      auto p0 = _targets[j];
      real sc = 0.5 * _scale;
      gg::geometry_logger::box(p0, vec3(sc, sc, sc), vec4(0.5, 0.8, 0.6, 1.0));
    }
    std::cout << "rendering debug" << std::endl;

    real tdist = _tdist;

    int i = 0;
    for (auto &w : _weights) {
      int ip = (i + 1) % _weights.size();
      coordinate_type t0 = _targets[i];
      coordinate_type t1 = _targets[ip];

      real sc = 1e-2 * log(1.0 / w);
      vec4 color = sc < 0 ? vec4(1.0, 0.0, 0.0, 1.0) : vec4(0.0, 1.0, 0.0, 1.0);
      if (sc > 0)
        gg::geometry_logger::box(_targets[i], vec3(sc, sc, sc), color);
      if ((t1 - t0).norm() < tdist)
        gg::geometry_logger::line(t0, t1, vec4(0.5, 0.5, 0.5, 1.0));
      i++;
    }
  }

  real _scale = 0.1;
  real _tdist = 0;
  real _totalWeight = 0.0;
  std::vector<real> _weights;
  std::vector<coordinate_type> _targets;
};

#if 0
template <typename SPACE>
class phase_point_target : public growth_target_base<SPACE> {

  virtual void init() {
    std::mt19937_64 rng;
    std::uniform_real_distribution<real> unif(-1.0, 1.0);

    for (int i = 0; i < 128; i++) {
      double c0 = 2.0 * unif(rng);
      double c1 = 2.0 * unif(rng);
      double c2 = 2.0 * unif(rng);
      coordinate_type c(c0, c1, c2);
      coordinate_type trans(2.0, 0.0, 0.0);

      _targets.push_back(c + trans);
    }
    _phase = std::vector<real>(_targets.size(), 1.0);
    for (auto &p : _phase) {
      p = unif(rng);
    }
  }

virtual real
weight(const coordinate_type &X, const coordinate_type &N, const real &ph) {
  std::mt19937_64 rng;
  std::normal_distribution<> norm(8.0, 2.0);

  double dX = X.norm();
  double eps = 0.001 * _scale;
  // double pdx = pow(dX, 6.0);
  double pdx = pow(dX, norm(rng));

  real cs = cos(dX * M_PI / (100.0 * _scale) + 0.005 * M_PI * this->t + ph);
  real ss = sin(dX * M_PI / (80.0 * _scale) + 0.006 * M_PI * this->t + ph);

  return (cs * cs + 0.1 * ss * ss) / (pdx + eps);
}

std::vector<coordinate_type> _targets;
std::vector<real> _phase;
};
#endif

template <typename SPACE> class mean_shift {
public:
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;

  typedef std::shared_ptr<mean_shift<SPACE>> ptr;
  static ptr create(real s) { return std::make_shared<mean_shift<SPACE>>(s); }

  mean_shift(real scale) : _scale(scale) { init_targets(); }

  void init_velocity(surf_ptr surf) {
    auto &verts = surf->get_vertices();
    coordinate_array act(verts.size(), coordinate_type::Zero());
    ci::set<SPACE, coordinate_type>(surf, act, SPACE::vertex_index::ACT);
  }

  void init_targets() {
    _targets = point_target<SPACE>::create();
    _targets->init(_scale);
  }

  coordinate_array step(surf_ptr surf, const coordinate_array &p0) {

    std::mt19937_64 rng;
    std::uniform_real_distribution<real> unif01(0.0, 1.0);
    std::uniform_real_distribution<real> unif11(-1.0, 1.0);

    vertex_array &verts = surf->get_vertices();
    coordinate_array normals = asawa::ci::get_vertex_normals<SPACE>(surf);

    int i = 0;
    coordinate_array p1(p0);
    coordinate_array vg(p0.size(), coordinate_type::Zero());

    // coordinate_array fN = calcPotential(surf, act, p0, 0.5 * _scale);
    _targets->pre();
    i = 0;
    for (auto &v : verts) {

      coordinate_type N = normals[i];
      coordinate_type x = p0[i];
      real area = ci::area<SPACE>(v);
      _targets->accumWeights(area, x, N);
      i++;
    }

    i = 0;
    calder::mesh_calculator<SPACE> calc;
    std::vector<typename SPACE::mat43> cov =
        calc.template covariance(surf, p0, 2.0 * _scale);

    for (auto &v : verts) {
      typename SPACE::mat3 U = cov[i].block(0, 0, 3, 3);
      coordinate_type s = cov[i].row(3);

      real s12 = s[1] / s[2] - 1.0;
      real s0 = 1e-2 * (s[0] / s[1] / s12);

      coordinate_type N = normals[i];
      coordinate_type x = p0[i];
      real area = ci::area<SPACE>(v);
      vg[i] += _targets->gradient(area, x, N, s0);
      i++;
    }

    coordinate_array v1 = calder::calcAvgVerts(surf, vg, p0, 2.0 * _scale);

    for (int i = 0; i < p0.size(); i++) {
      // gg::geometry_logger::line(p0[i], p0[i] + 10.0 * v1[i],
      //                           vec4(1.0, 0.0, 0.5, 1.0));
      p1[i] += 1.0 * v1[i];
    }

    render_targets();

    this->t += 1.0;

    return p1;
  }

  void render_targets() { _targets->render(); }

  typename growth_target_base<SPACE>::ptr _targets;
  real _scale = 0.1;
  real t = 0.0;
};

template <typename SPACE> class worms {
public:
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;

  typedef std::shared_ptr<worms<SPACE>> ptr;
  static ptr create(real s) { return std::make_shared<worms<SPACE>>(s); }

  worms(real scale) : _scale(scale) { init_targets(); }

  void init_velocity(surf_ptr surf) {
    auto &verts = surf->get_vertices();
    coordinate_array act(verts.size(), coordinate_type::Zero());
    ci::set<SPACE, coordinate_type>(surf, act, SPACE::vertex_index::ACT);
  }

  void init_targets() {
    _targets = point_target<SPACE>::create();
    _targets->init(_scale);
  }

  void calc_stdev(const std::vector<real> &K, real &mean, real &stdev) {

    real sum = std::accumulate(std::begin(K), std::end(K), 0.0);
    mean = sum / K.size();
    real accum = 0.0;
    std::for_each(std::begin(K), std::end(K),
                  [&](const double d) { accum += (d - mean) * (d - mean); });
    stdev = sqrt(accum / (K.size() - 1));

    std::cout << " K mean: " << mean << " stdev: " << stdev << std::endl;
  }

  void frame(vertex_ptr v, coordinate_type &s, mat3 &U) {

    U = mat3::Zero();
    coordinate_type p0 = ci::get_coordinate<SPACE>(v);
    for_each_vertex<SPACE>(v, [&U, p0](face_vertex_ptr fv) {
      coordinate_type p1 = ci::get_coordinate<SPACE>(fv->next()->vertex());
      coordinate_type N = ci::normal<SPACE>(fv->next()->face());
      real A = ci::area<SPACE>(fv->next()->face());

      coordinate_type dp = A * N;
      U += dp * dp.transpose();
    });

    Eigen::JacobiSVD<mat3> svd(U, Eigen::ComputeFullU);
    U = svd.matrixU();
    s = svd.singularValues();
  }

  real angle(const coordinate_type &a, const coordinate_type &b) {
    coordinate_type an = a.normalized();
    coordinate_type bn = b.normalized();
    coordinate_type x = an.cross(bn);
    return atan2(x.norm(), an.dot(bn));
  }

  bool inCone(const coordinate_type &p, const coordinate_type &heading,
              real thet, real r, const coordinate_type &t) {
    auto H = heading.normalized();
    bool inCone = false;

    coordinate_type dp = t - p;
    double d = dp.norm();
    dp /= d;
    real ti = angle(dp, H);
    if (ti < thet && d < r) {
      inCone = true;
    }

    return inCone;
    typedef Eigen::Matrix<real, 4, 1> vec4;
  }

  bool inCone(const coordinate_type &p, const coordinate_type &heading,
              real thet, real r, const std::vector<coordinate_type> &targets,
              coordinate_type &g) {
    auto H = heading.normalized();
    real dmin = 99999;
    bool inCone = false;
    real tmin = 99999;
    for (int j = 0; j < targets.size(); j++) {
      coordinate_type dp = targets[j] - p;
      double d = dp.norm();
      dp /= d;
      real ti = angle(dp, H);
      if (ti < thet && ti < tmin && d < r) {
        inCone = true;
        tmin = ti;
        g = dp;
      }
    }
    return inCone;
  }

  coordinate_array step(surf_ptr surf, const coordinate_array &p0) {

    std::mt19937_64 rng;
    std::uniform_real_distribution<real> unif01(0.0, 1.0);

    vertex_array &verts = surf->get_vertices();
    coordinate_array normals = asawa::ci::get_vertex_normals<SPACE>(surf);
    coordinate_array act =
        ci::get<SPACE, coordinate_type>(surf, SPACE::vertex_index::ACT);

    //////////////////////////
    // heat dist
    //////////////////////////

    std::vector<real> f(act.size());
    int i = 0;
    for (auto a : act) {
      real sg = va::sgn(act[i].dot(normals[i]));
      if (sg > 0)
        f[i] = act[i].norm();
      else
        f[i] = -act[i].norm();
      i++;
    }
    double dt = 0.5 * pow(_scale, 2);
    asawa::laplacian<SPACE, real> lap(surf);
    std::vector<real> dist = lap.heatDist(f, dt);
    //////////////////////////
    // curvature
    //////////////////////////

    asawa::cotan_curvature<SPACE> curve(surf);
    std::vector<real> K = curve();
    // std::for_each(K.begin(), K.end(), [](auto &k) { k = 1.0 / k; });
    real K_mean, K_stdev;
    calc_stdev(K, K_mean, K_stdev);

    //////////////////////////
    // probability
    //////////////////////////

    real N = double(verts.size());
    real d0 = 0.0;
    real dp = 0.0;
    real dn = 0.0;

    i = 0;
    for (auto &v : verts) {
      coordinate_type av = act[i];
      real a = av.norm();
      d0 += fabs(a) / N;
      i++;
    }

    std::cout << "current density: " << d0 << std::endl;
    double max_d = 0.0015;
    double prob_d = 0.25 / N / max_d;
    double d1 = d0;
    std::cout << "max density: " << max_d << ", prob density: " << prob_d
              << std::endl;

    i = 0;
    double theta = M_PI / 12.0;
    for (auto &v : verts) {
      coordinate_type av = act[i];
      real a = av.norm();

      if (a < 1e-1)
        a = 0.0;

      if (unif01(rng) < prob_d && d1 < max_d && a < 1e-4) {
        // if (K[i] > mean + 0.0 * stdev) {
        coordinate_type g;
        bool sees_target = inCone(p0[i], normals[i], 6.0 * theta,
                                  1000.0 * _scale, _targets, g);
        sees_target = sees_target || (!sees_target && unif01(rng) > 0.65);
        if (unif01(rng) < 0.75 && dist[i] < 0.30 && sees_target) {
          av = normals[i];
          a = 1.0;
        } else if (dist[i] > 0.70 && K[i] < K_mean - 0.5 * K_stdev) {
          a = 1.0;
          av = -normals[i];
        }
      }

      mat3 F;
      coordinate_type s;
      frame(v, s, F);

      if (fabs(a) > 0) {
        // if (unif01(rng) > 0.99)
        //   a = 0.0;
        // double theta = M_PI / 6.0 * unif01(rng);

        real sg = va::sgn(av.dot(normals[i]));
        coordinate_type g;
        bool sees_target =
            inCone(p0[i], av, 2.0 * theta, 10.0 * _scale, _targets, g);

        if (!sees_target) {

          a *= (1.0 - 0.04 * unif01(rng)); // cool a down

          coordinate_type uu = av.normalized();
          sees_target =
              inCone(p0[i], av, 4.0 * theta, 1000.0 * _scale, _targets, g);

          coordinate_type vv = uu.cross(normals[i]);
          coordinate_type ww = uu.cross(vv);

          uu.normalize();
          vv.normalize();
          ww.normalize();

          //   coordinate_type v = F.block(0, 2, 3, 1);
          if (sees_target && normals[i].dot(uu) > 0.0) {
            // if (normals[i].dot(uu) < 0.0) { theta *= -1;}
            vv = uu.cross(g);
            av = Eigen::AngleAxisd(theta, vv) * av;
          } else {

            vv = F.block(0, 1, 3, 1);
            ww = F.block(0, 2, 3, 1);

            if (unif01(rng) > 0.5 && !sees_target)
              theta *= -1;

            if (unif01(rng) > 0.5)
              av = Eigen::AngleAxisd(theta, vv) * av;
            else {
              av = Eigen::AngleAxisd(theta, ww) * av;
            }
          }
        }
      }

      d1 += a / N;
      act[i] = a * av.normalized();
      i++;
    }
    ci::set<SPACE, coordinate_type>(surf, act, SPACE::vertex_index::ACT);
    render_targets();
    coordinate_array p1(p0);
    p1 = update_position(surf, p0);
    p1 = decay(surf, p1, dist);
    return p1;
  }

  coordinate_array update_position(surf_ptr surf, const coordinate_array &p0) {
    std::mt19937_64 rng;
    std::uniform_real_distribution<real> unif(-1.0, 1.0);

    vertex_array &verts = surf->get_vertices();
    coordinate_array normals = asawa::ci::get_vertex_normals<SPACE>(surf);
    coordinate_array act =
        ci::get<SPACE, coordinate_type>(surf, SPACE::vertex_index::ACT);
    int i = 0;

#if 0
    //coordinate_array fN = calcPotential(surf, act, p0, 2.0 * _scale);
    i = 0;
    for (auto f : fN) {
      //      gg::geometry_logger::line(p0[i], p0[i] + 1.0 * f,
      //                                vec4(1.0, 1.0, 0.0, 1.0));
      act[i] += 1.0e-5 * f;
      i++;
    }
#endif
    coordinate_array p1(p0);

    double C_N = 1.0 * _scale;
    i = 0;
    for (auto v : verts) {
      // auto f = va::sgn(wN[i].dot(normals[i])) * wN[i].norm();
      coordinate_type vi = act[i];
      real sg = va::sgn(vi.dot(normals[i]));

      coordinate_type Ni = sg * normals[i];

      real vm = vi.norm();
      coordinate_type B = Ni.cross(vi.normalized());
      // double thet = atan2(vv.norm(), uu.dot(Ni));

      real thet = angle(vi, Ni);
      Eigen::AngleAxisd R = Eigen::AngleAxisd(thet, B);

      for_each_vertex<SPACE>(
          v, [i, &p0, &p1, &vm, sg, R, &unif, &rng, &normals,
              C_N](typename asawa::surf<SPACE>::face_vertex_ptr fv) {
            int j = fv->next()->vertex()->position_in_set();
            coordinate_type Nj = sg * normals[j];
            Nj = R * Nj;
            // p1[j] += 0.5 * C_N * unif(rng) * p1[j];
            p1[j] += 16.0 * C_N * vm * Nj;
          });

      // p1[i] += 0.5 * C_N * unif(rng) * p1[i];
      p1[i] -= 4.0 * C_N * act[i];

      i++;
    }

    return p1;
  }

  coordinate_array decay(surf_ptr surf, const coordinate_array &p0,
                         const vector<real> &dist) {

    coordinate_array p1(p0);

    coordinate_array normals = asawa::ci::get_vertex_normals<SPACE>(surf);

    int i = 0;
#if 1
    for (auto &p : p1) {
      real alpha = 1.0;
      // real d = cos(M_PI * alpha * dist[i]); // 2.0 * (0.5 - dist[i]);
      real d = 2.0 * (0.5 - dist[i]);
      coordinate_type dN = pow(d, 5.0) * normals[i];
      // gg::geometry_logger::line(p0[i], p0[i] + 0.1 * dN,
      //                           vec4(0.5, 0.0, 1.0, 1.0));
      //  if (d < 0) {
      p += 2.0e-2 * dN;
      //}
      i++;
    }
#endif

    calder::mesh_calculator<SPACE> calc;
#if 0
    std::vector<typename SPACE::mat43> covar =
        calc.template covariance(surf, p0, 2.0 * _scale);
    i = 0;
    for (auto &US : covar) {
      typename SPACE::mat3 U = US.block(0, 0, 3, 3);
      coordinate_type s = US.row(3);
      double s0 = s[0] / (s[0] + s[1] + s[2]);
      coordinate_type uc0 = U.col(0).transpose();
      coordinate_type N = normals[i];
      coordinate_type N0 = s0 * va::dot(uc0, N) * uc0;
      if (dist[i] < 0)
        p1[i] += 2e-2 * dist[i] * N0;
    }
#endif
#if 0
    vector<real> fAreas = asawa::ci::get_face_areas<SPACE>(surf);
    vector<real> vAreas = asawa::ci::get_vertex_areas<SPACE>(surf);
    vector<coordinate_type> G =
        calc.template gravitation<coordinate_type, real>(surf, fAreas, p0,
                                                         vAreas, 2.0 * _scale);
    i = 0;
    for (auto &g : G) {
      coordinate_type N = normals[i];
      coordinate_type N0 = va::dot(g, N) * g.normalized();
      // std::cout << g.transpose() << std::endl;
      real d = cos(M_PI * dist[i]); // 2.0 * (0.5 - dist[i]);

      if (d < 0) {
        // gg::geometry_logger::line(p0[i], p0[i] + 1.0e5 * d * N0,
        //                           vec4(0.5, 0.0, 1.0, 1.0));
        p1[i] -= 2e4 * d * N0;
      }
      i++;
    }
#endif

    return p1;
  }

  void render_targets() { _targets->render(); }
  real _scale = 0.1;
  typename growth_target_base<SPACE>::ptr _targets;
};

template <typename SPACE>
class mean_shift_experiment : public bend_experiment<SPACE> {
public:
  M2_TYPEDEFS;

  typedef std::shared_ptr<mean_shift_experiment<SPACE>> ptr;
  static ptr create() {
    return std::make_shared<mean_shift_experiment<SPACE>>();
  }

  virtual void _init() {
    asawa::affine<SPACE> aff;
    aff.move_to(this->_surf, typename SPACE::coordinate_type(-2.0, 0.0, 0.0));

    this->bend_experiment<SPACE>::_init();

    this->_integrator->add_vertex_policy(
        new action_policy<SPACE>(SPACE::vertex_index::ACT));
  };

  virtual void _step(const int &frame) {

    if (!_grower) {
      _grower = duchamp::mean_shift<SPACE>::create(this->_integrator->_max);
      _grower->init_velocity(this->_surf);
    }
    std::vector<typename SPACE::vec3> p0 =
        asawa::ci::get_coordinates<SPACE>(this->_surf);
    std::vector<typename SPACE::vec3> p1 = _grower->step(this->_surf, p0);

    this->willmore(p0);
    this->bending(p0, p1);
  };
  typename duchamp::mean_shift<SPACE>::ptr _grower;
};

} // namespace duchamp
#endif