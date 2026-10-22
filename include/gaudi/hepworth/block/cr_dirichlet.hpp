#ifndef GAUDI_HEPWORTH_CR_DIRICHLET_HPP
#define GAUDI_HEPWORTH_CR_DIRICHLET_HPP

#include <algorithm>
#include <array>
#include <limits>
#include <utility>
#include <vector>

#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/hepworth/block/cr_shape.hpp"
#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/hepworth/blocks/shell_cr_block.hpp"
#include "gaudi/kusama/vector_dirichlet.hpp"
#include "gaudi/kusama/vector_dirichlet_guided.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// One incident face, as vertex ids, plus the vertex triples of its 1-ring.
struct dirichlet_face {
  std::array<int, 3> v{-1, -1, -1};
  std::vector<std::array<int, 3>> stencil;
};

template <class At>
inline vec3 face_verts_cross(At at, const std::array<int, 3> &v) {
  if (v[0] < 0 || v[1] < 0 || v[2] < 0)
    return vec3::Zero();
  const vec3 x0 = at(v[0]);
  const vec3 x1 = at(v[1]);
  const vec3 x2 = at(v[2]);
  return (x1 - x0).cross(x2 - x0);
}

template <class At>
inline vec3 face_verts_center(At at, const std::array<int, 3> &v) {
  return (at(v[0]) + at(v[1]) + at(v[2])) / real(3);
}

template <class At>
inline bool face_tmin_of(At at, const std::array<int, 3> &v,
                         const std::vector<std::array<int, 3>> &stencil, vec3 &dir) {
  const vec3 cr = face_verts_cross(at, v);
  if (cr.squaredNorm() < real(1e-24))
    return false;
  const vec3 n = cr.normalized();
  const vec3 a = std::abs(n[0]) > real(0.5) ? vec3(0, 1, 0) : vec3(1, 0, 0);
  vec3 u = n.cross(a);
  if (u.squaredNorm() < real(1e-24))
    u = vec3(0, 0, 1).cross(n);
  if (u.squaredNorm() < real(1e-24))
    return false;
  u.normalize();
  const vec3 vv = n.cross(u);
  if (!(vv.squaredNorm() > real(1e-24)))
    return false;
  const vec3 foot = face_verts_center(at, v);
  const int m = static_cast<int>(stencil.size());
  if (m < 3)
    return false;
  Eigen::MatrixXd A(m, 3);
  Eigen::VectorXd b(m);
  for (int i = 0; i < m; ++i) {
    const vec3 d = face_verts_center(at, stencil[static_cast<size_t>(i)]) - foot;
    const real xi = d.dot(u);
    const real yi = d.dot(vv);
    A(i, 0) = xi * xi;
    A(i, 1) = xi * yi;
    A(i, 2) = yi * yi;
    b(i) = d.dot(n);
  }
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  const Eigen::Vector3d coef = svd.solve(b);
  if (!coef.allFinite())
    return false;
  Eigen::Matrix2d H;
  H(0, 0) = 2.0 * coef(0);
  H(0, 1) = coef(1);
  H(1, 0) = coef(1);
  H(1, 1) = 2.0 * coef(2);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(H);
  if (es.info() != Eigen::Success)
    return false;
  const Eigen::Vector2d ev = es.eigenvectors().col(0);
  dir = ev[0] * u + ev[1] * vv;
  if (!dir.allFinite() || dir.squaredNorm() < real(1e-24))
    return false;
  dir.normalize();
  return true;
}

/// Incident-face guidance of one edge, before the mesh-wide sign walk.
template <class At>
inline void unoriented_edge_guidance(At at, int v0, int v1,
                                    const std::vector<dirichlet_face> &faces,
                                    real &gp, real &gq) {
  gp = real(0);
  gq = real(0);
  std::vector<vec3> tmin(faces.size(), vec3::Zero());
  std::vector<char> have(faces.size(), 0);
  for (int f = 0; f < static_cast<int>(faces.size()); ++f) {
    vec3 d;
    if (face_tmin_of(at, faces[static_cast<size_t>(f)].v,
                     faces[static_cast<size_t>(f)].stencil, d)) {
      tmin[static_cast<size_t>(f)] = d;
      have[static_cast<size_t>(f)] = 1;
    }
  }
  vec3 acc = vec3::Zero();
  for (int f = 0; f < static_cast<int>(faces.size()); ++f) {
    const vec3 cr = face_verts_cross(at, faces[static_cast<size_t>(f)].v);
    const real len = cr.norm();
    if (!(len > real(1e-20)))
      continue;
    acc += real(0.5) * len * (cr / len);
  }
  const real s = acc.norm();
  vec3 n_avg = vec3::UnitZ();
  if (s > real(1e-20))
    n_avg = acc / s;
  vec3 e_par = vec3::UnitX();
  if (v0 >= 0 && v1 >= 0) {
    vec3 t = at(v1) - at(v0);
    const real len = t.norm();
    if (len > real(1e-20)) {
      const real sign = v0 <= v1 ? real(1) : real(-1);
      e_par = sign * t / len;
    }
  }
  const vec3 e_perp = cr_edge_perp(e_par, n_avg);
  real wsum = 0;
  real acc_p = 0;
  real acc_q = 0;
  for (int f = 0; f < static_cast<int>(faces.size()); ++f) {
    if (!have[static_cast<size_t>(f)])
      continue;
    vec3 d = tmin[static_cast<size_t>(f)] -
             n_avg * tmin[static_cast<size_t>(f)].dot(n_avg);
    const real dn = d.norm();
    if (!(dn > real(1e-14)))
      continue;
    d /= dn;
    const vec3 cr = face_verts_cross(at, faces[static_cast<size_t>(f)].v);
    const real w = real(0.5) * cr.norm();
    wsum += w;
    acc_p += w * d.dot(e_par);
    acc_q += w * d.dot(e_perp);
  }
  if (wsum > real(1e-20)) {
    gp = acc_p / wsum;
    gq = acc_q / wsum;
  }
}

/// Stein vector Dirichlet rows. One row per edge. The shell is walked while
/// the row is built. `project` reads positions from block 0 and writes this
/// edge's piece of λ M g. The coefficient unknown is block 1.
/// g is the incident-face guidance. The mesh-wide sign walk is not a row.
class cr_dirichlet : public block_constraint {
public:
  DEFINE_CREATE_FUNC(cr_dirichlet)

  cr_dirichlet(int dof, int v0, int v1, int nE, real mass_par, real mass_perp,
               std::vector<std::pair<int, real>> Lpar,
               std::vector<std::pair<int, real>> Lperp,
               std::vector<dirichlet_face> faces, const real &w, real lambda,
               std::vector<sim_block::ptr> blocks)
      : block_constraint({dof}, w, std::move(blocks)), _v0(v0), _v1(v1), _nE(nE),
        _mass_par(mass_par), _mass_perp(mass_perp), _Lpar(std::move(Lpar)),
        _Lperp(std::move(Lperp)), _faces(std::move(faces)), _lambda(lambda) {
    _rest_initialized = true;
  }

  std::string name() override { return "cr_dirichlet"; }

  int dof() const {
    return _ids.empty() ? -1 : static_cast<int>(_ids[0]);
  }

  void project(const vecX &q, vecX &p) override {
    if (_blocks.empty() || !(_w > real(0)) || !(_lambda > real(0)))
      return;
    auto at = [&](int i) { return _blocks[0]->get_vec3(i, q); };
    real gp = 0;
    real gq = 0;
    unoriented_edge_guidance(at, _v0, _v1, _faces, gp, gq);
    if (!std::isfinite(gp))
      gp = real(0);
    if (!std::isfinite(gq))
      gq = real(0);
    p(_id0) = _w * _lambda * _mass_par * gp;
    p(_id0 + 1) = _w * _lambda * _mass_perp * gq;
  }

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    if (_blocks.size() < 2 || !(_w > real(0)))
      return;
    const index_t off = _blocks[1]->get_offset_idx(0);
    for (const auto &t : _Lpar)
      triplets.push_back(trip(id0, off + t.first, _w * t.second));
    for (const auto &t : _Lperp)
      triplets.push_back(trip(id0 + 1, off + t.first, _w * t.second));
    if (_lambda > real(0)) {
      const int e = dof();
      if (e >= 0) {
        triplets.push_back(
            trip(id0, off + e, _w * _lambda * _mass_par));
        triplets.push_back(
            trip(id0 + 1, off + e + _nE, _w * _lambda * _mass_perp));
      }
    }
    id0 += 2;
  }

private:
  int _v0 = -1;
  int _v1 = -1;
  int _nE = 0;
  real _mass_par = 0;
  real _mass_perp = 0;
  std::vector<std::pair<int, real>> _Lpar;
  std::vector<std::pair<int, real>> _Lperp;
  std::vector<dirichlet_face> _faces;
  real _lambda = 0;
};


struct cr_complex {
  std::vector<cr_edge> edges;
  std::vector<cr_tri> tris;
  std::vector<int> slot_map;
  std::vector<int> face_local;
};

inline bool build_cr_complex(asawa::shell::shell &shell, cr_complex &out) {
  out = cr_complex{};
  if (!shell.verts_are_dense_packed())
    return false;
  const int nE = kusama::build_compact_edge_dof_map(shell, out.slot_map);
  if (nE <= 0)
    return false;
  std::vector<asawa::shell::CornerId> dof_to_corner;
  kusama::build_dof_to_corner(shell, out.slot_map, nE, dof_to_corner);

  out.face_local.assign(static_cast<size_t>(shell.face_count()), -1);
  std::vector<asawa::shell::FaceId> local_face;
  for (asawa::shell::FaceId fi : shell.get_face_range(true)) {
    if (shell.fsize(fi) != 3)
      continue;
    cr_tri tri;
    bool filled = false;
    shell.const_for_each_face_tri(
        fi, [&](asawa::shell::CornerId c0, asawa::shell::CornerId c1,
                asawa::shell::CornerId c2, const asawa::shell::shell &Ms) {
          // The walk invokes the callback again after the triangle is closed.
          // The later call repeats a corner. The first triple is the face.
          if (filled)
            return;
          filled = true;
          const asawa::shell::CornerId c[3] = {c0, c1, c2};
          for (int i = 0; i < 3; ++i) {
            tri.v[static_cast<size_t>(i)] = static_cast<int>(Ms.vert(c[i]));
            const int slot = static_cast<int>(c[i]) / 2;
            tri.e[static_cast<size_t>(i)] =
                slot >= 0 && slot < static_cast<int>(out.slot_map.size())
                    ? out.slot_map[static_cast<size_t>(slot)]
                    : -1;
          }
        });
    out.face_local[static_cast<size_t>(static_cast<int>(fi))] =
        static_cast<int>(out.tris.size());
    local_face.push_back(fi);
    out.tris.push_back(std::move(tri));
  }
  for (int i = 0; i < static_cast<int>(out.tris.size()); ++i) {
    for (asawa::shell::FaceId fj :
         shell.face_one_ring_face_ids(local_face[static_cast<size_t>(i)])) {
      const int lj = static_cast<int>(fj);
      if (lj < 0 || lj >= static_cast<int>(out.face_local.size()))
        continue;
      const int local = out.face_local[static_cast<size_t>(lj)];
      if (local >= 0)
        out.tris[static_cast<size_t>(i)].stencil.push_back(local);
    }
  }

  std::vector<std::vector<asawa::shell::FaceId>> incident;
  kusama::build_edge_incident_faces(shell, out.slot_map, nE, incident);
  std::vector<std::vector<kusama::edge_triangle_link>> adj;
  kusama::build_edge_triangle_adjacency(shell, out.slot_map, nE, adj);

  out.edges.assign(static_cast<size_t>(nE), cr_edge{});
  for (int e = 0; e < nE; ++e) {
    cr_edge &ed = out.edges[static_cast<size_t>(e)];
    const asawa::shell::CornerId c = dof_to_corner[static_cast<size_t>(e)];
    if (static_cast<int>(c) >= 0) {
      ed.v0 = static_cast<int>(shell.vert(c));
      ed.v1 = static_cast<int>(shell.vert(shell.next(c)));
    }
    for (asawa::shell::FaceId fi : incident[static_cast<size_t>(e)]) {
      const int lj = static_cast<int>(fi);
      if (lj < 0 || lj >= static_cast<int>(out.face_local.size()))
        continue;
      const int local = out.face_local[static_cast<size_t>(lj)];
      if (local >= 0)
        ed.faces.push_back(local);
    }
    for (const kusama::edge_triangle_link &lk : adj[static_cast<size_t>(e)]) {
      const int lj = static_cast<int>(lk.face);
      if (lj < 0 || lj >= static_cast<int>(out.face_local.size()))
        continue;
      const int local = out.face_local[static_cast<size_t>(lj)];
      if (local >= 0)
        ed.nbr.emplace_back(lk.other_edge, local);
    }
  }
  return true;
}

/// Incident edge dofs and local face indices for one vertex.
inline void stamp_cr_vertex(const asawa::shell::shell &shell,
                            const cr_complex &cx, int vertex,
                            std::vector<int> &faces, std::vector<int> &edges) {
  faces.clear();
  edges.clear();
  if (vertex < 0 || vertex >= shell.vert_count())
    return;
  const asawa::shell::VertId vi = asawa::shell::vert_id(vertex);
  if (static_cast<int>(shell.vbegin(vi)) < 0)
    return;
  shell.const_for_each_vertex(
      vi, [&](asawa::shell::CornerId ci, const asawa::shell::shell &Ms) {
        const int slot = static_cast<int>(ci) / 2;
        if (slot >= 0 && slot < static_cast<int>(cx.slot_map.size())) {
          const int e = cx.slot_map[static_cast<size_t>(slot)];
          if (e >= 0)
            edges.push_back(e);
        }
        const int lj = static_cast<int>(Ms.face(ci));
        if (lj >= 0 && lj < static_cast<int>(cx.face_local.size())) {
          const int local = cx.face_local[static_cast<size_t>(lj)];
          if (local >= 0)
            faces.push_back(local);
        }
      });
}

/// Copy one vertex fan out of a complex built in `init_*`.
/// `support` is the face across each link edge. That edge's derivative uses
/// both triangles, and the far one is not in the vertex fan.
inline void fan_from_ring(const cr_complex &cx, const std::vector<int> &face_ids,
                          const std::vector<int> &edge_ids,
                          std::vector<fan_face> &faces, std::vector<fan_face> &support,
                          std::vector<fan_edge> &edges) {
  faces.clear();
  support.clear();
  edges.clear();
  faces.reserve(face_ids.size());
  std::vector<char> have_face(cx.tris.size(), 0);
  auto push_face = [&](int f, std::vector<fan_face> &dst) {
    if (f < 0 || f >= static_cast<int>(cx.tris.size()))
      return;
    if (have_face[static_cast<size_t>(f)])
      return;
    have_face[static_cast<size_t>(f)] = 1;
    const cr_tri &t = cx.tris[static_cast<size_t>(f)];
    fan_face face;
    face.v = t.v;
    face.e = t.e;
    dst.push_back(face);
  };
  for (const int f : face_ids)
    push_face(f, faces);
  auto push_edge = [&](int dof) {
    if (dof < 0 || dof >= static_cast<int>(cx.edges.size()))
      return;
    for (const fan_edge &ed : edges) {
      if (ed.dof == dof)
        return;
    }
    const cr_edge &src = cx.edges[static_cast<size_t>(dof)];
    edges.push_back(fan_edge{dof, src.v0, src.v1});
  };
  for (const int e : edge_ids)
    push_edge(e);
  for (const fan_face &face : faces) {
    for (int i = 0; i < 3; ++i) {
      const int dof = face.e[static_cast<size_t>(i)];
      push_edge(dof);
      if (dof < 0 || dof >= static_cast<int>(cx.edges.size()))
        continue;
      for (const int f : cx.edges[static_cast<size_t>(dof)].faces)
        push_face(f, support);
    }
  }
  for (const fan_face &face : support) {
    for (int i = 0; i < 3; ++i)
      push_edge(face.e[static_cast<size_t>(i)]);
  }
}

inline std::vector<dirichlet_face> dirichlet_faces_of(const cr_complex &cx, int edge) {
  std::vector<dirichlet_face> faces;
  if (edge < 0 || edge >= static_cast<int>(cx.edges.size()))
    return faces;
  for (const int f : cx.edges[static_cast<size_t>(edge)].faces) {
    if (f < 0 || f >= static_cast<int>(cx.tris.size()))
      continue;
    dirichlet_face face;
    face.v = cx.tris[static_cast<size_t>(f)].v;
    for (const int s : cx.tris[static_cast<size_t>(f)].stencil) {
      if (s < 0 || s >= static_cast<int>(cx.tris.size()))
        continue;
      face.stencil.push_back(cx.tris[static_cast<size_t>(s)].v);
    }
    faces.push_back(std::move(face));
  }
  return faces;
}

inline cr_coeff_block::ptr cr_coeff_of(const std::vector<sim_block::ptr> &blocks,
                                       int *index = nullptr) {
  for (int i = 0; i < static_cast<int>(blocks.size()); ++i) {
    auto cr = std::dynamic_pointer_cast<cr_coeff_block>(blocks[static_cast<size_t>(i)]);
    if (cr) {
      if (index)
        *index = i;
      return cr;
    }
  }
  return nullptr;
}

inline void init_cr_dirichlet(asawa::shell::shell &shell,
                             std::vector<projection_constraint::ptr> &constraints,
                             const std::vector<vec3> &x, const real &w,
                             real lambda, std::vector<sim_block::ptr> blocks) {
  if (!(w > real(0)))
    return;
  if (!shell.verts_are_dense_packed())
    return;
  if (x.size() != static_cast<size_t>(shell.vert_count()))
    return;
  auto cr = cr_coeff_of(blocks);
  if (!cr || cr->nE <= 0)
    return;

  cr_complex cx;
  if (!build_cr_complex(shell, cx) ||
      static_cast<int>(cx.edges.size()) != cr->nE)
    return;
  const int nE = cr->nE;

  int nL = 0;
  const Eigen::SparseMatrix<real> Lmat =
      kusama::laplace::vector_dirichlet::energy(shell, x, &nL);
  if (nL != nE || Lmat.rows() != 2 * nE)
    return;
  std::vector<trip> L;
  L.reserve(static_cast<size_t>(Lmat.nonZeros()));
  for (int k = 0; k < Lmat.outerSize(); ++k) {
    for (Eigen::SparseMatrix<real>::InnerIterator it(Lmat, k); it; ++it)
      L.emplace_back(static_cast<int>(it.row()), static_cast<int>(it.col()),
                     it.value());
  }

  std::vector<real> mass = cr->mass;
  if (static_cast<int>(mass.size()) != 2 * nE)
    mass.assign(static_cast<size_t>(2 * nE), real(1e-8));

  int coeff = 1;
  cr_coeff_of(blocks, &coeff);
  std::vector<sim_block::ptr> row_blocks{blocks[0], blocks[static_cast<size_t>(coeff)]};
  std::vector<std::vector<std::pair<int, real>>> rows(static_cast<size_t>(2 * nE));
  for (const trip &t : L) {
    if (t.row() >= 0 && t.row() < 2 * nE)
      rows[static_cast<size_t>(t.row())].emplace_back(t.col(), t.value());
  }
  for (int e = 0; e < nE; ++e) {
    const cr_edge &ed = cx.edges[static_cast<size_t>(e)];
    const real mp = e < static_cast<int>(mass.size()) ? mass[static_cast<size_t>(e)]
                                                       : real(0);
    const real mq = e + nE < static_cast<int>(mass.size())
                        ? mass[static_cast<size_t>(e + nE)]
                        : real(0);
    constraints.push_back(cr_dirichlet::create(
        e, ed.v0, ed.v1, nE, mp, mq, rows[static_cast<size_t>(e)],
        rows[static_cast<size_t>(e + nE)], dirichlet_faces_of(cx, e), w, lambda,
        row_blocks));
  }
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
