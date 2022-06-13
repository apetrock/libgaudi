
#ifndef __TWOMANIFOLD_COORDINATE_INTERFACE__
#define __TWOMANIFOLD_COORDINATE_INTERFACE__
#include "m2.hpp"

namespace m2 {
namespace ci {

// pre declare function
template <typename SPACE>
typename SPACE::coordinate_type normal(typename surf<SPACE>::face_ptr f);

///////////////////////
// get/set coordinates
///////////////////////

template <typename SPACE> void dirty(typename surf<SPACE>::vertex_ptr v) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;

  v->set_dirty();
  for_each_vertex<SPACE>(v, [](face_vertex_ptr fv) {
    if (fv->face())
      fv->face()->set_dirty();
    if (fv->edge())
      fv->edge()->set_dirty();
  });
}

template <typename SPACE>
typename SPACE::coordinate_type
get_coordinate(const typename surf<SPACE>::vertex_ptr &v) {
  return v->template get<typename SPACE::coordinate_type>(
      SPACE::vertex_index::COORDINATE);
}

template <typename SPACE>
void set_coordinate(typename SPACE::coordinate_type p,
                    typename surf<SPACE>::vertex_ptr v) {
  v->template set<typename SPACE::coordinate_type>(
      SPACE::vertex_index::COORDINATE, p);

  dirty<SPACE>(v);
}

template <typename SPACE>
typename SPACE::coordinate_type
get_coordinate(typename surf<SPACE>::face_vertex_ptr fv) {
  return get_coordinate<SPACE>(fv->vertex());
}

template <typename SPACE>
void set_coordinate(typename SPACE::coordinate_type p,
                    typename surf<SPACE>::face_vertex_ptr fv) {
  set_coordinate<SPACE>(p, fv->vertex());
}

template <typename SPACE>
void set_dirty(typename surf<SPACE>::face_vertex_ptr fv) {
  set_dirty(fv->vertex());
}
///////////////////////
// Face Vertex Operations
///////////////////////

template <typename SPACE>
typename SPACE::real angle(typename surf<SPACE>::face_vertex_ptr fv) {
  using coordinate_type = typename SPACE::coordinate_type;

  coordinate_type ci = get_coordinate<SPACE>(fv);
  coordinate_type ca = get_coordinate<SPACE>(fv->next());
  coordinate_type cb = get_coordinate<SPACE>(fv->prev());
  coordinate_type cai = ca - ci;
  coordinate_type cbi = cb - ci;
  return va::angle_from_vectors(cai, cbi);
}

template <typename SPACE>
typename SPACE::real cotan(typename surf<SPACE>::face_vertex_ptr fv) {
  using coordinate_type = typename SPACE::coordinate_type;

  coordinate_type cp = get_coordinate<SPACE>(fv->prev());
  coordinate_type c0 = get_coordinate<SPACE>(fv);
  coordinate_type cn = get_coordinate<SPACE>(fv->next());
  return va::cotan(c0, cp, cn);
}

template <typename SPACE>
typename SPACE::real abs_cotan(typename surf<SPACE>::face_vertex_ptr fv) {
  using coordinate_type = typename SPACE::coordinate_type;

  coordinate_type cp = get_coordinate<SPACE>(fv->prev());
  coordinate_type c0 = get_coordinate<SPACE>(fv);
  coordinate_type cn = get_coordinate<SPACE>(fv->next());
  return va::abs_cotan(c0, cp, cn);
}

template <typename SPACE>
typename surf<SPACE>::coordinate_type
directed_vector(typename surf<SPACE>::face_vertex_ptr fv) {
  using coordinate_type = typename SPACE::coordinate_type;
  coordinate_type c0 = get_coordinate<SPACE>(fv);
  coordinate_type c1 = get_coordinate<SPACE>(fv->next());
  return c0 - c1;
}
///////////////////////
// Edge Operations
///////////////////////

template <typename SPACE>
typename SPACE::coordinate_type dir(typename surf<SPACE>::edge_ptr e) {
  using coordinate_type = typename SPACE::coordinate_type;
  coordinate_type c0 = get_coordinate<SPACE>(e->v1());
  coordinate_type c1 = get_coordinate<SPACE>(e->v2());
  return c0 - c1;
}

template <typename SPACE>
typename SPACE::real length(typename surf<SPACE>::edge_ptr e) {
  return dir<SPACE>(e).norm();
}

template <typename SPACE>
typename SPACE::coordinate_type normal(typename surf<SPACE>::edge_ptr e) {
  using coordinate_type = typename SPACE::coordinate_type;
  coordinate_type n1 = normal<SPACE>(e->v1()->face());
  coordinate_type n2 = normal<SPACE>(e->v2()->face());
  coordinate_type na = n1 + n2;
  na.normalize();
  return na;
}

template <typename SPACE>
typename SPACE::real cotan(typename surf<SPACE>::edge_ptr e) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;
  face_vertex_ptr fv1p = e->v1()->prev();
  face_vertex_ptr fv2p = e->v2()->prev();
  return cotan<SPACE>(fv1p) + cotan<SPACE>(fv2p);
}

template <typename SPACE>
typename surf<SPACE>::line_type get_line(typename surf<SPACE>::edge_ptr e) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;
  using line_type = typename surf<SPACE>::line_type;
  using coordinate_type = typename SPACE::coordinate_type;

  coordinate_type v0 = get_coordinate<SPACE>(e->v1());
  coordinate_type v1 = get_coordinate<SPACE>(e->v2());
  int i0 = e->v1()->vertex()->position_in_set();
  int i1 = e->v2()->vertex()->position_in_set();

  return line_type(v0, v1, i0, i1, e->position_in_set());
}

///////////////////////
// Face Operations
///////////////////////

template <typename SPACE>
typename SPACE::real area(typename surf<SPACE>::face_ptr f) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;
  using coordinate_type = typename SPACE::coordinate_type;

  typename SPACE::real a = 0.0;
  if (f->size() < 3)
    return 0.0;
  coordinate_type c0 = get_coordinate<SPACE>(f->fbegin());

  for_each_face<SPACE>(f, [&a, c0](face_vertex_ptr itb) {
    coordinate_type c1 = get_coordinate<SPACE>(itb);
    coordinate_type c2 = get_coordinate<SPACE>(itb->next());
    coordinate_type c10 = c1 - c0;
    coordinate_type c20 = c2 - c0;
    coordinate_type n = va::cross(c10, c20);
    a += n.norm();
    if (isnan(a)) {
      std::cout << c1.transpose() << "-" << c2.transpose() << std::endl;
      assert(!isnan(a));
    }
  });
  return 0.5 * a;
}

template <typename SPACE>
typename SPACE::coordinate_type
point_to_bary(typename surf<SPACE>::face_ptr f,
              typename SPACE::coordinate_type c) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;
  using coordinate_type = typename SPACE::coordinate_type;

  std::vector<coordinate_type> vertices;
  for_each_face<SPACE>(f, [&vertices](face_vertex_ptr itb) mutable {
    vertices.push_back(get_coordinate<SPACE>(itb));
  });
  return va::calc_bary(c, vertices);
  ;
}

template <typename SPACE>
typename SPACE::coordinate_type normal(typename surf<SPACE>::face_ptr f) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;
  using coordinate_type = typename SPACE::coordinate_type;

  coordinate_type n = coordinate_type(0, 0, 0);
  coordinate_type c0 = get_coordinate<SPACE>(f->fbegin());

  for_each_face_except_0<SPACE>(f, [&n, c0](face_vertex_ptr itb) {
    coordinate_type c1 = get_coordinate<SPACE>(itb);
    coordinate_type c2 = get_coordinate<SPACE>(itb->next());
    coordinate_type c10 = c1 - c0;
    coordinate_type c20 = c2 - c0;
    coordinate_type ni = va::cross(c10, c20);
    n += ni;
  });

  n.normalize();
  return n;
}

template <typename SPACE>
typename SPACE::coordinate_type center(typename surf<SPACE>::face_ptr f) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;
  using coordinate_type = typename SPACE::coordinate_type;

  coordinate_type cen = coordinate_type(0, 0, 0);
  typename SPACE::real n = 0.0;
  for_each_face<SPACE>(f, [&n, &cen](face_vertex_ptr itb) {
    cen += get_coordinate<SPACE>(itb);
    // std::cout << "      <<" << get_coordinate<SPACE>(itb).transpose()
    //          << std::endl;
    n += 1.0;
  });

  // std::cout << "      <<" << n << std::endl;

  cen *= 1.0 / n;
  return cen;
}

template <typename SPACE>
typename SPACE::box_type bound(typename surf<SPACE>::face_ptr f) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;
  using coordinate_type = typename SPACE::coordinate_type;

  coordinate_type min = center(f), max = min;
  for_each_face<SPACE>(f, [&min, &max](face_vertex_ptr itb) {
    coordinate_type p = get_coordinate(itb);
    max = va::max(p, max);
    min = va::min(p, max);
  });

  return makeBoxMinMax(min, max);
}

template <typename SPACE>
typename SPACE::real thinness(typename surf<SPACE>::face_ptr f) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;
  using T = typename SPACE::face_vertex_ptr;

  face_vertex_ptr fv0 = f->fbegin();
  face_vertex_ptr fv1 = fv0->next();
  face_vertex_ptr fv2 = fv1->next();
  T d0 = length(fv0->edge());
  T d1 = length(fv1->edge());
  T d2 = length(fv2->edge());
  T tEps = 0.95;
  T thin0 = d0 / (d1 + d2);
  T thin1 = d1 / (d0 + d2);
  T thin2 = d2 / (d0 + d1);
  T thin = thin0 > thin1 ? thin0 : thin1;
  thin = thin > thin2 ? thin : thin2;
  return thin;
}

template <typename SPACE>
std::vector<typename surf<SPACE>::triangle_type>
get_tris(typename surf<SPACE>::face_ptr f) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;
  using triangle_type = typename surf<SPACE>::triangle_type;
  using coordinate_type = typename SPACE::coordinate_type;

  std::vector<triangle_type> tris;

  face_vertex_ptr itb = f->fbegin();
  face_vertex_ptr c0 = f->fbegin();
  face_vertex_ptr c1 = c0->next();
  face_vertex_ptr c2 = c1->next();
  face_vertex_ptr ite = f->fbegin();

  coordinate_type v0 = get_coordinate<SPACE>(c0);
  int i0 = c0->vertex()->position_in_set();
  int fid = f->position_in_set();
  // std::cout << "face size: " << mSize << std::endl;
  bool at_head = false;
  int k = 0;
  while (!at_head) {

    if (k > 10)
      std::cout << k << " " << c0->vertex() << " " << c1->vertex()
                << c2->vertex() << " | ";
    if (k > 20)
      break;

    at_head = ite == c2;
    coordinate_type v1 = get_coordinate<SPACE>(c1);
    coordinate_type v2 = get_coordinate<SPACE>(c2);
    int i1 = c1->vertex()->position_in_set();
    int i2 = c2->vertex()->position_in_set();
    c1 = c2;
    c2 = c2->next();

    if ((v1 - v0).norm() == 0)
      continue;
    if ((v2 - v0).norm() == 0)
      continue;

    tris.push_back(triangle_type(v0, v1, v2, i0, i1, i2, fid));
    k++;
  }
  if (k > 10)
    std::cout << endl;
  // std::cout << "tri size: " << tris.size() << std::endl;
  return tris;
}

///////////////////////
// Vertex Operations
///////////////////////

template <typename SPACE>
typename SPACE::real norm(typename surf<SPACE>::vertex_ptr v0,
                          typename surf<SPACE>::vertex_ptr v1) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;
  using coordinate_type = typename SPACE::coordinate_type;
  using T = typename SPACE::real;

  coordinate_type c0 = get_coordinate<SPACE>(v0);
  coordinate_type c1 = get_coordinate<SPACE>(v1);
  return m2::va::norm(coordinate_type(c0 - c1));
}

template <typename SPACE>
typename SPACE::coordinate_type normal(typename surf<SPACE>::vertex_ptr v) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;
  using coordinate_type = typename SPACE::coordinate_type;
  using T = typename SPACE::real;

  coordinate_type N(0, 0, 0);
  for_each_vertex<SPACE>(v, [&N](face_vertex_ptr itb) {
    coordinate_type Ni = normal<SPACE>(itb->face()) * area<SPACE>(itb->face());
    N += Ni;
  });
  N.normalize();
  return N;
}

template <typename SPACE>
typename SPACE::real area(typename surf<SPACE>::vertex_ptr v) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;
  using real = typename SPACE::real;
  using coordinate_type = typename SPACE::coordinate_type;
  using T = typename SPACE::real;

  real a = 0.0;
  for_each_vertex<SPACE>(v, [&a](face_vertex_ptr itb) {
    a += 0.333 * area<SPACE>(itb->face()); // assuming triangle...
  });
  return a;
}

template <typename SPACE>
typename SPACE::real thinness(typename surf<SPACE>::face_vertex_ptr v) {
  // this isn't a terrible idea, but should get replaced with an SVD
  using vertex_ptr = typename surf<SPACE>::vertex_ptr;
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;
  using T = typename SPACE::real;

  // this would be better with an svd?
  if (v->size() == 0)
    return 0.0;

  face_vertex_ptr fvb = v->fbegin();
  face_vertex_ptr fve = v->fend();
  bool iterating = true;
  T vthin = 0;
  int s = 0;
  while (iterating) {
    iterating = fvb != fve;
    vthin += fvb->face()->thinness();
    fvb = fvb->vnext();
    s++;
  }
  vthin /= (T)s;

  return vthin;
}

template <typename SPACE>
typename surf<SPACE>::coordinate_array
coordinate_trace(typename surf<SPACE>::face_ptr f) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;

  typename surf<SPACE>::coordinate_array array_out;
  for_each_face<SPACE>(f, [array_out](face_vertex_ptr itb) {
    array_out.push_back(get_coordinate(itb));
  });

  return array_out;
}

template <typename SPACE>
typename surf<SPACE>::coordinate_array
flagged_coordinate_trace(typename surf<SPACE>::face_ptr f,
                         unsigned int flag_num) {
  using face_vertex_ptr = typename surf<SPACE>::face_vertex_ptr;

  typename surf<SPACE>::coordinate_array array_out;
  for_each_face<SPACE>(f, [array_out, flag_num](face_vertex_ptr itb) {
    int cflag = itb->vertex()->flag;
    if (cflag == flag_num) {
      array_out.push_back(get_coordinate(itb));
    }
  });
  return array_out;
}

///////////////////////
// Surf Operations
///////////////////////

template <typename SPACE>
typename SPACE::real geometric_mean_length(typename surf<SPACE>::surf_ptr s) {
  using real = typename SPACE::real;
  real sum = 0;
  int N = 0;
  for (auto e : s->get_edges()) {
    if (!e)
      continue;
    sum += log(m2::ci::length<SPACE>(e));
    N++;
  }
  sum /= real(N);
  return exp(sum);
}

template <typename SPACE>
typename SPACE::real mean_length(typename surf<SPACE>::surf_ptr s) {
  using real = typename SPACE::real;
  real sum = 0;
  int N = 0;
  for (auto e : s->get_edges()) {
    real dist = m2::ci::length<SPACE>(e);
    if (!e)
      continue;
    // std::cout << N << " " << dist << std::endl;
    if (dist == 0)
      std::cout << "zero" << std::endl;
    sum += m2::ci::length<SPACE>(e);
    N++;
  }

  // std::cout << sum << " " << N << " " << sum / real(N) << std::endl;
  sum /= real(N);
  return sum;
}

template <typename SPACE>
typename SPACE::real area(typename surf<SPACE>::surf_ptr s) {
  using real = typename SPACE::real;
  real ar = 0;
  auto verts = s->get_faces();
  for (auto f : s->get_faces()) {
    if (!f)
      continue;
    ar += area<SPACE>(f);
  }
  return ar;
}

template <typename SPACE>
typename SPACE::coordinate_type center(typename surf<SPACE>::surf_ptr s) {
  using T = typename SPACE::real;
  using coordinate_type = typename SPACE::coordinate_type;
  typename SPACE::coordinate_type avg(0, 0, 0);

  auto verts = s->get_vertices();
  for (int i = 0; i < verts.size(); i++) {
    if (!verts[i])
      continue;
    avg += get_coordinate<SPACE>(verts[i]);
  }

  avg /= (T)verts.size();
  return avg;
}

template <typename SPACE>
typename SPACE::box_type bound(typename surf<SPACE>::surf_ptr s) {
  using T = typename SPACE::real;
  using coordinate_type = typename SPACE::coordinate_type;

  coordinate_type max = center<SPACE>(s);
  coordinate_type min = max;
  auto verts = s->get_vertices();
  for (int i = 0; i < verts.size(); i++) {
    coordinate_type p = get_coordinate<SPACE>(verts[i]);
    for (int j = 0; j < 3; j++) {
      max[j] = max[j] > p[j] ? max[j] : p[j];
      min[j] = min[j] < p[j] ? min[j] : p[j];
    }
  }

  return makeBoxMinMax<T, coordinate_type>(min, max);
}

template <typename SPACE, typename TYPE>
std::vector<TYPE> get(m2::surf<SPACE> *surf, typename SPACE::vertex_index id) {
  auto vertices = surf->get_vertices();
  std::vector<TYPE> vals(surf->get_vertices().size());
  int i = 0;
  for (auto v : surf->get_vertices()) {
    vals[i] = v->template get<TYPE>(id);
    i++;
  }
  return vals;
}

template <typename SPACE, typename TYPE>
void set(m2::surf<SPACE> *surf, std::vector<TYPE> vals,
           typename SPACE::vertex_index id) {
  int i = 0;
  for (auto v : surf->get_vertices()) {
    v->template set<TYPE>(id, vals[i]);
    i++;
  }
}

template <typename SPACE>
std::vector<typename SPACE::coordinate_type>
get_coordinates(typename surf<SPACE>::surf_ptr s) {

  std::vector<typename SPACE::coordinate_type> coords;

  auto verts = s->get_vertices();
  int i = 0;
  for (auto v : s->get_vertices()) {
    if (!s->has_vertex(i++))
      continue;
    coords.push_back(get_coordinate<SPACE>(v));
  }

  return coords;
}

template <typename SPACE>
void set_coordinates(const std::vector<typename SPACE::coordinate_type> &coords,
                     typename surf<SPACE>::surf_ptr s) {
  using T = typename SPACE::real;

  auto verts = s->get_vertices();
  int i = 0;
  int j = 0;

  for (auto c : coords) {
    if (!s->has_vertex(j++))
      continue;

    set_coordinate<SPACE>(c, verts[i++]);
  }
}

template <typename SPACE>
std::vector<typename SPACE::coordinate_type>
get_centers(typename surf<SPACE>::surf_ptr s) {
  using T = typename SPACE::real;
  using coordinate_type = typename SPACE::coordinate_type;
  std::vector<coordinate_type> coords;
  for (auto f : s->get_faces()) {
    coords.push_back(center<SPACE>(f));
  }
  return coords;
}

template <typename SPACE>
std::vector<typename SPACE::coordinate_type>
get_vertex_normals(typename surf<SPACE>::surf_ptr s) {

  std::vector<typename SPACE::coordinate_type> coords;

  auto verts = s->get_vertices();
  for (auto v : s->get_vertices()) {
    coords.push_back(normal<SPACE>(v));
  }

  return coords;
}

template <typename SPACE>
std::vector<typename surf<SPACE>::line_type>
get_lines(std::vector<typename surf<SPACE>::edge_ptr> edges) {
  using line_type = typename surf<SPACE>::line_type;
  std::vector<line_type> lines;

  for (auto e : edges) {
    if (!e)
      continue;
    line_type line = m2::ci::get_line<SPACE>(e);
    lines.push_back(line);
  }
  return lines;
}

template <typename SPACE>
std::vector<typename surf<SPACE>::triangle_type>
get_lines(typename surf<SPACE>::surf_ptr s) {
  return get_lines(s->get_edges());
}

template <typename SPACE>
std::vector<typename surf<SPACE>::triangle_type>
get_tris(std::vector<typename surf<SPACE>::face_ptr> faces) {
  using triangle_type = typename surf<SPACE>::triangle_type;
  std::vector<triangle_type> tris;

  for (auto f : faces) {
    if (!f)
      continue;
    std::vector<triangle_type> ftris = m2::ci::get_tris<SPACE>(f);
    tris.insert(tris.end(), ftris.begin(), ftris.end());
  }
  return tris;
}

template <typename SPACE, typename V>
std::vector<V> verts_to_faces(std::vector<V> vert_vals,
                              typename surf<SPACE>::surf_ptr s) {
  M2_TYPEDEFS;

  auto &faces = s->get_faces();
  std::vector<V> face_vals(faces.size());

  for (int i = 0; i < faces.size(); i++) {
    if (!s->has_face(i))
      continue;
    if (faces[i]->size() < 3)
      continue;

    V r = z::zero<V>();

    m2::for_each_face<SPACE>(faces[i], [&r, &vert_vals](face_vertex_ptr fv) {
      int j = fv->vertex()->position_in_set();
      real l = fv->template get<real>(SPACE::face_vertex_index::BARY);
      r += l * vert_vals[j];
    });
    // std::cout << r << " " << m2::va::norm(r) << std::endl;
    face_vals[i] = r;
  }
  return face_vals;
}

template <typename SPACE, typename V>
std::vector<V> faces_to_verts(std::vector<V> face_vals,
                              typename surf<SPACE>::surf_ptr s) {
  M2_TYPEDEFS;

  auto &verts = s->get_vertices();
  std::vector<V> vert_vals(verts.size());

  for (int i = 0; i < verts.size(); i++) {
    if (!s->has_vertex(i))
      continue;

    V r = z::zero<V>();

    m2::for_each_vertex<SPACE>(verts[i], [&r, &face_vals](face_vertex_ptr fv) {
      int j = fv->face()->position_in_set();
      real l = fv->template get<real>(SPACE::face_vertex_index::BARY);
      r += l * face_vals[j];
    });
    // std::cout << r << " " << m2::va::norm(r) << std::endl;
    
    vert_vals[i] = r;
  }
  return vert_vals;
}

} // namespace ci

template <typename SPACE> class default_interface {
  M2_TYPEDEFS;

public:
  coordinate_type coordinate(vertex_ptr fv) {
    return ci::get_coordinate<SPACE>(fv);
  }

  void coordinate(typename SPACE::coordinate_type p, vertex_ptr fv) {
    return ci::set_coordinate<SPACE>(p, fv);
  }

  coordinate_type coordinate(face_vertex_ptr fv) {
    return ci::get_coordinate<SPACE>(fv);
  }

  void coordinate(typename SPACE::coordinate_type p, face_vertex_ptr fv) {
    return ci::set_coordinate<SPACE>(p, fv);
  }

  coordinate_type center(face_ptr f) { return ci::center<SPACE>(f); }
  coordinate_type normal(face_ptr f) { return ci::normal<SPACE>(f); }
  box_type bound(face_ptr f) { return ci::bound<SPACE>(f); }
  T area(face_ptr f) { return ci::area<SPACE>(f); }

  T cotan(face_vertex_ptr fv) { return ci::cotan<SPACE>(fv); }

  T cotan(edge_ptr e) { return ci::cotan<SPACE>(e); }
  T length(edge_ptr e) { return ci::length<SPACE>(e); }
  coordinate_type dir(edge_ptr e) { return ci::dir<SPACE>(e); }

  box_type bound(surf_ptr s) { return ci::bound<SPACE>(s); }
};

} // namespace m2

#endif