#include "m2.hpp"
#include <zlib.h>

#ifndef __TWOMANIFOLDFLAT__
#define __TWOMANIFOLDFLAT__

namespace asawa {
/////////////////////////////////////////////////////////////////////////////
// flat
/////////////////////////////////////////////////////////////////////////////

template <typename SPACE, typename ITYPE> struct flattened_data_vector {
public:
  flattened_data_vector(){};
  virtual ~flattened_data_vector(){};

  virtual void push_data(
      std::vector<typename surf<SPACE>::template data_node<ITYPE> *> nodes,
      const ITYPE &i) = 0;

  virtual void apply_data(
      std::vector<typename surf<SPACE>::template data_node<ITYPE> *> nodes,
      const ITYPE &i) = 0;

  virtual void write(FILE *file) const = 0;
  virtual void read(FILE *file) = 0;
  virtual void clear() = 0;
};

template <typename SPACE, typename ITYPE, typename TYPE>
struct flattened_data_vector_t : public flattened_data_vector<SPACE, ITYPE> {
public:
  flattened_data_vector_t(){};
  virtual ~flattened_data_vector_t(){};

  virtual void push_data(
      std::vector<typename surf<SPACE>::template data_node<ITYPE> *> nodes,
      const ITYPE &i) {
    for (auto n : nodes) {
      TYPE data = n->template get<TYPE>(i);
      _data.push_back(data);
    }
  };

  virtual void apply_data(
      std::vector<typename surf<SPACE>::template data_node<ITYPE> *> nodes,
      const ITYPE &i) {
    int idx = 0;
    for (auto n : nodes) {
      TYPE data = _data[idx++];
      n->template set<TYPE>(i, data);
    }
  };

  virtual void write(FILE *file) const {
    size_t nData = _data.size();
    size_t e;
    e = fwrite((void *)&nData, sizeof(size_t), 1, file);
    e = fwrite((void *)_data.data(), sizeof(TYPE), nData, file);
  }

  virtual void read(FILE *file) {
    size_t nData;
    size_t e;
    e = fread((void *)&nData, sizeof(size_t), 1, file);
    _data.resize(nData);
    e = fread((void *)_data.data(), sizeof(TYPE), nData, file);
  }

  virtual void clear() { _data.clear(); }

  vector<TYPE> _data;
};

template <typename SPACE, typename ITYPE> struct flattened_data {
public:
  int _size = static_cast<int>(ITYPE::MAXINDEX);
  flattened_data() {}

  ~flattened_data() {
    for (auto d : _data) {
      delete d;
    }
  }

  virtual void init_data() {
    using flat_data = flattened_data_vector<SPACE, ITYPE>;
    for (int i = 0; i < static_cast<int>(ITYPE::MAXINDEX); i++) {
      ITYPE ii = static_cast<ITYPE>(i);

      storage_type type = SPACE::get_type(ii);

      flat_data *data;
      switch (type) {
      case storage_type::REAL:
        data =
            new flattened_data_vector_t<SPACE, ITYPE, typename SPACE::real>();
        break;
      case storage_type::VEC3:
        data = new flattened_data_vector_t<SPACE, ITYPE,
                                           typename SPACE::coordinate_type>();
        break;
      default:;
        // do nothing
      }
      _data.push_back(data);
    }
  }

  virtual void
  flatten(const std::vector<typename surf<SPACE>::template data_node<ITYPE> *>
              &nodes) {
    using flat_data = flattened_data_vector<SPACE, ITYPE>;
    _data = std::vector<flat_data *>();
    this->init_data();
    for (int i = 0; i < static_cast<int>(ITYPE::MAXINDEX); i++) {
      ITYPE ii = static_cast<ITYPE>(i);
      _data[i]->push_data(nodes, ii);
    }
  };

  void flatten(const std::vector<typename surf<SPACE>::face_ptr> &faces) {
    std::vector<typename asawa::surf<SPACE>::template data_node<
        typename SPACE::face_index> *>
        face_nodes;
    for (auto f : faces)
      face_nodes.push_back(f);
    std::cout << " flatten face" << std::endl;
    this->flatten(face_nodes);
  };

  void flatten(const std::vector<typename surf<SPACE>::edge_ptr> &edges) {
    std::vector<typename asawa::surf<SPACE>::template data_node<
        typename SPACE::edge_index> *>
        edge_nodes;
    for (auto e : edges)
      edge_nodes.push_back(e);
    this->flatten(edge_nodes);
  };

  void flatten(const std::vector<typename surf<SPACE>::vertex_ptr> &vertices) {

    std::vector<typename asawa::surf<SPACE>::template data_node<
        typename SPACE::vertex_index> *>
        vertex_nodes;
    for (auto v : vertices)
      vertex_nodes.push_back(v);
    this->flatten(vertex_nodes);
  };

  virtual void
  inflate(const std::vector<typename surf<SPACE>::template data_node<ITYPE> *>
              &nodes) {

    for (int i = 0; i < static_cast<int>(ITYPE::MAXINDEX); i++) {
      ITYPE ii = static_cast<ITYPE>(i);
      flattened_data_vector<SPACE, ITYPE> *data = _data[i];
      data->apply_data(nodes, ii);
    }
  };

  void inflate(const std::vector<typename surf<SPACE>::face_ptr> &faces) {
    std::vector<typename asawa::surf<SPACE>::template data_node<
        typename SPACE::face_index> *>
        face_nodes;
    for (auto f : faces)
      face_nodes.push_back(f);

    std::cout << "face nodes: " << face_nodes.size() << std::endl;

    this->inflate(face_nodes);
  };

  void inflate(const std::vector<typename surf<SPACE>::edge_ptr> &edges) {
    std::vector<typename asawa::surf<SPACE>::template data_node<
        typename SPACE::edge_index> *>
        edge_nodes;

    for (auto e : edges)
      edge_nodes.push_back(e);

    this->inflate(edge_nodes);
  };

  void inflate(const std::vector<typename surf<SPACE>::vertex_ptr> &vertices) {

    std::vector<typename asawa::surf<SPACE>::template data_node<
        typename SPACE::vertex_index> *>
        vertex_nodes;
    for (auto v : vertices)
      vertex_nodes.push_back(v);

    this->inflate(vertex_nodes);
  };

  virtual void write(FILE *file) const {
    for (auto datum : _data)
      datum->write(file);
  }

  virtual void read(FILE *file) {
    this->init_data();
    for (auto datum : _data)
      datum->read(file);
  }

  virtual void clear() {
    for (auto datum : _data)
      datum->clear();
  }

  std::vector<flattened_data_vector<SPACE, ITYPE> *> _data;
};

template <typename SPACE> struct flattened_surf {
  // class used for serialization of surfaces and associated data
public:
  M2_TYPEDEFS;
  flattened_surf() {}

  flattened_surf(const surf_ptr surf) { this->from_surf(surf); }

  ~flattened_surf() {}

  void from_surf(const surf_ptr surf) {

    nFaces = surf->get_faces().size();
    nVertices = surf->get_vertices().size();
    nCorners = 2 * surf->get_edges().size();
    corner_verts = std::vector<size_t>(2 * surf->get_edges().size(), 0);
    corner_faces = std::vector<size_t>(2 * surf->get_edges().size(), 0);
    corner_next = std::vector<size_t>(2 * surf->get_edges().size(), 0);

    std::vector<edge_ptr> edges = surf->get_edges();

    std::vector<typename asawa::surf<SPACE>::template data_node<
        typename SPACE::face_vertex_index> *>
        corners;

    auto push_fv = [&](face_vertex_ptr fv) {
      size_t c = fv->position_in_flat_set();
      vertex_ptr v = fv->vertex();
      face_ptr f = fv->face();
      face_vertex_ptr fvn = fv->next();
      size_t cn = fvn->position_in_flat_set();

      corner_next[c] = cn;
      corner_faces[c] = f->position_in_set();
      corner_verts[c] = v->position_in_set();
      corners.push_back(fv);
    };

    for (auto e : edges) {
      push_fv(e->v1());
      push_fv(e->v2());
    }

    corner_data.flatten(corners);
    edge_data.flatten(surf->get_edges());
    face_data.flatten(surf->get_faces());
    vertex_data.flatten(surf->get_vertices());
  }

  surf_ptr to_surf() {
    surf_ptr surf = new surf_type();

    vertex_array verts;
    for (int i = 0; i < nVertices; i++) {
      vertex_ptr v = new vertex_type();
      verts.push_back(v);
      surf->push_vertex(v);
    }

    face_array faces;
    for (int i = 0; i < nFaces; i++) {
      face_ptr f = new face_type();
      faces.push_back(f);
      surf->push_face(f);
    }

    std::vector<typename asawa::surf<SPACE>::template data_node<
        typename SPACE::face_vertex_index> *>
        corners;

    std::vector<face_vertex_ptr> fvs;
    for (int i = 0; i < corner_verts.size(); i += 2) {
      face_vertex_ptr fv0 = new face_vertex_type();
      face_vertex_ptr fv1 = new face_vertex_type();

      edge_ptr e = new edge_type();
      e->set(fv0, fv1);
      surf->push_edge(e);
      verts[corner_verts[i]]->set_front(fv0);
      verts[corner_verts[i + 1]]->set_front(fv1);
      faces[corner_faces[i]]->set_front(fv0);
      faces[corner_faces[i + 1]]->set_front(fv1);
      fvs.push_back(fv0);
      fvs.push_back(fv1);
      corners.push_back(fv0);
      corners.push_back(fv1);
    }

    int i = 0;

    for (auto fv : fvs) {
      fv->set_next(fvs[corner_next[i++]]);
    }

    corner_data.inflate(corners);
    edge_data.inflate(surf->get_edges());
    face_data.inflate(surf->get_faces());
    vertex_data.inflate(surf->get_vertices());

    surf->update_all();
    return surf;
  }

  ///////////////////////////////////////////////////////////////////////
  // write out a field to a file stream
  ///////////////////////////////////////////////////////////////////////
  void clear() {
    nFaces = 0;
    nVertices = 0;
    corner_verts.clear();
    corner_faces.clear();
    corner_next.clear();

    corner_data.clear();
    edge_data.clear();
    face_data.clear();
    vertex_data.clear();
  }

  ///////////////////////////////////////////////////////////////////////
  // write out a field to a file stream
  ///////////////////////////////////////////////////////////////////////
  void write(FILE *file) const {
    size_t e;
    e = fwrite((void *)&nFaces, sizeof(size_t), 1, file);
    e = fwrite((void *)&nVertices, sizeof(size_t), 1, file);
    e = fwrite((void *)&nCorners, sizeof(size_t), 1, file);

    // always write out as a double

    e = fwrite((void *)corner_verts.data(), sizeof(size_t), nCorners, file);
    e = fwrite((void *)corner_faces.data(), sizeof(size_t), nCorners, file);
    e = fwrite((void *)corner_next.data(), sizeof(size_t), nCorners, file);
    corner_data.write(file);
    edge_data.write(file);
    face_data.write(file);
    vertex_data.write(file);
  }

  ///////////////////////////////////////////////////////////////////////
  // read in a field from a file stream
  ///////////////////////////////////////////////////////////////////////
  void read(FILE *file) {
    // read dimensions
    size_t e;
    e = fread((void *)&nFaces, sizeof(size_t), 1, file);
    e = fread((void *)&nVertices, sizeof(size_t), 1, file);
    e = fread((void *)&nCorners, sizeof(size_t), 1, file);
    corner_verts.resize(nCorners);
    corner_faces.resize(nCorners);
    corner_next.resize(nCorners);

    e = fread((void *)corner_verts.data(), sizeof(size_t), nCorners, file);
    e = fread((void *)corner_faces.data(), sizeof(size_t), nCorners, file);
    e = fread((void *)corner_next.data(), sizeof(size_t), nCorners, file);
    corner_data.read(file);
    edge_data.read(file);
    face_data.read(file);
    vertex_data.read(file);
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  void write(string filename) const {
    FILE *file;
    file = fopen(filename.c_str(), "wb");
    if (file == NULL) {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : "
           << endl;
      cout << " FIELD_3D write failed! " << endl;
      cout << " Could not open file " << filename.c_str() << endl;
      exit(0);
    }

    cout << " Writing file " << filename.c_str() << " ... ";
    flush(cout);

    // write to the stream
    write(file);

    // close the stream
    fclose(file);

    cout << " done. " << endl;
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  void read(string filename) {
    int size = filename.size();

    /*
    if (filename[size - 1] == 'z' && filename[size - 2] == 'g') {
      readGz(filename);
      return;
    }
    */

    FILE *file;
    file = fopen(filename.c_str(), "rb");
    if (file == NULL) {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : "
           << endl;
      cout << " FIELD_3D read failed! " << endl;
      cout << " Could not open file " << filename.c_str() << endl;
      exit(0);
    }
    // read from the stream
    read(file);

    // close the file
    fclose(file);
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  size_t nCorners = 0;
  size_t nFaces = 0;
  size_t nVertices = 0;
  std::vector<size_t> corner_verts;
  std::vector<size_t> corner_faces;
  std::vector<size_t> corner_next;

  flattened_data<SPACE, typename SPACE::face_vertex_index> corner_data;
  flattened_data<SPACE, typename SPACE::edge_index> edge_data;
  flattened_data<SPACE, typename SPACE::face_index> face_data;
  flattened_data<SPACE, typename SPACE::vertex_index> vertex_data;
};
} // namespace asawa
#endif