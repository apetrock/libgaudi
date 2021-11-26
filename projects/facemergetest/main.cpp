//#include "nanoguiincludes.h"

#if defined(WIN32)
#include <windows.h>
#endif

#include <stdio.h> /* defines FILENAME_MAX */
// #define WINDOWS  /* uncomment this line to use it for windows.*/
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include <iostream>
#include <limits>
#include <string>

#include "GaudiGraphics/buffers.hpp"
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"

#include "manifold/bins.hpp"
#include "manifold/m2Includes.h"
#include "manifold/m2Operators.h"
#include "manifold/make.hpp"
#include "manifold/objloader.hpp"

#include "manifold/moving_mesh.hpp"

#include "manifold/coordinate_interface.hpp"

#include "manifold/vec_addendum.h"

//#include "m2Operators.h"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;

using namespace GaudiMath;

namespace m2 {

template <typename SPACE> class triangle_ops {
public:
  M2_TYPEDEFS;

  std::vector<edge_ptr> getSharedEdges(face_ptr fA, face_ptr fB) {

    std::vector<edge_ptr> edges;

    face_vertex_ptr fvA = fA->fbegin();
    face_vertex_ptr fvAe = fA->fend();
    bool itA = true;
    while (itA) {
      itA = fvA != fvAe;
      edge_ptr eA = fvA->edge();
      face_vertex_ptr fvB = fB->fbegin();
      face_vertex_ptr fvBe = fB->fend();
      bool itB = true;
      while (itB) {
        itB = fvB != fvBe;
        fvB = fvB->next();
        edge_ptr eB = fvB->edge();
        if (eA == eB)
          edges.push_back(eA);
      }
      fvA = fvA->next();
    }
    return edges;
  }

  void joinFaceOneSharedEdge(face_ptr fA, face_ptr fB, surf_ptr mesh) {
    std::vector<edge_ptr> edges = getSharedEdges(fA, fB);
    if (edges.empty())
      return;
    edge_ptr e = edges[0];

    std::cout << "shared_edge: " << fA->position_in_set() << " "
              << fB->position_in_set() << std::endl;
    face_vertex_ptr fvA0 = e->v1();
    face_vertex_ptr fvA1 = fvA0->next();
    face_vertex_ptr fvA2 = fvA1->next();

    face_vertex_ptr fvB0 = e->v2();
    face_vertex_ptr fvB1 = fvB0->prev();
    face_vertex_ptr fvB2 = fvB1->prev();

    fvA0->vertex()->remove_face_vertex(fvA0);
    fvA1->vertex()->remove_face_vertex(fvA1);
    fvA2->vertex()->remove_face_vertex(fvA2);

    fvB0->vertex()->remove_face_vertex(fvB0);
    fvB1->vertex()->remove_face_vertex(fvB1);
    fvB2->vertex()->remove_face_vertex(fvB2);

    vertex_ptr vA0 = fvA0->vertex();
    vertex_ptr vA2 = fvA2->vertex();
    vertex_ptr vB0 = fvB0->vertex();
    vertex_ptr vB1 = fvB1->vertex();

    vA0->front() = fvA0->vprev();
    vB0->front() = fvB0->vprev();

    vA2->add_face_vertex(fvB2->coedge());
    vA2->update_all();

    edge_ptr eA1 = fvA1->edge();
    edge_ptr eA2 = fvA2->edge();
    edge_ptr eB1 = fvB1->edge();
    edge_ptr eB2 = fvB2->edge();

    eA1->set(fvA1->coedge(), fvB1->coedge());
    eA2->set(fvA2->coedge(), fvB2->coedge());

    delete fvA0;
    delete fvA1;
    delete fvA2;
    delete fvB0;
    delete fvB1;
    delete fvB2;

    mesh->remove_edge(e->position_in_set());
    mesh->remove_edge(eB1->position_in_set());
    mesh->remove_edge(eB2->position_in_set());
    mesh->remove_face(fA->position_in_set());
    mesh->remove_face(fB->position_in_set());
    mesh->remove_vertex(vB1->position_in_set());

    vA0->verify();
    vA2->verify();
    vB0->verify();
    eA1->v1()->face()->verify();
    eA1->v2()->face()->verify();
    eA2->v1()->face()->verify();
    eA2->v2()->face()->verify();
  }

  void joinFaceNoSharedEdges(face_ptr fA, face_ptr fB, surf_ptr mesh) {
    std::cout << "no shared_edge: " << fA->position_in_set() << " "
              << fB->position_in_set() << std::endl;

    face_vertex_ptr fvA0 = fA->fbegin();
    face_vertex_ptr fvA1 = fvA0->next();
    face_vertex_ptr fvA2 = fvA1->next();

    face_vertex_ptr fvB0 = fB->fbegin()->prev();
    face_vertex_ptr fvB1 = fvB0->prev();
    face_vertex_ptr fvB2 = fvB1->prev();

    vertex_ptr vA0 = fvA0->vertex();
    vertex_ptr vA1 = fvA1->vertex();
    vertex_ptr vA2 = fvA2->vertex();

    vertex_ptr vB0 = fvB0->vertex();
    vertex_ptr vB1 = fvB1->vertex();
    vertex_ptr vB2 = fvB2->vertex();

    vA0->remove_face_vertex(fvA0);
    vA1->remove_face_vertex(fvA1);
    vA2->remove_face_vertex(fvA2);

    vB0->remove_face_vertex(fvB0);
    vB1->remove_face_vertex(fvB1);
    vB2->remove_face_vertex(fvB2);

    vA0->add_face_vertex(fvB0->coedge());
    vA1->add_face_vertex(fvB1->coedge());
    vA2->add_face_vertex(fvB2->coedge());

    vA0->update_all();
    vA1->update_all();
    vA2->update_all();

    edge_ptr eA0 = fvA0->edge();
    edge_ptr eA1 = fvA1->edge();
    edge_ptr eA2 = fvA2->edge();

    edge_ptr eB0 = fvB0->edge();
    edge_ptr eB1 = fvB1->edge();
    edge_ptr eB2 = fvB2->edge();

    eA0->set(fvA0->coedge(), fvB0->coedge());
    eA1->set(fvA1->coedge(), fvB1->coedge());
    eA2->set(fvA2->coedge(), fvB2->coedge());

    delete fvA0;
    delete fvA1;
    delete fvA2;
    delete fvB0;
    delete fvB1;
    delete fvB2;

    mesh->remove_edge(eB0->position_in_set());
    mesh->remove_edge(eB1->position_in_set());
    mesh->remove_edge(eB2->position_in_set());

    mesh->remove_vertex(vB0->position_in_set());
    mesh->remove_vertex(vB1->position_in_set());
    mesh->remove_vertex(vB2->position_in_set());

    mesh->remove_face(fA->position_in_set());
    mesh->remove_face(fB->position_in_set());

    eA0->v1()->face()->verify();
    eA0->v2()->face()->verify();
    eA1->v1()->face()->verify();
    eA1->v2()->face()->verify();
    eA2->v1()->face()->verify();
    eA2->v2()->face()->verify();
  }

  void alignFaces(face_ptr fA, face_ptr fB) {
    fA->update_all();
    fB->update_all();

    face_vertex_ptr fvA0 = fA->fbegin();
    face_vertex_ptr fvA1 = fvA0->next();
    face_vertex_ptr fvA2 = fvA1->next();

    face_vertex_ptr fvB0 = fB->fbegin();
    face_vertex_ptr fvB1 = fvB0->prev();
    face_vertex_ptr fvB2 = fvB1->prev();

    fA->print_vert_ids();
    fB->print_vert_ids(true);

    coordinate_type ca0 = fvA0->coordinate();
    coordinate_type ca1 = fvA1->coordinate();
    coordinate_type ca2 = fvA2->coordinate();
    // assume the other coordinate rotates in opposite direction since they
    // are facing opposite directions
    coordinate_type cb0 = fvB0->coordinate();
    coordinate_type cb1 = fvB1->coordinate();
    coordinate_type cb2 = fvB2->coordinate();

    T d0 = 1.0 / 3.0 *
           ((cb0 - ca0).squaredNorm() + (cb1 - ca1).squaredNorm() +
            (cb2 - ca2).squaredNorm());
    T d1 = 1.0 / 3.0 *
           ((cb0 - ca1).squaredNorm() + (cb1 - ca2).squaredNorm() +
            (cb2 - ca0).norm());
    T d2 = 1.0 / 3.0 *
           ((cb0 - ca2).squaredNorm() + (cb1 - ca0).squaredNorm() +
            (cb2 - ca1).squaredNorm());

    std::cout << " distance: " << d0 << " " << d1 << " " << d2 << std::endl;

    face_vertex_ptr fvAb = fvA0;
    fvAb = (d1 < d0) && (d1 < d2) ? fvA1 : fvAb;
    fvAb = (d2 < d0) && (d2 < d1) ? fvA2 : fvAb;

    fA->front() = fvAb;
    fB->front() = fvB0;
  }

  void averageVerts(face_ptr fA, face_ptr fB) {
    face_vertex_ptr fvA = fA->fbegin();
    face_vertex_ptr fvAe = fA->fend();
    face_vertex_ptr fvB = fB->fbegin();
    bool itf = true;
    while (itf) {
      itf = fvA != fvAe;
      vertex_ptr vA = fvA->vertex();
      vertex_ptr vB = fvB->vertex();
      if (vA != vB) {
        coordinate_type cA = vA->coordinate();
        coordinate_type cB = vB->coordinate();
        vA->coordinate() = 0.5 * (cA + cB);
      }
      fvA = fvA->next();
      fvB = fvB->prev();
    }
  }

  void getNearest(surf_ptr mesh) {

    // using namespace nanogui;

    std::vector<face_ptr> faces = mesh->get_faces();
    face_ptr f = faces[0];

    std::vector<triangle_type> tris;

    std::for_each(faces.begin(), faces.end(), [&tris](const face_ptr &f) {
      std::vector<triangle_type> ftris = f->get_tris();
      tris.insert(tris.end(), ftris.begin(), ftris.end());
    });

    m2::aabb_tree<SPACE, triangle_type> tree(tris);

    std::cout << "tree.nodes.size(): " << tree.nodes.size() << std::endl;

    auto triDist = [](const triangle_type &A, const triangle_type &B) -> T {
      return A.dist(B);
    };

    T tol = 0.025;

    vector<triangle_pair> collected;
    // for (int i = 0; i < 1; i++) {
    // std::cout << " foo !" << std::endl;
    for (int i = 0; i < tris.size(); i++) {
      triangle_type triNear =
          m2::getNearest<SPACE, triangle_type, triangle_type>(
              tris[i], tree, tris, triDist, tol);
      if (triNear.faceId > -1) {
        collected.push_back(triangle_pair(tris[i], triNear));
      }
    }

    std::sort(collected.begin(), collected.end(), std::less<triangle_pair>());

    std::cout << " full list: " << collected.size() << std::endl;
    auto it = std::unique(collected.begin(), collected.end());

    collected.erase(it, collected.end());
    std::cout << " unique list: " << collected.size() << std::endl;

    bool trimming = true;
    while (trimming) {
      T d = collected.back().dist();
      if (d > tol)
        collected.pop_back();
      else
        trimming = false;
    }
    std::cout << " trimmed list: " << collected.size() << std::endl;

    collected.erase(std::remove_if(collected.begin(), collected.end(),
                                   [faces](auto p) {
                                     face_ptr fA = faces[p.A.faceId];
                                     face_ptr fB = faces[p.B.faceId];

                                     if (fA->flag > 0)
                                       return true;
                                     if (fB->flag > 0)
                                       return true;
                                     fA->flag = 1;
                                     fB->flag = 1;
                                     return false;
                                   }),
                    collected.end());

    std::cout << " flagged duplicates list: " << collected.size() << std::endl;

    for (int i = 0; i < collected.size(); i++) {
      // /for (int i = 0; i < 3; i++) {
      auto p = collected[i];
      face_ptr fA = faces[p.A.faceId];
      face_ptr fB = faces[p.B.faceId];

      std::cout << i << " " << fA->size() << " " << fB->size() << std::endl;

      int shared = fA->count_shared_vertices(fB);

      if (shared == 1) {
        std::cout << " shared == 1, punting!" << std::endl;
        continue;
      }
      std::cout << " shared: " << shared << std::endl;

      alignFaces(fA, fB);
      averageVerts(fA, fB);

      // continue;

      fA->print_vert_ids();
      fB->print_vert_ids(true);
      fA->print_shared_edge();
      fB->print_shared_edge();

      fA->flag_vertices();
      if (shared == 0) {
        int adjacent_shared = fA->count_adjacent_shared_vertices(fB);
        std::cout << " adjacent shared: " << shared << std::endl;
        if (adjacent_shared == 0)
          joinFaceNoSharedEdges(fA, fB, mesh);

      } else if (shared == 2)
        joinFaceOneSharedEdge(fA, fB, mesh);
    };

    mesh->verify();
    mesh->update_all();
    mesh->pack();
  }
};

} // namespace m2
template <typename SPACE>
m2::aabb_tree<SPACE, typename m2::face_triangle<SPACE>>
build_tree(m2::surf<SPACE> *mesh) {
  using namespace nanogui;
  M2_TYPEDEFS;
  std::cout << "  building tree" << std::endl;

  std::vector<face_ptr> faces = mesh->get_faces();
  std::vector<triangle_type> tris;

  std::for_each(faces.begin(), faces.end(), [&tris](const face_ptr &f) {
    std::vector<triangle_type> ftris = f->get_tris();
    tris.insert(tris.end(), ftris.begin(), ftris.end());
  });

  m2::aabb_tree<SPACE, triangle_type> tree(tris);

  using tree_type = m2::aabb_tree<SPACE, triangle_type>;
  std::cout << "tree.nodes.size(): " << tree.nodes.size() << std::endl;

  return tree;
};
template <typename SPACE>
void debugTree(
    const m2::aabb_tree<SPACE, typename m2::face_triangle<SPACE>> &tree,
    gg::DebugBufferPtr debug) {
  M2_TYPEDEFS;
  using namespace nanogui;
  for (int i = 0; i < tree.nodes.size(); i++) {
    auto &node = tree.nodes[i];
    if (node.size == 0)
      continue;

    Vec4 col(0.5, 0.5, 0.5, 1.0);
    box_type box = node.bbox;

    if (node.isLeaf()) {
      box.inflate(coordinate_type(0.11, 0.11, 0.11));
      coordinate_type half = box.half;
      coordinate_type cen = box.center;
      debug->pushBox(Vec4(cen[0], cen[1], cen[2], 1.0),
                     Vec4(half[0], half[1], half[2], 0.0),
                     Vec4(0.5, 0.75, 1.0, 1.0));
    }
  }

  debug->renderLines();
};

template <typename SPACE>
void debugLines(m2::surf<SPACE> *mesh, gg::DebugBufferPtr debug) {
  using namespace nanogui;
  M2_TYPEDEFS;
  std::cout << "  rendering vorticity" << std::endl;

  for (auto e : mesh->get_edges()) {
    coordinate_type c0 = e->v1()->coordinate();
    coordinate_type c1 = e->v2()->coordinate();
    coordinate_type N0 = e->v1()->face()->normal();
    coordinate_type N1 = e->v2()->face()->normal();
    c0 += 0.025 * N0;
    c1 += 0.025 * N1;

    gg::Vec4 col(0.5, 0.5, 0.5, 1.0);
    // gg::Vec4 col(0.5, 0.5, 0.5, 1.0);
    debug->pushLine(gg::Vec4(c0[0], c0[1], c0[2], 1.0),
                    gg::Vec4(c1[0], c1[1], c1[2], 1.0), col);
  }

  debug->renderLines();
};

template <typename SPACE>
void debugNormals(m2::surf<SPACE> *mesh, gg::DebugBufferPtr debug) {
  using namespace nanogui;
  M2_TYPEDEFS;
  std::cout << "  rendering vorticity" << std::endl;

  for (auto f : mesh->get_faces()) {
    coordinate_type c = f->center();
    coordinate_type N = f->normal();
    T area = f->calc_area();
    coordinate_type c0 = c;
    coordinate_type c1 = c + sqrt(area) * N;
    gg::Vec4 col(1.0, 0.0, 0.0, 1.0);
    // gg::Vec4 col(0.5, 0.5, 0.5, 1.0);
    debug->pushLine(gg::Vec4(c0[0], c0[1], c0[2], 1.0),
                    gg::Vec4(c1[0], c1[1], c1[2], 1.0), col);
  }

  debug->renderLines();
};

class FaceMergeScene;
using FaceMergeScenePtr = std::shared_ptr<FaceMergeScene>;

class FaceMergeScene : public gg::Scene {

public:
  static FaceMergeScenePtr create() {
    return std::make_shared<FaceMergeScene>();
  }

  FaceMergeScene() : gg::Scene() {
    using namespace nanogui;
    _obj = gg::BufferObject::create();
    _obj->init();
    mSceneObjects.push_back(_obj);

    _debugLines = gg::DebugBuffer::create();
    _debugLines->init();
    mSceneObjects.push_back(_debugLines);

    testMergeLow();
  }

  void testMergeLow() {

    m2::obj_loader<space3> load;
    m2::subdivide<space3> sub;
    m2::make<space3> mk;
    m2::convex_hull<space3> ch;
    m2::add_handle<space3> ah;
    m2::remesh<space3> rem;

    m2::construct<space3> bevel;
    m2::modify<space3> mod;

    _meshGraph = mk.two_tets(1.0, -0.095);

    _meshGraph = &sub.subdivide_control(*_meshGraph);
    _meshGraph = &sub.subdivide_control(*_meshGraph);
    _meshGraph = &sub.subdivide_control(*_meshGraph);
    _meshGraph = &sub.subdivide_control(*_meshGraph);
    std::vector<m2::surf<space3>::vertex_ptr> verts =
        _meshGraph->get_vertices();
    int start_s = clock();
    int i;
    space3::coordinate_type bar = verts[0]->coordinate();
    for (i = 0; i < 5000; i++) {
      std::for_each(
          verts.begin(), verts.end(), [](m2::surf<space3>::vertex_ptr v) {
            const space3::coordinate_type &bar = v->coordinate();
            m2::set_coordinate<space3>(bar, v);
            space3::coordinate_type foo = m2::get_coordinate<space3>(v);
            foo *= 1.01;
            m2::set_coordinate<space3>(foo, v);
          });
    }
    int stop_s = clock();
    cout << "time 0: " << (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000
         << std::endl;

    start_s = clock();
    i;
    for (i = 0; i < 5000; i++) {
      std::for_each(verts.begin(), verts.end(),
                    [](m2::surf<space3>::vertex_ptr v) {
                      const space3::coordinate_type &bar = v->coordinate();
                      v->coordinate() = bar;
                      space3::coordinate_type foo = v->coordinate();
                      foo *= 1.01;
                      v->coordinate() = foo;
                    });
    }
    stop_s = clock();
    cout << "time 1: " << (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000
         << std::endl;

    _meshGraph->print();
    rem.triangulate(_meshGraph);
    std::cout << "--center" << std::endl;
    std::cout << "--update_all" << std::endl;
    _meshGraph->update_all();
    std::cout << "--pack" << std::endl;
    _meshGraph->pack();
    std::cout << "--build sheet" << std::endl;

    mod.centerGeometry(*_meshGraph);
    std::cout << "creating buffer" << std::endl;
    m2::triangle_ops<space3> nearest;
    nearest.getNearest(_meshGraph);

    _meshGraph->update_all();
    _meshGraph->pack();
    mod.centerGeometry(*_meshGraph);

    std::cout << "creating buffer" << std::endl;

    // auto tree = build_tree(_meshGraph);
    // debugTree(tree, _debugLines);

    debugLines(_meshGraph, _debugLines);
    debugNormals(_meshGraph, _debugLines);
    gg::fillBuffer(_meshGraph, _obj);
  }

  virtual void onAnimate() { gg::fillBuffer(_meshGraph, _obj); }

  virtual void onDraw(gg::Viewer &viewer) {
    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });
    //_debugLines->clear();
  }

private:
  gg::DebugBufferPtr _debugLines = NULL;
  std::vector<gg::DrawablePtr> mSceneObjects;

  gg::BufferObjectPtr _obj = NULL;
  m2::surf<space3> *_meshGraph;
};

std::string GetCurrentWorkingDir(void) {
  char buff[FILENAME_MAX];
  GetCurrentDir(buff, FILENAME_MAX);
  std::string current_working_dir(buff);
  return current_working_dir;
}

int main(int argc, char *argv[]) {
  try {
    cout << "You have entered " << argc << " arguments:"
         << "\n";

    for (int i = 0; i < argc; ++i)
      cout << argv[i] << "\n";

    nanogui::init();

    gg::SimpleAppPtr app = gg::SimpleApp::create(std::string(argv[0]));

    app->setScene(FaceMergeScene::create());
    app->drawAll();
    app->setVisible(true);
    nanogui::mainloop();
    // delete app;
    nanogui::shutdown();

  }

  catch (const std::runtime_error &e) {
    std::string error_msg =
        std::string("Caught a fatal error: ") + std::string(e.what());

#if defined(WIN32)
    MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#else
    std::cerr << error_msg << endl;
#endif

    return -1;
  }

  return 0;
}
