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
#include <random>
#include <string>

#include "GaudiGraphics/buffers.hpp"
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"

#include "manifold/bins.hpp"
#include "manifold/m2Includes.h"
#include "manifold/m2Operators.h"
#include "manifold/make.hpp"
#include "manifold/objloader.hpp"

#include "manifold/harmonic_integrators.hpp"
#include "manifold/tree_code.hpp"

#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
#include <iostream>
//#include "m2Operators.h"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;

using namespace GaudiMath;
// typedef  euclidean_space<double> space;

struct CSC {
  vector<int> irow; // Row index of all nonzero elements of A.
  vector<int> pcol; // Pointer to the beginning of each column (in irow and A).
  vector<double> A; // Nonzero elements of A.
};

template <typename SPACE> void buildSymMatrix(m2::surf<SPACE> &in, CSC &m) {
  M2_TYPEDEFS;

  vector<vertex_ptr> &tverts = in.get_vertices();
  for (int i = 0; i < tverts.size(); i++) {
    vertex_ptr v = tverts[i];
    face_vertex_ptr itb = v->fbegin();
    face_vertex_ptr ite = itb->vprev();
    float n = 0.0;
    bool iterating = true;
    T mii = 0.0;
    int sz = 0;

    typedef std::pair<float, int> ColumnEntry;
    std::vector<ColumnEntry> adjacency;
    m2::mesh_calculator<SPACE> mcalc;
    while (iterating) {
      iterating = itb != ite;
      sz++;
      int j = itb->next()->vertex()->position_in_set();
      // double area = itb->next()->face()->calc_area();
      // double area = mcalc.cotan(itb);
      // double area = mcalc.baryArea(itb);
      double area = 1.0;

      mii += area;
      if (j > i) {
        // ColumnEntry adj(-1.0 + double(rand()%10)/100.0, j);
        ColumnEntry adj(-area, j);

        adjacency.push_back(adj);
      }

      itb = itb->vnext();
    }

    m.pcol.push_back(m.irow.size());
    m.irow.push_back(i);
    m.A.push_back(1.0 * mii);

    std::sort(adjacency.begin(), adjacency.end(),
              [](ColumnEntry a, ColumnEntry b) { return b.second > a.second; });

    for (int j = 0; j < adjacency.size(); j++) {
      m.A.push_back(adjacency[j].first);
      m.irow.push_back(adjacency[j].second);
    }
  }
  m.pcol.push_back(m.irow.size());
  /*
  std::cout << "A: ";
  for(int i = 0; i < m.A.size(); i++){
    std::cout << m.A[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "irow: ";
  for(int i = 0; i < m.irow.size(); i++){
    std::cout << m.irow[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "pcol: ";
  for(int i = 0; i < m.pcol.size(); i++){
    std::cout << m.pcol[i] << " ";;
  }
  std::cout << std::endl;
  std::cout << "A[pcol]: ";
  for(int i = 0; i < m.pcol.size()-1; i++){
    int b = m.pcol[i];
    int e = m.pcol[i+1];
    for(int j = b; j < e; j++)
      std::cout << m.A[j] << " ";;
  }
  std::cout << std::endl;
  */
};

class MainScene;
using MainScenePtr = std::shared_ptr<MainScene>;

class MainScene : public gg::Scene {

public:
  static MainScenePtr create() { return std::make_shared<MainScene>(); }

  MainScene() : gg::Scene() {
    using namespace nanogui;

    _debugLines = gg::DebugBuffer::create();
    _debugLines->init();

    mSceneObjects.push_back(_debugLines);

    createGeometry();
    _debugLines->renderLines();
  }

  space3::vec3 rainbow(double d) {
    double r = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.000));
    double g = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.333));
    double b = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.666));
    return space3::vec3(r, g, b);
  }

  void createGeometry() {
    std::cout << "creating buffer" << std::endl;
    m2::obj_loader<space3> load;
    m2::modify<space3> mod;

    std::cout << "  loading assets" << std::endl;
    _meshGraph = &load("assets/messer.obj");
    _meshGraph->update_all();
    std::cout << "--pack" << std::endl;
    _meshGraph->pack();

    std::cout << "  init buffer" << std::endl;

    mod.centerGeometry(*_meshGraph);
    //_obj = gg::BufferObject::create();
    //_obj->init();
    // gg::fillBuffer(_meshGraph, _obj);
    // mSceneObjects.push_back(_obj);
#if 0
    typedef m2::vertex<space3> *vertex_ptr;
    std::vector<vertex_ptr> tverts = _meshGraph->get_vertices();
    int n = m.pcol.size() - 1;
    int nnz = m.A.size();
    char uplo = 'L';
    int nSolutions = 500;
    int nconv;                               // Number of converged eigenvalues.
    double *EigVal = new double[nSolutions]; // Real part of the eigenvalues.
    double *EigVec = new double[nSolutions * n]; // Eigenvectors.

    int ncv = 0;
    double tol = 0.01;
    int maxit = 0;
    double *resid = NULL;
    bool autoShift = true;
    // Finding the five eigenvalues with largest magnitude
    // and the related eigenvectors.
    std::cout << n << " " << nnz << std::endl;

    ARluSymMatrix<double> matrix(n, nnz, &m.A.front(), &m.irow.front(),
                                 &m.pcol.front(), uplo);
    ARluSymStdEig<double> prob(nSolutions, matrix, "SA", ncv, tol, maxit, resid,
                               autoShift);

    nconv = prob.EigenValVectors(EigVec, EigVal);
    /*
    nconv = AREig(EigVal, EigVec, n, nnz,
      &m.A.front(), &m.irow.front(), &m.pcol.front(),
      uplo, nSolutions, "SM");
    */
    // Printing eigenvalues.

    std::cout << "Eigenvalues (" << nconv << "):" << std::endl;
    for (int i = 0; i < nconv; i++) {
      std::cout << "  lambda[" << (i + 1) << "]: " << EigVal[i];
    }

    for (int i = 0; i < tverts.size(); i++) {
      vertex_ptr v = tverts[i];

      float eigVali = EigVec[499 * n + i];
      // float eigVali = EigVec[i];

      // std::cout << eigVali << " "  << std::endl;
      tverts[i]->coordinate() += 5.0 * eigVali * tverts[i]->normal();
      // tverts[i]->coordinate() += 0.5*tverts[i]->normal();
    }
    std::cout << std::endl;

#endif
  }

  virtual void onAnimate() {

    m2::modify<space3> mod;
    gg::fillBuffer(_meshGraph, _obj);
  }

  virtual void onDraw(gg::Viewer &viewer) {
    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });
  }

private:
  std::vector<gg::DrawablePtr> mSceneObjects;

  gg::BufferObjectPtr _obj;
  gg::PointBufferPtr _points;
  gg::DebugBufferPtr _debugLines;

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

    app->setScene(MainScene::create());
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
