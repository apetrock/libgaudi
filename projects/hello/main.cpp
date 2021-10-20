//#include "nanoguiincludes.h"

#if defined(WIN32)
#include <windows.h>
#endif

#include <stdio.h>  /* defines FILENAME_MAX */
// #define WINDOWS  /* uncomment this line to use it for windows.*/ 
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif
#include<iostream>

//#include <nanogui/glutil.h>
#include <nanogui/nanogui.h>

#include <iostream>
#include <string>


#include "GaudiGraphics/viewer.hpp"
#include "GaudiGraphics/buffers.hpp"

#include "manifold/m2Includes.h"
#include "manifold/objloader.hpp"
#include "manifold/bins.hpp"
#include "manifold/fdt.hpp"
#include "manifold/make.hpp"
#include "manifold/m2Operators.h"

#include "manifold/moving_mesh.hpp"
//#include "m2Operators.h"

#define TRACKBALLSIZE  (0.8f)
#define RENORMCOUNT 97

using std::cout;
using std::cerr;
using std::endl;

//typedef  euclidean_space<double> space;  


template<typename SPACE>
void centerGeometry(m2::control<SPACE>& in){
  M2_TYPEDEFS
    vector<vertex_ptr>& tVerts = in.get_vertices();
  int fc = 0;
  coordinate_type cen = in.calc_center();
  m2::bounding_box<SPACE> bb = in.calc_bbox();
      
  coordinate_type dl = bb.max - bb.min;
  T maxl = dl[0];
  maxl = maxl > dl[1] ? maxl : dl[1]; 
  maxl = maxl > dl[2] ? maxl : dl[2]; 
  T s = 2.0/maxl;
  m2::modify<SPACE> mod;
  std::cout << "offset: "<< cen.transpose() << std::endl;
  std::cout << "scale: "<< s << std::endl;
  mod.translate(&in,-cen[0],-cen[1],-cen[2]);
  mod.scale(&in,s,s,s);
};



template <typename SPACE>
void insertDebugTri(m2::control<SPACE > & mesh, int i, GaudiGraphics::ImmediateLines * buffer){  
  using namespace nanogui;
  typedef Eigen::Vector3f Vec3;
  typedef Eigen::Vector3f Vec4;
  M2_TYPEDEFS;

  std::vector<face_ptr>   faces = mesh.get_faces();
  std::vector<vertex_ptr> verts = mesh.get_vertices();
  face_ptr fi = faces[i];
  coordinate_type c0 = fi->fbegin()->coordinate();
  coordinate_type c1 = fi->fbegin()->next()->coordinate();
  coordinate_type c2 = fi->fbegin()->prev()->coordinate();
  buffer->pushLine(xyz(c0), xyz(c1));
  buffer->pushLine(xyz(c1), xyz(c2));
  buffer->pushLine(xyz(c2), xyz(c0));
};

template <typename SPACE>
void fillDebugLines(m2::control<SPACE > & mesh, GaudiGraphics::ImmediateLines * buffer){  
  using namespace nanogui;
  typedef Eigen::Vector3f Vec3;
  typedef Eigen::Vector3f Vec4;
  M2_TYPEDEFS;

  std::vector<edge_ptr>   edges = mesh.get_edges();
  std::vector<vertex_ptr> verts = mesh.get_vertices();
  for(int i = 0; i < edges.size(); i++){
    edge_ptr e = edges[i];
    coordinate_type c1 = e->v1()->coordinate();
    coordinate_type c2 = e->v2()->coordinate();
    buffer->pushLine(xyz(c1), xyz(c2));
  }  
};

template <typename SPACE>
void fillBuffer(m2::control<SPACE > * mesh, GaudiGraphics::BufferObjectPtr obj){  
  using namespace nanogui;
  M2_TYPEDEFS;
  std::vector<face_ptr>   faces = mesh->get_faces();
  std::vector<vertex_ptr> verts = mesh->get_vertices();
  
  int numVerts = 0;
  int numIndices = 0;

  for(int i = 0; i < faces.size(); i++){
    numIndices += faces[i]->size();
  }
  numVerts = verts.size();

  std::cout << numIndices << " " << numVerts << std::endl;
  obj->fillBuffer([&](GaudiGraphics::BufferObject & o)->void {
          o.allocateVerts(faces.size(), numVerts);
          auto & indices = o.indices();
          auto & positions = o.positions();

          std::cout << " face size: " << faces.size() << std::endl;
          for(int i = 0; i < faces.size(); i++){
            face_vertex_ptr fvb = faces[i]->fbegin();
            face_vertex_ptr fve = faces[i]->fend();
            bool it = true;
            int j = 0;
                        
            while(it){ 
              it = fvb != fve;
              indices.col(i)[j] = fvb->vertex()->position_in_set();
              //std::cout << fvb->vertex()->position_in_set() << " ";
              //                    std::cout << fvb->vertex()->position_in_set() << std::endl;
              j++; fvb = fvb->next();
              if(j==3) break;
            }
                              //std::cout << std::endl;
            //std::cout << indices.col(i) << std::endl;
          }
                
          for(int i = 0; i < verts.size(); i++){
            for(int j = 0; j < 3; j++){
              //std::cout << verts[i]->coordinate()[j] << std::endl;
              positions.col(i)[j] = verts[i]->coordinate()[j];
            }
          }
                //std::cout << indices << std::endl;
         }
      );
};




class HelloApplication : public nanogui::Screen {
public:
  
  typedef double Real;
  
  typedef Eigen::Matrix3f Mat3;
  typedef Eigen::Matrix4f Mat4;
  typedef Eigen::Vector2i Vec2i;
  typedef Eigen::Vector3f Vec3;
  typedef Eigen::Vector4f Vec4;
  typedef Eigen::Vector4d Vec4d;
  typedef Eigen::Quaternion<float> Quat;
  typedef euclidean_space<Real> space;
  typedef m2::BlockTree<space, 8,4,4, std::vector<m2::vertex<space>*> > BlockTree;
  typedef m2::BlockTree<space, 8,4,4, std::vector<m2::face<space>*> > BlockTriTree;
  
  
  HelloApplication(std::string file) : nanogui::Screen(Eigen::Vector2i(800, 800), "NanoGUI Test") {
    using namespace nanogui;
    
    //now for GUI
    Window *window = new Window(this, "coordinates");
    window->setPosition(Vector2i(15, 15));
    window->setLayout(new GroupLayout());

    performLayout(mNVGContext);

    _viewer = GaudiGraphics::Viewer::create(mSize);

    m2::obj_loader<space3>	load;
    m2::subdivide<space3> sub;
    m2::make<space3>	    mk;
    m2::convex_hull<space3>   ch;
    m2::add_handle<space3>    ah;


    m2::construct<space3>	    bevel;
    m2::modify<space3> mod;  

    _meshGraph = &load("assets/bunny.obj");
    //std::cout << "--make cube" << std::endl;
    //_meshGraph =  mk.cube(1.0,1.0,1.0);  
    //_meshGraph =  &sub.subdivide_control(*_meshGraph);
    
    //_meshGraph =  &sub.subdivide_control(*_meshGraph);  
    //_meshGraph =  &sub.subdivide_control(*_meshGraph);
    //_meshGraph =  &sub.subdivide_control(*_meshGraph);

    //_meshGraph->print_stack();
    
    
    std::cout << "--center" << std::endl;
    std::cout << "--update_all" << std::endl;
    _meshGraph->update_all();
    std::cout << "--pack" << std::endl;
    _meshGraph->pack();
    std::cout << "--build sheet" << std::endl;
    _vortex = new m2::vortex_sheet<space>(_meshGraph);

    centerGeometry(*_meshGraph);

    _obj = GaudiGraphics::BufferObject::create();
    _obj->initBuffer();

    std::cout << _obj->scale() << std::endl;

    mSceneObjects.push_back(_obj);
     
    /*
      ImmediateLines * blockLines = new ImmediateLines();
      ImmediateLines * triLines = new ImmediateLines();
      ImmediateLines * stepLines = new ImmediateLines();
    
      blockLines->buildBuffer();
      blockLines->initDisplayShader();
      triLines->buildBuffer();
      triLines->initDisplayShader();
      stepLines->buildBuffer();
      stepLines->initDisplayShader();

      triLines->setColor(Vec3(0,1,0));
      stepLines->setColor(Vec3(0,0,1));
    
      mDebug.push_back(blockLines);
      mDebug.push_back(triLines);
      mDebug.push_back(stepLines);
    */
    // Enable depth test
    glEnable(GL_DEPTH_TEST);
    // Accept fragment if it closer to the camera than the former one
    glDepthFunc(GL_LESS); 
    // Cull triangles which normal is not towards the camera
    glEnable(GL_CULL_FACE);
     

  }
    
  ~HelloApplication() {
  }
    
  ///////////////////////////
  //  events
  ///////////////////////////
  
    
  virtual bool mouseButtonEvent(const Eigen::Vector2i &p, 
				int button, bool down, int modifiers) {
    _viewer->onMouseButton(p, button, down, modifiers);
    Screen::mouseButtonEvent(p, button, down,  modifiers);
    return true;
  }
    
    
    
  virtual bool mouseMotionEvent(const Eigen::Vector2i &p, 
				const Eigen::Vector2i &rel, 
				int button, int modifiers) {

    _viewer->onMouseMotion(p, rel, button, modifiers);
    Screen::mouseMotionEvent(p, rel, button, modifiers);
    return true;
  }



  ///////////////////////////
  //  Draw
  ///////////////////////////

  virtual void draw(NVGcontext *ctx) {
    /* Draw the user interface */
    double max = 0.05;
    double min = 0.01;
    double dt = 0.1;

    _vortex->minLength = max;
    _vortex->minCollapseLength = min;
    _vortex->integrateBaroclinity(dt);
    //_vortex->addCurveVorticity(0.05, 1);
    _vortex->updateCirculation();
    _vortex->integrateVelocityRK2(dt);
    _vortex->remesh();
    _vortex->updateVorticity();
    m2::remesh<space3>	    rem;
    rem.triangulate(_meshGraph);

    _meshGraph->update_all();
    _meshGraph->pack();
    centerGeometry(*_meshGraph);
    fillBuffer(_meshGraph, _obj);
    if(false)
      Screen::draw(ctx);
  }

  virtual void drawContents() {
    using namespace nanogui;
    glfwGetTime();
    std::for_each(mSceneObjects.begin(), mSceneObjects.end(), 
		  [&](GaudiGraphics::BufferObjectPtr obj) mutable {
		    obj->onFrame(_frame, *obj, *_viewer, nullptr );
		    if(obj->isVisible)
		      obj->draw(_viewer->getProjection(), _viewer->getModelView());
		  });
    /*
      std::for_each(mDebug.begin(), mDebug.end(), 
      [&](BufferObject * obj) mutable {
      //std::cout << obj->isVisible << std::endl;
      if(obj->isVisible)
      obj->draw(mProject, mModelView);
      });
    */
    _frame++;
  }


private:
  unsigned int _frame = 0;
  std::vector<GaudiGraphics::BufferObjectPtr> mDebug;
  std::vector<GaudiGraphics::BufferObjectPtr> mSceneObjects;
  GaudiGraphics::ViewerPtr _viewer;

  GaudiGraphics::BufferObjectPtr _obj; 
  m2::control<space >	* _meshGraph;
  m2::vortex_sheet<space > * _vortex;
};

std::string GetCurrentWorkingDir( void ) {
  char buff[FILENAME_MAX];
  GetCurrentDir( buff, FILENAME_MAX );
  std::string current_working_dir(buff);
  return current_working_dir;
}

int main(int argc, char *argv[]) {
  try {
    cout << "You have entered " << argc 
         << " arguments:" << "\n"; 
  
    for (int i = 0; i < argc; ++i) 
        cout << argv[i] << "\n"; 
  
    nanogui::init();
    HelloApplication *app = new HelloApplication(std::string(argv[0]));
    app->drawAll();
    app->setVisible(true);
    nanogui::mainloop();
    delete app;
    nanogui::shutdown();
  }

  catch (const std::runtime_error &e) {
    std::string error_msg = std::string("Caught a fatal error: ")
      + std::string(e.what());

#if defined(WIN32)
    MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#else
    std::cerr << error_msg << endl;
#endif

    return -1;
  }

  return 0;
}
