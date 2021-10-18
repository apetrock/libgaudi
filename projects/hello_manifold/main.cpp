#include <nanogui/glutil.h>
#include <nanogui/nanogui.h>

#include <iostream>
#include <string>

//#include "meshHandler.hpp"
//#include "movingMesh.hpp"
//#include "make.hpp"


#include "buffers.hpp"

#include "manifold/m2Includes.h"
#include "manifold/objloader.hpp"
#include "manifold/movingMesh.hpp"
#include "manifold/bins.hpp"
#include "manifold/fdt.hpp"
#include "manifold/make.hpp"
#include "manifold/m2Operators.h"


#include <cmath>

//#include "m2Operators.h"

#define TRACKBALLSIZE  (0.8f)
#define RENORMCOUNT 97

using std::cout;
using std::cerr;
using std::endl;

//typedef  euclidean_space<double> space;  

struct CSC {
  vector<int>    irow;  // Row index of all nonzero elements of A.
  vector<int>    pcol;  // Pointer to the beginning of each column (in irow and A).
  vector<double> A;     // Nonzero elements of A.

};

template<typename SPACE>
void buildSymMatrix(m2::control<SPACE>& in, CSC& m){
  M2_TYPEDEFS;
  
  vector<vertex_ptr>& tverts = in.get_vertices();
  for(int i = 0; i < tverts.size(); i++){
    vertex_ptr v = tverts[i];
    face_vertex_ptr itb = v->fbegin();
    face_vertex_ptr ite = itb->vprev(); 
    float n = 0.0;
    bool iterating = true;
    T   mii = 0.0;
    int sz = 0;
    
    
    typedef std::pair<float, int> ColumnEntry;
    std::vector<ColumnEntry> adjacency;
    m2::mesh_calculator<SPACE> mcalc;
    while (iterating) {
      iterating = itb != ite;
      sz++;
      int j = itb->next()->vertex()->position_in_set();
      //double area = itb->next()->face()->calc_area();
      //double area = mcalc.cotan(itb);
      // double area = mcalc.baryArea(itb);
     double area = 1.0;
      
      mii+=area; 
      if(j > i){
        //ColumnEntry adj(-1.0 + double(rand()%10)/100.0, j);
        ColumnEntry adj(-area, j);
        
        adjacency.push_back(adj);
      }

      itb = itb->vnext();
    }

    m.pcol.push_back(m.irow.size());
    m.irow.push_back(i);
    m.A.push_back(1.0*mii);

    std::sort(adjacency.begin(), adjacency.end(), 
              [](ColumnEntry a, ColumnEntry b) {
                return b.second > a.second;   
              });

    for(int j=0; j<adjacency.size();j++){
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

template<typename SPACE>
void centerGeometry(m2::control<SPACE>& in){
  M2_TYPEDEFS;
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
void insertDebugTri(m2::control<SPACE > & mesh, int i, ImmediateLines * buffer){  
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
void fillDebugLines(m2::control<SPACE > & mesh, ImmediateLines * buffer){  
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
BufferObject * fillBuffer(m2::control<SPACE > & mesh){  
  using namespace nanogui;
  M2_TYPEDEFS;
  std::vector<face_ptr>   faces = mesh.get_faces();
  std::vector<vertex_ptr> verts = mesh.get_vertices();
  
  int numVerts = 0;
  int numIndices = 0;

  for(int i = 0; i < faces.size(); i++){
      numIndices += faces[i]->size();
  }
  numVerts = verts.size();
  
  BufferObject * cube = new BufferObject();
  std::cout << numIndices << " " << numVerts << std::endl;
  cube->buildBuffer(faces.size(), numVerts,
		   [&](nanogui::MatrixXu & indices, 
		      nanogui::MatrixXf & positions)->void {

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
                          //std::cout << fvb->vertex()->position_in_set() << std::endl;
			  j++; fvb = fvb->next();
			}
                        //std::cout << std::endl;
			//std::cout << indices.col(i) << std::endl;
		      }
		      
		      for(int i = 0; i < verts.size(); i++){
			for(int j = 0; j < 3; j++){
			  positions.col(i)[j] = verts[i]->coordinate()[j];
			}
		      }
		      //std::cout << indices << std::endl;
		   });
  cube->initDisplayShader();
  //cube->initPickingShader();
  
  return cube;
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

  HelloApplication() : nanogui::Screen(Eigen::Vector2i(1600, 1600), "NanoGUI Test") {
    
    using namespace nanogui;
    
    //now for GUI
    Window *window = new Window(this, "coordinates");
    window->setPosition(Vector2i(15, 15));
    window->setLayout(new GroupLayout());
    //performLayout(mNVGContext);

      
    textBoxX = new TextBox(window);
    textBoxY = new TextBox(window);
    textBoxZ = new TextBox(window);
    textBoxX->setFixedSize(Vector2i(120, 25));
    textBoxY->setFixedSize(Vector2i(120, 25));
    textBoxZ->setFixedSize(Vector2i(120, 25));
      
    textBoxX->setValue("");
    textBoxY->setValue("");
    textBoxZ->setValue("");
      
    performLayout(mNVGContext);

    mProject = this->makeProjectionMatrix(60.0*M_PI/180.0, 1.0, 0.1, 20);
    mModelView.setIdentity();
    mModelViewOld.setIdentity();
    
    mModelRotOld.setIdentity();
    mModelRotNew.setIdentity();

    mDragging = false;
    mCurrentHover = -1;
    mCurrentGroup = 0;
    mCurrentSelect.push_back(-1);
    mCurrentSelect.push_back(-1);
 
    ball = new nanogui::Arcball();
    ball->setSize(mSize);

    mPosition = Vec4(0,0,4,0);
    initObjects();
      
      
    m2::obj_loader<space3>	load;
    m2::subdivide<space3> sub;
    m2::make<space3>	    mk;
    m2::convex_hull<space3>   ch;
    m2::add_handle<space3>    ah;
    m2::remesh<space3>	    rem;

    m2::construct<space3>	    bevel;
    m2::modify<space3> mod;  
    
    _meshGraph = &load("assets/bunny.obj");
    //_meshGraph = &load("assets/models/donut.obj");
    //_meshGraph =  mk.cube(1.0,1.0,1.0);  
    //_meshGraph =  &sub.subdivide_control(*_meshGraph);
    //_meshGraph =  &sub.subdivide_control(*_meshGraph);  
    //_meshGraph =  &sub.subdivide_control(*_meshGraph);
    //_meshGraph =  &sub.subdivide_control(*_meshGraph);
    //_meshGraph =  &sub.subdivide_control(*_meshGraph);
    //_meshGraph =  &sub.subdivide_control(*_meshGraph);
    
    //_meshGraph->print_stack();
    rem.triangulate(_meshGraph);
    

    centerGeometry(*_meshGraph);
    _meshGraph->update_all();
    _meshGraph->pack();
    

    CSC arp_matrix;
    CSC & m = arp_matrix;
    
    buildSymMatrix(*_meshGraph, arp_matrix);
    
    BufferObject * obj = fillBuffer(*_meshGraph);
    std::cout << obj->scale() << std::endl;
    mSceneObjects.push_back(obj);
     
    ImmediateLines * blockLines = new ImmediateLines();
    ImmediateLines * triLines = new ImmediateLines();
    ImmediateLines * stepLines = new ImmediateLines();
    
    /*
    blocks = makeBlock(_meshGraph);
    std::cout << "bsize: "<< blocks->zRes() << " " 
	      << blocks->zRes() << " " 
	      << blocks->zRes() << std::endl; 
    for(int z = 0; z < blocks->zRes(); z++){
      for(int y = 0; y < blocks->yRes(); y++){
	for(int x = 0; x < blocks->xRes(); x++){
	  if(blocks->mask(x,y,z)){
	    
	    Vec4 cen = blocks->center(x,y,z).cast<float>();
	    double scale = blocks->dx();
	    Vec3 scale3 = Vec3(0.5*scale, 0.5*scale, 0.5*scale);
	    blockLines->pushBox(xyz(cen), scale3);
	    //BufferObject * cell = makeCube();
	      
	  }
	  }
	  }
	  }
	  blockLines->buildBuffer();
	  blockLines->initDisplayShader();
	  mDebug.push_back(blockLines);
    */

    //tree = m2::makeTree<space,8,8,8>(_meshGraph);
    
    //fillDebugLines(*_meshGraph, blockLines);  

#if 0
    tritree = m2::makeTriTree<space,8,4,4>(_meshGraph);
    //tritree = m2::makeTriTree<space,8,4,4>(_meshGraph2, stepLines, triLines);
    
    BlockTriTree & root = *tritree;
    for(int z0 = 0; z0 < root.zRes(); z0++){
      for(int y0 = 0; y0 < root.yRes(); y0++){
	for(int x0 = 0; x0 < root.xRes(); x0++){
	  if(!root.mask(x0,y0,z0)) continue;
          std::cout << " " << x0 << " " << y0 << " " << z0 << std::endl;
	  Vec4 cen0 = root.center(x0,y0,z0).cast<float>();
	  double scale0 = root.dx();
	  Vec3 scale03 = Vec3(0.5*scale0, 0.5*scale0, 0.5*scale0);
	  blockLines->pushBox(xyz(cen0), scale03);
	  BlockTriTree::Block1 * b1 = tritree->get(x0,y0,z0);
	  for(int z1 = 0; z1 < b1->zRes(); z1++){
	    for(int y1 = 0; y1 < b1->yRes(); y1++){
	      for(int x1 = 0; x1 < b1->xRes(); x1++){
		if(!b1->mask(x1,y1,z1)) continue;
		Vec4 cen1 = b1->center(x1,y1,z1).cast<float>();
		double scale1 = b1->dx();
		Vec3 scale13 = Vec3(0.5*scale1, 0.5*scale1, 0.5*scale1);
		blockLines->pushBox(xyz(cen1), scale13);
		BlockTriTree::Block2 * b2 = b1->get(x1,y1,z1);	   
		for(int z2 = 0; z2 < b2->zRes(); z2++){
		  for(int y2 = 0; y2 < b2->yRes(); y2++){
		    for(int x2 = 0; x2 < b2->xRes(); x2++){
		      if(!b2->mask(x2,y2,z2)) continue;
		      Vec4 cen2 = b2->center(x2,y2,z2).cast<float>();
		      double scale2 = b2->dx();
		      Vec3 scale23 = Vec3(0.5*scale2, 0.5*scale2, 0.5*scale2);
		      blockLines->pushBox(xyz(cen2), scale23);
		      //BlockTree::Block1 & b1 = tree->get(x1,y1,z1);	   
		
		    }}}
	      }}}
	  
	}}}
    #endif

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
 
    // Enable depth test
    glEnable(GL_DEPTH_TEST);
    // Accept fragment if it closer to the camera than the former one
    glDepthFunc(GL_LESS); 
    // Cull triangles which normal is not towards the camera
    glEnable(GL_CULL_FACE);
     

  }
    
    ~HelloApplication() {
    }
    

    void updateTextBox(Vec4 pos){
      textBoxX->setValue(realToString(pos[0]));
      textBoxY->setValue(realToString(pos[1]));
      textBoxZ->setValue(realToString(pos[2]));
      
      textBoxX->setUnits("x");
      textBoxY->setUnits("y");
      textBoxZ->setUnits("z");
    
    }
    
    void initCube(BufferObject & cube){
      //init a cube with random size, and shapes
      //next use rotations, but haven't tested those
      auto fRand = [](double fMin, double fMax)->double{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
      };

      cube.setColor(Vec3(fRand(.1,.9), fRand(.1,.9), fRand(.1,.9)));
      cube.translate(Vec4(fRand(-1.6,1.6), 
			  fRand(-1.6,1.6), 
			  fRand(-1.6,1.6),1));
      cube.setScale(Vec4(fRand(.1,.25), fRand(.1,.25), fRand(.1,.25), 1.0));
      
    }

    void initObjects(){      
      
      for(int i=0; i < 0; i++){
	BufferObject * cube = makeCube();	
	initCube(*cube);
	mSceneObjects.push_back(cube);
      }
    }

    virtual int pickObject(const Eigen::Vector2i &p, Vec4 & hitPoint){
      int isx = -1; Real tmin = 999;
      Vec4 rayDir = castRay(p);
      //mSceneObjects should get abstracted away to be a container
      //that holds objects and can perform scene hit testing.  This way
      //you can add an acceleration structure to speed up testing, and
      //each object can be tested by more than just a broadphase test
      //...if I have time
      
      for(int i = 0; i < mSceneObjects.size(); i++){
	Real tn, tf;
	Vec4 hn, hf;
	if(mSceneObjects.at(i)->intersectBbox(mPosition, rayDir, tn, tf, hn, hf)){
	  if(tn < tmin){
	    hitPoint = hn;
	    tmin = tn;
	    isx = i;
	  }
	}
      }
      
      //if(isx > -1){
      //	std::cout << tmin << std::endl;
      //}
      return isx;
    }

    ///////////////////////////
    //  events
    ///////////////////////////
  
    
    virtual bool mouseButtonEvent(const Eigen::Vector2i &p, 
				  int button, bool down, int modifiers) {
      ball->button(p,down);
      
      Vec4 hit;
      int oldSelection = mCurrentSelect[mCurrentGroup];
      int newSelection = pickObject(p, hit);
      int curSelection;
      
      Screen::mouseButtonEvent(p, button, down,  modifiers);
      return true;
    }
    
    
    
    virtual bool mouseMotionEvent(const Eigen::Vector2i &p, 
          const Eigen::Vector2i &rel, 
          int button, int modifiers) {
      ball->motion(p);
      
      Screen::mouseMotionEvent(p, rel, button, modifiers);
      return true;
    }

    ///////////////////////////
    //  Draw
    ///////////////////////////


    void updateFrame(){
      //cache the old rotation
      mModelRotOld = mModelRotNew;
      ball->setState(Quat::Identity()); //reset ball to identity
    }

    void updatePosition(){
      //This is my camera model, its not a seperat object, but doesn't necessarily
      //warrant one.  All the fancy stuff was implemented in the arcball class
      Mat4 r = ball->matrix();
      //accumulate rotations... this likely could be done with the quaternions in the 
      //arcball, but this works after some fiddling
      mModelRotNew = r*mModelRotOld;
      //set the modelView equal to the new rotation
      mModelView = mModelRotNew;
      mModelView.col(3) << Vec4(0,0,-3,1); 
      //rotate the world, then take a step back
    
      //get the position from the inverse  of the matrix
      mPosition = mModelRotNew.transpose()*Vec4(0,0,4,1);
    }

    /*
      void frameBufferPick(const Eigen::Vector2i &p){
      //color picking seems like the simplest way to go, since Ray tracing can suffer
      //from precision errors, and color picking is quick to implement.  If I have time,
      //I'll try both
      mPickerBuffer.bind();
    
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glEnable(GL_DEPTH_TEST);

      cube.pickingShader().bind();
    
      Mat4 matrix = cube.matrix();
      Mat4 mvp = mProject*mModelView*matrix;
    
      cube.pickingShader().setUniform("MVP", mvp);
      cube.pickingShader().drawIndexed(GL_TRIANGLES, 0, 12);
    
      //mPickerBuffer.blit();
    
      std::vector<std::uint8_t> data(4*4*4);
      glReadBuffer(GL_BACK);
      glReadPixels(p[0], p[1], 4, 4, GL_RGBA8, GL_UNSIGNED_BYTE, &data[0]);
      for(int i = 0; i < data.size(); i++){
      std::cout << data[i]/255 << std::endl;

      }
    
      mPickerBuffer.release();
      };
    */

    virtual void draw(NVGcontext *ctx) {
      /* Draw the user interface */
      if(false)
	Screen::draw(ctx);
    }

    virtual void drawContents() {
      using namespace nanogui;
      glfwGetTime();

      std::for_each(mSceneObjects.begin(), mSceneObjects.end(), 
		    [&](BufferObject * obj) mutable {
		      //std::cout << obj->isVisible << std::endl;
		      if(obj->isVisible)
			obj->draw(mProject, mModelView);
		    });

      std::for_each(mDebug.begin(), mDebug.end(), 
		    [&](BufferObject * obj) mutable {
		      //std::cout << obj->isVisible << std::endl;
		      if(obj->isVisible)
			obj->draw(mProject, mModelView);
		    });
    
    }


  
    Vec3 unproject(const Vec3 &win,
		   const Mat4 &model,
		   const Mat4 &proj,
		   const Vec2i &viewportSize) {
      //the nanogui::unproject appeared to be backwards in the z-direction
      //so I copied it over and changed it here:
      Mat4 Inverse = (proj * model).inverse();

      Vec4 tmp;
      tmp << win, 1;
      tmp(0) = tmp(0) / viewportSize.x();
      tmp(1) = tmp(1) / viewportSize.y();
      //tmp = tmp.array() * 2.0f - 1.0f;
      tmp(0) = tmp(0) * 2.0f - 1.0f;
      tmp(1) = 1.0f - tmp(1) * 2.0f;
    
      Vec4 obj = Inverse * tmp;
      obj /= obj(3);

      return obj.head(3);
    }
  
    Vec4 castRay(const Eigen::Vector2i &p){
      
      Real width  = this->mSize[0];
      Real height = this->mSize[1];
      /* ---this is just to check my position.  If I'm doing it right,
	 ---the back unprojected vector should be very close to my 
	 ---position
	 
	Vec4 iwPosition;
	iwPosition << unproject(Vec3(p[0], p[1], -1), 
	mModelView, mProject, mSize), 1;
      
	//just make sure that 
	//assert((mPosition - iwPosition).norm() < .25);
      */
      Vec4 wPosition(0,0,0,1);
      wPosition << unproject(Vec3(p[0], p[1], 1), 
			     mModelView, mProject, mSize), 1;
     
      Vec4 ray = wPosition - mPosition;
      ray.normalize();
      //rayWorld.normalize();
      //std::cout << rayWorld.transpose() << "\n" << mPosition.transpose() << std::endl;
      return ray;
    
    }
  
    Mat4 makeProjectionMatrix(Real fovY, Real aspectRatio, 
			      Real zNear, Real zFar){
      Real yScale = tan(M_PI_2 - fovY/2.0);
      Real xScale = yScale/aspectRatio;
      Mat4 M;
      M.col(0) << xScale, 0, 0, 0;
      M.col(1) << 0, yScale, 0, 0;
      M.col(2) << 0, 0, -(zFar+zNear)/(zFar-zNear), -1;
      M.col(3) << 0, 0, -2*zNear*zFar/(zFar-zNear), 0;
      return M;
      //return nanogui::frustum(-.1,.1,-.1,.1,.1,20);
    }

    std::vector<BufferObject*> & getHelpers(){
      return mHelpers;
    }

  private:
    std::vector<BufferObject*> mHelpers;
    std::vector<BufferObject*> mDebug;
    std::vector<BufferObject*> mSceneObjects;

  m2::control<space >	* _meshGraph;
    
  m2::Block<space, 24, std::vector<m2::vertex<space>*> >* blocks;
  BlockTree* tree;
  BlockTriTree* tritree;
  //nanogui::GLFramebuffer mPickerBuffer;
     

    //camera variables, this could be in its own class
    Vec4 mPosition;
    Mat4 mProject;
    Mat4 mModelView;
    Mat4 mModelRotNew;
    Mat4 mModelRotOld;
    Mat4 mModelViewOld;
    
    //variables for selection and dragging, I maintain two selection groups
    //one for widgets and one for objects
    bool mDragging;
    Vec4 mObjCenCache;
    Vec4 mDragStart;
    int mCurrentHover; //we'll maintain one hover group
    int mCurrentGroup; //maintain the currently selected group
    std::vector<int> mCurrentSelect; //a vector showing which object is in each group
    
    //GUI
    nanogui::TextBox *textBoxX;
    nanogui::TextBox *textBoxY;
    nanogui::TextBox *textBoxZ;
      

    nanogui::Arcball * ball;

  };

  int main(int argc, char *argv[]) {

    // Printing solution.
    std::cout << argc << " - "<< argv[0] << std::endl;
    try {
      nanogui::init();
      HelloApplication *app = new HelloApplication();
      app->drawAll();
      app->setVisible(true);

      nanogui::mainloop();

      delete app;

      nanogui::shutdown();
    } catch (const std::runtime_error &e) {
      std::string error_msg = std::string("Caught a fatal error: ") + 
        std::string(e.what());

#if defined(WIN32)
      MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#else
      std::cerr << error_msg << endl;
#endif
      return -1;
    }

    return 0;
  }
