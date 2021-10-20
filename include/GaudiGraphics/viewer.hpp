#ifndef __APP_VIEWER__
#define __APP_VIEWER__

#include <nanogui/glutil.h>
#include <GaudiMath/typedefs.h>

#include <iostream>
#include <string>
#include <memory>

namespace GaudiGraphics{
  using namespace GaudiMath;
  class Viewer;
  using ViewerPtr = std::shared_ptr<Viewer>;

  class Viewer {
  public:
    static ViewerPtr  create(Vec2i sz){
      return
	std::make_shared<Viewer>(sz);
    }

    Viewer(){
      //init();
    }
    
    Viewer(Vec2i sz){
      init(sz);
    }
    
    ~Viewer(){}


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

    
    void init(Vec2i size){
      
      mProject = this->makeProjectionMatrix(60.0*M_PI/180.0, 1.0, 0.1, 20);
      mModelView.setIdentity();
      mModelViewOld.setIdentity();
    
      mModelRotOld.setIdentity();
      mModelRotNew.setIdentity();

      mDragging = false;
    
      ball = new nanogui::Arcball();
      ball->setSize(size);
      mPosition = Vec4(0,0,4,0);
      
    }

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

    virtual bool onMouseButton(const Eigen::Vector2i &p, 
			       int button, bool down, int modifiers) {
      updateFrame(); 
      ball->button(p,down);
      mDragging = false;
    	  
      return true;
    }
    
    
    virtual bool onMouseMotion(const Eigen::Vector2i &p, 
			       const Eigen::Vector2i &rel, 
			       int button, int modifiers) {

      ball->motion(p);
      updatePosition();
      return true;
    }

    Mat4 & getModelView(){ return mModelView;}
  Mat4 & getProjection(){ return mProject;}

  protected:
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

    nanogui::Arcball * ball;

  };
}
#endif
