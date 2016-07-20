#include <nanogui/screen.h>
#include <nanogui/window.h>
#include <nanogui/layout.h>
#include <nanogui/label.h>
#include <nanogui/checkbox.h>
#include <nanogui/button.h>
#include <nanogui/toolbutton.h>
#include <nanogui/popupbutton.h>
#include <nanogui/combobox.h>
#include <nanogui/progressbar.h>
#include <nanogui/entypo.h>
#include <nanogui/messagedialog.h>
#include <nanogui/textbox.h>
#include <nanogui/slider.h>
#include <nanogui/imagepanel.h>
#include <nanogui/imageview.h>
#include <nanogui/vscrollpanel.h>
#include <nanogui/colorwheel.h>
//#include <trackball.h>
#if defined(WIN32)
#include <windows.h>
#endif
#include <nanogui/glutil.h>
#include <iostream>
#include <string>


#define TRACKBALLSIZE  (0.8f)
#define RENORMCOUNT 97

using std::cout;
using std::cerr;
using std::endl;


//nano gui deals mostly with vec3s, which is annoying
//since I'm adhering to using homogeneous coordinates
//so far as I can tell, Eigen can grab blocks, but
//that syntax is laborious, so here as a few helper functions
Eigen::Vector3f xyz(Eigen::Vector4f in){
  return Eigen::Vector3f(in[0],in[1],in[2]);
};

Eigen::Vector4f xyzw(Eigen::Vector3f in){
  return Eigen::Vector4f(in[0],in[1],in[2], 1);
};

Eigen::Vector4f normalize(Eigen::Vector4f in){
  Eigen::Vector4f out = in;
  out.normalize();
  return out;
};

Eigen::Vector4f hnormalize(Eigen::Vector4f in){
  //this function never gets used.  As it turns out
  //most vectors that need to get normalized are already 
  //nonhomogenous
  in[3] = 0;
  in.normalize();
  in[3] = 1;
  return in;
};

Eigen::Vector4f rayPlaneIntersect(Eigen::Vector4f l0, 
				  Eigen::Vector4f r, 
				  Eigen::Vector4f N, 
				  Eigen::Vector4f p){
  //projects a pt onto a plane with a  normal, N and point c
  //pC[3] = 0;
  // pt[3] = 0;
  N[3] = 0; N.normalize();
  Eigen::Vector4f dpc = p - l0;
  hnormalize(N); //makeSure its Normalized
  Eigen::Vector4f itx = N.dot(dpc)/(N.dot(r))*r; //project vector onto the normal
  return l0 + itx; //return the addition of pt + projected vector
}

Eigen::Vector4f projectToPlane(Eigen::Vector4f pt, 
			       Eigen::Vector4f N, 
			       Eigen::Vector4f c){
  N[3] = 0; N.normalize();
  Eigen::Vector4f dpc = c - pt;
  hnormalize(N); //makeSure its Normalized
  Eigen::Vector4f projpc = N.dot(dpc)*N; //project vector onto the normal
  return pt + projpc; //return the addition of pt + projected vector
}
 
Eigen::Vector4f projectToLine(Eigen::Vector4f p, 
			      Eigen::Vector4f l0, 
			      Eigen::Vector4f l1){
  //very ver similar to project to plane...
 
  Eigen::Vector4f dl = l0 - l1;
  Eigen::Vector4f dlp = l0 - p;
  
  Eigen::Vector4f projdlpdl = dl.dot(dlp)*normalize(dl); //project vector onto the normal
  return l0 + projdlpdl; //return the addition of pt + projected vector
}

Eigen::Matrix4f RandomRotation(){
  auto fRand = [](double fMin, double fMax)->double{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
  };
   Eigen::Quaternion<float> Q = 
     Eigen::Quaternion<float> (Eigen::AngleAxisf(fRand(0, 2.0)*M_PI, 
						 Eigen::Vector3f(fRand(-1,1),
								 fRand(-1,1),
								 fRand(-1,1))));
  Q.normalize();
  Eigen::Matrix3f R3 = Q.matrix();
  Eigen::Matrix4f R4; R4.setIdentity(); R4.block<3,3>(0,0) = R3;
  return R4;
}


std::string realToString(float f){
  std::ostringstream convert;   // stream used for the conversion
  convert << f;      // insert the textual representation of 'Number' in the 
  return convert.str(); // set 'Result' to the contents of the stream
}

//forward declare the application so I can use it as a pointer later,
//without having to use a static cast
class GnomonApplication;

class Drawable {
  //Really simple hierarchy, all classes derive from drawable, and thats really
  //the only method they have, so that we can collect all drawables 
public: 
  typedef double Real;
  
  typedef Eigen::Matrix3f Mat3;
  typedef Eigen::Matrix4f Mat4;
  typedef Eigen::Vector2i Vec2i;
  typedef Eigen::Vector3f Vec3;
  typedef Eigen::Vector4f Vec4;
  typedef Eigen::Quaternion<float> Quat;
  typedef nanogui::MatrixXu MatrixXu;
  typedef nanogui::MatrixXf MatrixXf;
  
  Drawable(){}
  ~Drawable(){}
  virtual void draw(){}
  virtual bool intersectBbox(Vec4 r0, Vec4 r1, Real & tnear, Real & tfar ){};
  
  virtual bool rayBoxIntersect(Vec4 r_o, Vec4 r_d, Vec4 boxmin, Vec4 boxmax,
			       Real & tnear, Real & tfar){
    
    Vec4 invR = Vec4(1.0/r_d[0],1.0/r_d[1],1.0/r_d[2],0.0);
    
    Vec4 tbot, ttop;
    for(int i = 0; i < 4; i++){
      tbot[i] = invR[i] * (boxmin[i] - r_o[i]);
      ttop[i] = invR[i] * (boxmax[i] - r_o[i]);
    }
    // re-order intersections to find smallest and largest on each axis
    Vec4 tmin, tmax;
    for(int i = 0; i < 4; i++){
      tmin[i] =  std::min(ttop[i], tbot[i]);
      tmax[i] =  std::max(ttop[i], tbot[i]);
    }
    
    // find the largest tmin and the smallest tmax
    Real largest_tmin  
      = std::max(std::max(tmin[0], tmin[1]), std::max(tmin[0], tmin[2]));
    Real smallest_tmax 
      = std::min(std::min(tmax[0], tmax[1]), std::min(tmax[0], tmax[2]));

    tnear = largest_tmin;
    tfar = smallest_tmax;
    
    return(smallest_tmax > largest_tmin);
  }
};

class BufferObject : public Drawable{
public:
  BufferObject(){};
  ~BufferObject(){
    //free everything here
    
    mDispShader.free();
    mPickShader.free();
  };
  
    void buildBuffer(int nInd, int nVerts, 
		     std::function<void(nanogui::MatrixXu&, 
					nanogui::MatrixXf&)> buildFunc){
      mIndices   = nanogui::MatrixXu(3, nInd);
      mPositions = nanogui::MatrixXf(3, nVerts);
      buildFunc(mIndices, mPositions);
      computeNormals();
      calcBbox();
      mMatrix.setIdentity();
      mRot.setIdentity();
      mColor = Vec3(0.75, 0.75, 1.0);
      mCen = Vec4(0.0,0.0,0.0,1.0);//add the homogeneous 1
      mScale = Vec4(1.0,1.0,1.0,1.0);
      mOffset = Vec4(0.0,0.0,0.0,1.0);
      selectionGroup = 0;
      isVisible = true;
    };
  
  virtual bool intersectBbox(Vec4 r0, Vec4 r1, 
			     Real & tnear, Real & tfar, 
			     Vec4 & hNear, Vec4 & hFar){
      //the model may have a local transform, so rotate the rays
      //so that they are in model space
    
      Mat4 iMatrix = mMatrix.inverse();
      Vec4 r0p = iMatrix*r0; //p for prime
      Vec4 r1p = iMatrix*r1;
      
      Real tpnear, tpfar;
      bool hit = rayBoxIntersect(r0p,r1p, bbmin, bbmax, tpnear, tpfar);
      Vec4 hpNear = (r0p + tpnear*r1p); 
      Vec4 hpFar  = (r0p + tpfar*r1p);  
                                        
      
      hNear = mMatrix*hpNear;
      hFar = mMatrix*hpFar;
      //we need to reclaim our ts in the new coordinates
      tnear = (hNear - r0).norm()/r1.norm(); //these should be the same
      tfar = (hFar - r0).norm()/r1.norm();   //as the old ones
#if 0
	std::cout << " ---- " << std::endl;
      
	std::cout << iMatrix << std::endl;
	std::cout << iMatrix.inverse() << std::endl;
	std::cout << mMatrix << std::endl;
	
	std::cout << "hnp: " << hNear.transpose() << std::endl;
	std::cout << "hfp: " << hFar.transpose() << std::endl;
	std::cout << "hn: " << (r0 + tnear*r1).transpose() << std::endl;
	std::cout << "hf: " << (r0 + tfar*r1).transpose() << std::endl;

	std::cout << "tpn: " << tpnear << std::endl;
	std::cout << "tpf: " << tpfar << std::endl;
	std::cout << "tn: " << tnear << std::endl;
	std::cout << "tf: " << tfar << std::endl;
#endif
      return hit;
      
    }
  
    void calcBbox(){
      mVertNormals = MatrixXf(mPositions.rows(), mPositions.cols());
      mVertNormals.setZero();
      //in order to calculate the vertex normals to send to the shader
      //the face normals have to be calculated
      bbmin = Vec4(99999,99999,99999,1);
      bbmax = Vec4(-99999,-99999,-99999,1);
      for(int i = 0; i < mPositions.cols(); i++){
	Vec3 v = mPositions.col(i);
	for(int j = 0; j < 3; j++){
	  bbmin[j] = std::min(bbmin[j], v[j]);
	  bbmax[j] = std::max(bbmax[j], v[j]);
      
	}
      }
    }

    void computeNormals(){
      //this isn't really that necessary, now, normals are computed
      //per fragment on the shader
      mVertNormals = MatrixXf(mPositions.rows(), mPositions.cols());
      mVertNormals.setZero();
      //in order to calculate the vertex normals to send to the shader
      //the face normals have to be calculated
      for(int i = 0; i < mIndices.cols(); i++){
	Vec3 v0 = mPositions.col(mIndices(0,i));
	Vec3 v1 = mPositions.col(mIndices(1,i));
	Vec3 v2 = mPositions.col(mIndices(2,i));
      
	Vec3 N = (v1-v0).cross(v2-v0);
	N.normalize();
	mVertNormals.col(mIndices(0,i)) += N;
	mVertNormals.col(mIndices(1,i)) += N;
	mVertNormals.col(mIndices(2,i)) += N;
      }
    
      for(int i = 0; i < mVertNormals.cols(); i++){
	mVertNormals.col(i).normalize();
      }
    
    }
  
    void initDisplayShader(){
      //std::cout << " system: ";
      //system("less src/shaders/standard.vs");
      mId = rand();
      int Number = mId;       // number to be converted to a string
      std::string Result;          // string which will contain the result
      std::ostringstream convert;   // stream used for the conversion
      convert << "a_standard_shader: " << Number;      // insert the textual representation of 'Number' in the characters in the stream
      Result = convert.str();
      mDispShader.initFromFiles(
				/* An identifying name */
				convert.str(),
				"src/shaders/standard.vs",
				"src/shaders/standard.fs");

      mDispShader.bind();
      mDispShader.uploadIndices(mIndices);
      mDispShader.uploadAttrib("vertexPosition_modelspace", mPositions);
      mDispShader.setUniform("uColor", mColor);
    };

    void initPickingShader(){
      //I couldn't get backbuffer style picking to work in short time, so I abandoned it,
      //issues with opengl set or nanoguis wrapper wasn't working, though I tried
      //using raw openGL as well. I imagine that its my fault somewhere, but I didn't
      //want to spend too much longer debugging it
      //here is the legacy code:
      mPickShader.init(
		       /* An identifying name */
		       "a_picking_shader",

		       /* Vertex shader */
		       "#version 330\n"
		       "uniform mat4 MVP;\n"
		       "in vec3 vertexPosition_modelspace;\n"
		       "void main() {\n"
		       "    gl_Position = MVP * vec4(vertexPosition_modelspace, 1.0);\n"
		       "}",

		       /* Fragment shader */
		       "#version 330\n"
		       "out vec4 color;\n"
		       "uniform vec3 ID;\n"
		       "void main() {\n"
		       "    color = vec4(ID, 1.0);\n"
		       "}"
		       );

      mPickShader.bind();
      mPickShader.shareAttrib(mDispShader, "indices");
      mPickShader.shareAttrib(mDispShader, "vertexPosition_modelspace");
      mPickShader.setUniform("ID", mId);
      //mPickShader.release();
    };


  void updateModel(){
    //first rotation 
    mMatrix = mRot;
    Vec3 cen;
    //then translate into world coordinates
    
    //std::cout << mRot << std::endl;
    Vec4 nCen = mCen + mOffset;
    
    nCen[0] /= mScale[0];
    nCen[1] /= mScale[1];
    nCen[2] /= mScale[2];
    nCen = mRot.transpose()*nCen;
    mMatrix = nanogui::scale(mMatrix, xyz(mScale));
    mMatrix = nanogui::translate(mMatrix, xyz(nCen));
    
    //std::cout << mMatrix << std::endl;
  }
  
  void bind(){
      mDispShader.bind();
  }
  
  void highlightColor(){
    Vec3 hColor = mColor + Vec3(0.25,0.25, 0.25);
    mDispShader.bind();
    mDispShader.setUniform("uColor", hColor);
  }
  
  void selectColor(){
    Vec3 hColor(1.0,0,0);
    mDispShader.setUniform("uColor", hColor);
    mDispShader.bind();
  }

  void resetColor(){
    mDispShader.bind();
    mDispShader.setUniform("uColor", mColor);
  }

  void setColor(Vec3 s){
    mColor = s;
    mDispShader.bind();
    mDispShader.setUniform("uColor", mColor);
  }

  void setScale(Vec4 t){
    mScale = t;
    //mMatrix = nanogui::scale(mMatrix, xyz(t));
    updateModel();
  }

  void translate(Vec4 t){
    t(0) /= mScale[0];
    t(1) /= mScale[1];
    t(2) /= mScale[2];
    
    //t = mRot.transpose()*t;
    
    mCen(0) += t(0);
    mCen(1) += t(1);
    mCen(2) += t(2);
    //mMatrix = nanogui::translate(mMatrix, xyz(t));
    updateModel();
  }

  void setOffset(Vec4 nOff){
    mOffset = nOff;
    updateModel();
  }
  
  void setCenter(Vec4 nCen){
    mCen = nCen;
    updateModel();
  }
  
  void applyRotation(Mat4 r){
    mRot = r*mRot;
    updateModel();
  }

  void setRotation(Mat4 r){
    mRot = r;
    updateModel();
  }


  void draw(Mat4 & mProject, Mat4 & mModelView){
    this->displayShader().bind();
    Mat4 matrix = this->matrix();
    Mat4 mvp = mProject*mModelView*matrix;
    this->displayShader().setUniform("MVP", mvp);
    this->displayShader().setUniform("V", mModelView);
    this->displayShader().setUniform("M", this->matrix());
    this->displayShader().setUniform("LightPosition_worldspace", Vec3(3,3.,5.));
    
    /* Draw 2 triangles starting at index 0 */
    this->displayShader().drawIndexed(GL_TRIANGLES, 0, mIndices.cols());
  }
  //getters and setters.  Some variables can be gotten here, but only if they don'te
  //require special handling
  int   ID() const {return mId;}
  int & ID()       {return mId;}
  Mat4   matrix() const {return mMatrix;}
  Mat4 & matrix()       {return mMatrix;}
  
  Mat4   rmatrix() const {return mRot;}
  Mat4 & rmatrix()       {return mRot;}

  Vec4   center() const {return mCen;} 
  Vec4   scale()  const {return mScale;} 
  Vec4   offset()  const {return mOffset;} 
    
  nanogui::GLShader& displayShader(){ return mDispShader; }
  nanogui::GLShader& pickingShader(){ return mPickShader; }

  //make these public just to make life easy, not sure I'd do this in full production
  //though I like public variables.  Right now, the context
  std::function<void(Vec4 hitPoint, BufferObject & buffer, GnomonApplication & context, void * data)>  onSelect;
  std::function<void(Vec4 hitPoint, BufferObject & buffer, GnomonApplication & context, void * data)>  onDeselect;
  std::function<void(Vec4 dragStart, Vec4 dragDest, Vec4 cenStart, 
			BufferObject & buffer, GnomonApplication & context, 
			void * data)>  onDrag;
  Vec4 bbmin;
  Vec4 bbmax;
  int selectionGroup;
  bool isVisible;
private:
  int mId; //this is the color we use for picking
  Vec3 mColor; //this is the color we use for picking
  nanogui::MatrixXu mIndices;
  nanogui::MatrixXf mPositions;
  nanogui::MatrixXf mVertNormals;
  nanogui::GLShader mDispShader;
  nanogui::GLShader mPickShader;
  
  bool isHovered;
  bool isSelected;
  Vec4 mScale;
  Vec4 mCen;
  Mat4 mRot;
  Vec4 mOffset;
  Mat4 mMatrix;
};

BufferObject * makeCube(){  
  using namespace nanogui;
  nanogui::GLShader shader;
  BufferObject * cube = new BufferObject();
  cube->buildBuffer(12, 8,
		   [](nanogui::MatrixXu & indices, 
		      nanogui::MatrixXf & positions)->void {
		     indices.col(0)  << 0, 1, 2;
		     indices.col(1)  << 2, 3, 0;
		     indices.col(2)  << 0, 3, 4;
		     indices.col(3)  << 4, 5, 0;
		     indices.col(4)  << 0, 5, 6;
		     indices.col(5)  << 6, 1, 0;

		     indices.col(6)  << 1, 6, 7;
		     indices.col(7)  << 7, 2, 1;

		     indices.col(8)  << 7, 4, 3;
		     indices.col(9)  << 3, 2, 7;

		     indices.col(10) << 4, 7, 6;
		     indices.col(11) << 6, 5, 4;

		     positions.col(0) <<  1, 1, 1;
		     positions.col(1) << -1, 1, 1;
		     positions.col(2) << -1,-1, 1;
		     positions.col(3) <<  1,-1, 1;

		     positions.col(4) <<  1,-1,-1;
		     positions.col(5) <<  1, 1,-1;
		     positions.col(6) << -1, 1,-1;
		     positions.col(7) << -1,-1,-1;
		   });
  cube->initDisplayShader();
  cube->initPickingShader();
  
  return cube;
};


BufferObject * makeCone(){  
  using namespace nanogui;
  nanogui::GLShader shader;
  BufferObject * cone = new BufferObject();
  int nVerts = 16 + 1;
  int nIndices = nVerts - 1;
  cone->buildBuffer(nIndices, nVerts,
		   [&](nanogui::MatrixXu & indices, 
		       nanogui::MatrixXf & positions)->void {
		     
		     float radius = 0.1;
		     float height = 0.5;
		     for(int i = 0; i < nVerts - 1; i++){
		       float phi = 2.0*float(i)/float(nVerts-1)*M_PI;
		       float x = radius*cos(phi);
		       float y = radius*sin(phi);
		       positions.col(i) << x, y, 0;
		     }
		     //add the tip
		     positions.col(nVerts - 1) << 0,0,height;
		     
		     for(int i = 0; i < nIndices; i++){
		       //wrap to the beginning
		       indices.col(i) << i, (i+1)%nIndices, nVerts-1;
		     }
		   });
  cone->initDisplayShader();
  cone->initPickingShader();
  
  return cone;
};




  class GnomonApplication : public nanogui::Screen {
  public:
    typedef double Real;
  
    typedef Eigen::Matrix3f Mat3;
    typedef Eigen::Matrix4f Mat4;
    typedef Eigen::Vector2i Vec2i;
    typedef Eigen::Vector3f Vec3;
    typedef Eigen::Vector4f Vec4;
    typedef Eigen::Quaternion<float> Quat;

    GnomonApplication() : nanogui::Screen(Eigen::Vector2i(1600, 1600), "NanoGUI Test") {
      using namespace nanogui;
      performLayout(mNVGContext);
      mProject = this->makeProjectionMatrix(60.0*M_PI/180.0, 1.0, 0.1, 20);
      mModelView.setIdentity();
      mModelViewOld.setIdentity();
    
      mModelRotOld.setIdentity();
      mModelRotNew.setIdentity();

      //initialize selection indices
      mCurrentHover = -1;
      mCurrentGroup = 0;
      mCurrentSelect.push_back(-1);
      mCurrentSelect.push_back(-1);

      mPosition = Vec4(0,0,4,0);
      initObjects();
      ball.setSize(mSize);

      // Enable depth test
      glEnable(GL_DEPTH_TEST);
      // Accept fragment if it closer to the camera than the former one
      glDepthFunc(GL_LESS); 
      // Cull triangles which normal is not towards the camera
      glEnable(GL_CULL_FACE);

      //now for GUI
      Window *window = new Window(this, "coordinates");
      window->setPosition(Vector2i(15, 15));
      window->setLayout(new GroupLayout());
      
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


    }

    

    ~GnomonApplication() {
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
      
      //////////////////////////////////////////
      //ON DRAG: shape moves itself and helpers
      //////////////////////////////////////////
      cube.onDrag 
	= [&](Vec4 dragStart, Vec4 dragDest, Vec4 cenStart, 
	      BufferObject & buffer, GnomonApplication & context, 
	      void * data){
	Vec4 cen = buffer.center();
	Vec4 scale = buffer.scale();
	Vec4 dp = dragDest - dragStart;
	Vec4 nCen = cenStart + dp;
	buffer.setCenter(nCen);

	std::vector<BufferObject*> & helpers = context.getHelpers();
	for(int i = 0; i < helpers.size(); ++i){
	  helpers[i]->setCenter(nCen);
	  helpers[i]->isVisible = true;
	}
	context.updateTextBox(nCen);
      };

      
      //////////////////////////////////////////
      //ON SELECT: initializes helpers
      //////////////////////////////////////////

      cube.onSelect = [&](Vec4 hitPoint, BufferObject & buffer, 
			  GnomonApplication & context, void * data){
	Vec4 cen = buffer.center();
	Vec4 scale = buffer.scale();
	Real scaleAvg = 1.0/2.0*(scale[0] + scale[1] + scale[2]);
	  
	buffer.displayShader().bind();
	buffer.selectColor();

	std::vector<BufferObject*> & helpers = context.getHelpers();
	for(int i = 0; i < helpers.size(); ++i){
	  
	  helpers[i]->setScale(Vec4(scaleAvg,scaleAvg,scaleAvg,1));
	  
	  //now we'll randomize the direction vectors
	  //because hey, they translate along their orientation
	  //they can be in any direction they want
	  Mat4 R = RandomRotation();  
	  Vec4 off = R*Vec4(0,0,1,1);
	
	  helpers[i]->setRotation(R);
	  helpers[i]->setCenter(cen);
	  helpers[i]->setOffset(scaleAvg*off);
	  helpers[i]->setColor(xyz(off));
	  helpers[i]->isVisible = true;
	}
      };

      //////////////////////////////////////////
      //ON DESELECT: hides helpers
      //////////////////////////////////////////

      cube.onDeselect = [&](Vec4 hitPoint, BufferObject & buffer, 
			    GnomonApplication & context, void * data){
	
	buffer.displayShader().bind();
	buffer.resetColor();
      
	std::vector<BufferObject*> & helpers = context.getHelpers();
	for(int i = 0; i < helpers.size(); i++)
	  helpers[i]->isVisible = false;
	context.updateTextBox(Vec4(0,0,0,0));
	
      };
      
    }

    BufferObject * translationHelper(Mat4 R){
      
      
      auto onDrag = [](Vec4 dragStart, Vec4 dragDest, Vec4 cenStart, 
			BufferObject & hBuffer, GnomonApplication & context, 
			void * data){
	
	BufferObject & oBuffer = *static_cast<BufferObject*> (data);
	Vec4 hcen = hBuffer.center();
	Mat4 R = hBuffer.rmatrix();
	//the drag direction is defined by the widgets rotation... these
	//could be arbitrary... hmmmm
	Vec4 dragDir = R*Vec4(0,0,-1,1);
	
	Vec4 pDest = projectToLine(dragDest, hcen, hcen + dragDir);
	Vec4 pStart = projectToLine(dragStart, hcen, hcen + dragDir);

	//this is a little hacky the widget just calls the 'parent' buffers onDrag
	//the parent then moves all the widgets
	oBuffer.onDrag(pDest, pStart, cenStart, oBuffer, context, NULL);
      };

      BufferObject * helper = makeCone();
      R.setIdentity();
      helper->setRotation(R);
      //Mat4 R = hBuffer.rmatrix();
      Vec4 axis = R*Vec4(0,0,1,1);
      helper->setColor(xyz(axis));
	
      helper->setOffset(axis);
      helper->isVisible = false;
      //widgets will be handled by seperating out the selection groups
      helper->selectionGroup = 1;
      
      helper->onDrag = onDrag;
      return helper;
    }

    void initHelperWidgets(){
      
      auto fRand = [](double fMin, double fMax)->double{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
      };
      
      //Just for fun, we'll initialize 24 random widgets helper widgets
      //who says orthoganal directions are best.
      for(int i = 0; i < 24; i++){
	
	Mat4 R4 = RandomRotation();
	BufferObject * helper = translationHelper(R4);
	
	mHelpers.push_back(helper);
	mSceneObjects.push_back(helper);
      }
    }

    void initObjects(){      
      
      for(int i=0; i < 64; i++){
	BufferObject * cube = makeCube();	
	initCube(*cube);
	mSceneObjects.push_back(cube);
      }
      initHelperWidgets();
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
      ball.button(p,down);
      
      Vec4 hit;
      int oldSelection = mCurrentSelect[mCurrentGroup];
      int newSelection = pickObject(p, hit);
      int curSelection;
      
      if(down){
	if(newSelection != -1){
	  int oldGroup = mCurrentGroup;
	  mCurrentGroup = mSceneObjects[newSelection]->selectionGroup;
	  
	  curSelection = newSelection;
	  
	  mCurrentSelect[mCurrentGroup] = curSelection;

	  if(oldSelection != curSelection){
	    BufferObject & buffer = *mSceneObjects[curSelection];
	    
	    //If we've chosen a new object of the same group, then deselect the old
	    //object... Make sure that the old selection is valid!!!
	    if(oldSelection != mCurrentSelect[mCurrentGroup] &&
	       oldGroup == mCurrentGroup &&
	       oldSelection != -1){
	      BufferObject & oldBuffer = *mSceneObjects[oldSelection];
	      if(oldBuffer.onDeselect)
		oldBuffer.onDeselect(hit, oldBuffer, *this, NULL);
	    }

	    //If we've chose a new object, select it and then choose the 
	    //proper event callback invocation for the group.
	    //in this case, group 1 is widgets and group 0 is objects.
	    switch(mCurrentGroup){
	    case 0:
	      if(buffer.onSelect)
		buffer.onSelect(hit, buffer, *this, NULL);
	      break;
	    case 1: //interaction widgets need to send over a pointer to the currently
	      //selected buffer
	      if(buffer.onSelect)
		buffer.onSelect(hit, buffer, *this, mSceneObjects[mCurrentSelect[0]]);
	      break;
	    }
	  }
	 	  
	
	  //this is for dragging the object later
	  //Need to cache the object's center so we can apply our dragged offsets
	  //to it during dragging
	  mObjCenCache = mSceneObjects[mCurrentSelect[mCurrentGroup]]->center();
	  
	  mDragStart = hit;
	  
	  mDragging = true;
	 
	}
	//we also need to deselect if we don't touch anything
	//std::cout << newSelection << " " << oldSelection << " " << mCurrentSelect[0] << std::endl;
	if(newSelection == -1 &&
	   oldSelection != -1){
	  BufferObject & oldBuffer = *mSceneObjects[mCurrentSelect[0]];
	  if(oldBuffer.onDeselect)
	    oldBuffer.onDeselect(hit, oldBuffer, *this, NULL);
	  
	  mCurrentSelect[0] = -1;
	  mCurrentSelect[1] = -1;
	  mDragging = false;
	}

      } else{
	
	//BufferObject & buffer = *mSceneObjects[mCurrentSelect];
	
	updateFrame(); 
	mDragging = false;
      }
      
      Screen::mouseButtonEvent(p, button, down,  modifiers);
      return true;
    }
    
    
    
    virtual bool mouseMotionEvent(const Eigen::Vector2i &p, 
				  const Eigen::Vector2i &rel, 
				  int button, int modifiers) {
      ball.motion(Vec2i(p[0], p[1]));
      
      //perform hovering
      int oldHover = mCurrentHover;
      Vec4 hit;
      //handling hover should probably be accomplished with events, but I've explicitely
      //handled it here to reduce programming
      
      //Prevent the camera from moving during a drag
      if(!mDragging)
	mCurrentHover = pickObject(p, hit);
      //Hovering logic:
      if(mCurrentHover != -1 && mCurrentSelect[0]== -1){
	mSceneObjects[mCurrentHover]->displayShader().bind();
	mSceneObjects[mCurrentHover]->highlightColor();
      }
      if(mCurrentHover != oldHover && oldHover != -1
	 && mCurrentSelect[mCurrentGroup] == -1){
	mSceneObjects[oldHover]->displayShader().bind();
	mSceneObjects[oldHover]->resetColor();
      }

      if(mDragging && mCurrentSelect[mCurrentGroup] > -1 && 
	 oldHover == mCurrentHover){
	
	BufferObject & buffer = *mSceneObjects[mCurrentSelect[mCurrentGroup]];
	Vec4 wCamNormal = castRay(Vec2i(mSize[0]/2.0, 
						mSize[1]/2.0));
	//we'll project the ray onto the camera plane, 
	//in order to move the object on the plane
	Vec4 dragDest= rayPlaneIntersect(mPosition, 
					 castRay(p), 
					 wCamNormal, 
					 mDragStart);
	
	switch(mCurrentGroup){
	    case 0:
	      if(buffer.onDrag)
		buffer.onDrag(mDragStart,dragDest,mObjCenCache, 
			      buffer, *this, NULL);
	      break;
	    case 1: //interaction widgets need to send over a pointer to the currently
	      //selected buffer
	      if(buffer.onDrag)
	      	buffer.onDrag(mDragStart,dragDest,mObjCenCache, 
			      buffer, *this, mSceneObjects[mCurrentSelect[0]]);
	      
	      break;
	    }
	
      } else {
	updatePosition();
      }
      
      Screen::mouseMotionEvent(p, rel, button, modifiers);
      return true;
    }

    ///////////////////////////
    //  Draw
    ///////////////////////////


    void updateFrame(){
      //cache the old rotation
      mModelRotOld = mModelRotNew;
      ball.setState(Quat::Identity()); //reset ball to identity
    }

    void updatePosition(){
      //This is my camera model, its not a seperat object, but doesn't necessarily
      //warrant one.  All the fancy stuff was implemented in the arcball class
      Mat4 r = ball.matrix(Mat4());
      //accumulate rotations... this likely could be done with the quaternions in the 
      //arcball, but this works after some fiddling
      mModelRotNew = r*mModelRotOld;
      //set the modelView equal to the new rotation
      mModelView = mModelRotNew;
      mModelView.col(3) << Vec4(0,0,-4,1); 
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
#if GUI
      Screen::draw(ctx);
#endif
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
      Vec4 wPosition;
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
    std::vector<BufferObject*> mSceneObjects;
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
      

    nanogui::Arcball ball;

  };

  int main(int argc, char *argv[]) {
    try {
      nanogui::init();

      GnomonApplication *app = new GnomonApplication();
      app->drawAll();
      app->setVisible(true);

      nanogui::mainloop();

      delete app;

      nanogui::shutdown();
    } catch (const std::runtime_error &e) {
      std::string error_msg = std::string("Caught a fatal error: ") + std::string(e.what());
#if defined(WIN32)
      MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#else
      std::cerr << error_msg << endl;
#endif
      return -1;
    }

    return 0;
  }
