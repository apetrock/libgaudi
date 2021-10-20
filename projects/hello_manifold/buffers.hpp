#ifndef __BUFFEROBJECT__
#define __BUFFEROBJECT__

//#include <nanogui/glutil.h>
#include <iostream>
#include <string>


//nano gui deals mostly with vec3s, which is annoying
//since I'm adhering to using homogeneous coordinates
//so far as I can tell, Eigen can grab blocks, but
//that syntax is laborious, so here as a few helper functions
Eigen::Vector3f xyz(Eigen::Vector4d in){
  return Eigen::Vector3f(in[0],in[1],in[2]);
};

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
  return convert.str(); // set 'Result' to the contents of the 
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
  typedef Eigen::Vector3i Vec3i;
  
  typedef Eigen::Vector3f Vec3;
  typedef Eigen::Vector4f Vec4;
  typedef Eigen::Quaternion<float> Quat;
  typedef nanogui::MatrixXu MatrixXu;
  typedef nanogui::MatrixXf MatrixXf;
  
  Drawable(){}
  ~Drawable(){}
  virtual void draw(){}
  virtual bool intersectBbox(Vec4 r0, Vec4 r1, Real & tnear, Real & tfar ){return true;};
  
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
  BufferObject(){
  };
  ~BufferObject(){
    //free everything here
    
    mDispShader->free();
    mPickShader->free();
  };
  
    void buildBuffer(int nInd, int nVerts, 
		     std::function<void(nanogui::MatrixXu&, 
					nanogui::MatrixXf&)> buildFunc){
      
      mDispShader = new nanogui::GLShader;
      mPickShader = new nanogui::GLShader;
  
      mIndices   = nanogui::MatrixXu(3, nInd);
      mPositions = nanogui::MatrixXf(3, nVerts);
      
      buildFunc(mIndices, mPositions);
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
      //std::cout << mIndices.transpose() << std::endl;
      for(int i = 0; i < mIndices.cols(); i++){
        //std::cout << mIndices.size() << " " << mIndices(0,i) << std::endl;
        //std::cout << mPositions.col(mIndices(0,i)).transpose() << std::endl;
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
  

    virtual void initDisplayShader(){
      //std::cout << " system: ";
      //system("less src/shaders/standard.vs");
      mId = rand();
      int Number = mId;       // number to be converted to a string
      std::string Result;          // string which will contain the result
      std::ostringstream convert;   // stream used for the conversion
      convert << "a_standard_shader: " << Number;      // insert the textual representation of 'Number' in the characters in the stream
      Result = convert.str();
      
      mDispShader->initFromFiles(
				/* An identifying name */
				convert.str(),
				"src/shaders/standard.vs",
				"src/shaders/standard.fs");
      mDispShader->bind();
      computeNormals();
      
      mDispShader->uploadAttrib("vertexNormal_modelspace", mVertNormals);
      mDispShader->uploadIndices(mIndices);
      mDispShader->uploadAttrib("vertexPosition_modelspace", mPositions);
      mDispShader->setUniform("uColor", mColor);
    };
  
  virtual void initPickingShader(){
      //I couldn't get backbuffer style picking to work in short time, so I abandoned it,
      //issues with opengl set or nanoguis wrapper wasn't working, though I tried
      //using raw openGL as well. I imagine that its my fault somewhere, but I didn't
      //want to spend too much longer debugging it
      //here is the legacy code:
      mPickShader->init(
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

      mPickShader->bind();
      mPickShader->shareAttrib(*mDispShader, "indices");
      mPickShader->shareAttrib(*mDispShader, "vertexPosition_modelspace");
      mPickShader->setUniform("ID", mId);
      //mPickShader->release();
    };


  virtual void updateModel(){
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
    mMatrix = mMatrix * nanogui::scale(xyz(mScale));
    mMatrix = mMatrix * nanogui::translate(xyz(nCen));
    
    //std::cout << mMatrix << std::endl;
  }
  
  void bind(){
      mDispShader->bind();
  }
  
  virtual void highlightColor(){
    Vec3 hColor = mColor + Vec3(0.25,0.25, 0.25);
    mDispShader->bind();
    mDispShader->setUniform("uColor", hColor);
  }
  
  virtual void selectColor(){
    Vec3 hColor(1.0,0,0);
    mDispShader->setUniform("uColor", hColor);
    mDispShader->bind();
  }

  virtual void resetColor(){
    mDispShader->bind();
    mDispShader->setUniform("uColor", mColor);
  }

  virtual void setColor(Vec3 s){
    mColor = s;
    mDispShader->bind();
    mDispShader->setUniform("uColor", mColor);
  }

  virtual void setScale(Vec4 t){
    mScale = t;
    //mMatrix = nanogui::scale(mMatrix, xyz(t));
    updateModel();
  }

  virtual void translate(Vec4 t){
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

  virtual void setOffset(Vec4 nOff){
    mOffset = nOff;
    updateModel();
  }
  
  virtual void setCenter(Vec4 nCen){
    mCen = nCen;
    updateModel();
  }
  
  virtual void applyRotation(Mat4 r){
    mRot = r*mRot;
    updateModel();
  }

  virtual void setRotation(Mat4 r){
    mRot = r;
    updateModel();
  }


  virtual void draw(Mat4 & mProject, Mat4 & mModelView){
    this->displayShader().bind();
    Mat4 matrix = this->matrix();
    Mat4 mvp = mProject*mModelView*matrix;
    std::cout << mvp << std::endl;
    this->displayShader().setUniform("MVP", mvp);
    this->displayShader().setUniform("V", mModelView);
    this->displayShader().setUniform("M", this->matrix());
    this->displayShader().setUniform("LightPosition_worldspace", Vec3(3,3.,5.));
    
    /* Draw 2 triangles starting at index 0 */
    this->displayShader().drawIndexed(GL_TRIANGLES, 0, mIndices.cols());
  }

  //getters and setters.  Some variables can be gotten here, but only if they don'te
  //require special handling
  virtual int   ID() const {return mId;}
  virtual int & ID()       {return mId;}
  virtual Mat4   matrix() const {return mMatrix;}
  virtual Mat4 & matrix()       {return mMatrix;}
  
  virtual Mat4   rmatrix() const {return mRot;}
  virtual Mat4 & rmatrix()       {return mRot;}

  virtual Vec4   center() const {return mCen;} 
  virtual Vec4   scale()  const {return mScale;} 
  virtual Vec4   offset()  const {return mOffset;} 
    
  nanogui::GLShader& displayShader(){ return *mDispShader; }
  nanogui::GLShader& pickingShader(){ return *mPickShader; }

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
protected:
  int mId; //this is the color we use for picking
  Vec3 mColor; //this is the color we use for picking
  nanogui::MatrixXu mIndices;
  nanogui::MatrixXf mPositions;
  nanogui::MatrixXf mVertNormals;

  nanogui::GLShader * mDispShader;
  nanogui::GLShader * mPickShader;
  
  bool isHovered;
  bool isSelected;
  Vec4 mScale;
  Vec4 mCen;
  Mat4 mRot;
  Vec4 mOffset;
  Mat4 mMatrix;
};


class MeshObject : public BufferObject {
public:

};

class ImmediateLines : public BufferObject{
public:
  ImmediateLines(){
  };
  ~ImmediateLines(){
    //free everything here
  };
  

  void buildBuffer(){
    //don't need to override default buildBuffer, since we do all our counting internally
    mDispShader = new nanogui::GLShader;
    
    int ni = 0, nv = 0;
    for(int i = 0; i < mImIndices.size(); i++){
      ni += mImIndices[i].cols();
    }

    for(int i = 0; i < mImPositions.size(); i++){
      nv += mImPositions[i].cols();
    }

    mIndices   = nanogui::MatrixXu(2, ni);
    mPositions = nanogui::MatrixXf(3, nv);
    
    int iI = 0;
    int iP = 0;
    for(int i = 0; i < mImIndices.size(); i++){
      for(int j = 0; j < mImIndices[i].cols(); j++){
	
	mIndices.col(iI + j) << (Vec2i(iP,iP).cast<unsigned int>() + mImIndices[i].col(j));
      }
      iP += mImPositions[i].cols();
      iI += mImIndices[i].cols();
    }

    
    iP = 0;
    for(int i = 0; i < mImPositions.size(); i++){
      for(int j = 0; j < mImPositions[i].cols(); j++){
	mPositions.col(iP + j) << mImPositions[i].col(j);
      }
      iP += mImPositions[i].cols();
    }
    
    //std::cout << mIndices << std::endl;
    //std::cout << mPositions << std::endl;
    calcBbox();
      
    mMatrix.setIdentity();
    mRot.setIdentity();
    mColor = Vec3(1.00, 0.0, 0.0);
    mCen = Vec4(0.0,0.0,0.0,1.0);//add the homogeneous 1
    mScale = Vec4(1.0,1.0,1.0,1.0);
    mOffset = Vec4(0.0,0.0,0.0,1.0);
    selectionGroup = 0;
    isVisible = true;
  };
  
    void initDisplayShader(){
      //std::cout << " system: ";
      //system("less src/shaders/standard.vs");
      mDispShader->init(
			/* An identifying name */
			"a_line_shader",

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
			"uniform vec3 uColor;\n"
			"void main() {\n"
			"    color = vec4(uColor, 1.0);\n"
			"}"
			);

      mDispShader->bind();
      mDispShader->uploadIndices(mIndices);
      mDispShader->uploadAttrib("vertexPosition_modelspace", mPositions);
      mDispShader->setUniform("uColor", mColor);
    };


  void pushBox(Vec3 cen, Vec3 h){
    nanogui::MatrixXu indices   = nanogui::MatrixXu(2, 12);
    nanogui::MatrixXf positions = nanogui::MatrixXf(3, 8);

    indices.col(0)  << 0, 1;
    indices.col(1)  << 1, 2;
    indices.col(2)  << 2, 3;
    indices.col(3)  << 3, 0;
		     
    indices.col(4)  << 0, 5;
    indices.col(5)  << 1, 6;
    indices.col(6)  << 2, 7;
    indices.col(7)  << 3, 4;

    indices.col(8)   << 4, 5;
    indices.col(9)   << 5, 6;
    indices.col(10)  << 6, 7;
    indices.col(11)  << 7, 4;

    positions.col(0) << cen + Vec3( h[0], h[1], h[2]);
    positions.col(1) << cen + Vec3(-h[0], h[1], h[2]);
    positions.col(2) << cen + Vec3(-h[0],-h[1], h[2]);
    positions.col(3) << cen + Vec3( h[0],-h[1], h[2]);

    positions.col(4) << cen + Vec3( h[0],-h[1],-h[2]);
    positions.col(5) << cen + Vec3( h[0], h[1],-h[2]);
    positions.col(6) << cen + Vec3(-h[0], h[1],-h[2]);
    positions.col(7) << cen + Vec3(-h[0],-h[1],-h[2]);
    
    mImPositions.push_back(positions);
    mImIndices.push_back(indices);
  };

  void pushLine(Vec3 c0, Vec3 c1){
    nanogui::MatrixXu indices   = nanogui::MatrixXu(2, 1);
    nanogui::MatrixXf positions = nanogui::MatrixXf(3, 2);

    indices.col(0)  << 0, 1;
    
    positions.col(0) << c0;
    positions.col(1) << c1;
    mImPositions.push_back(positions);
    mImIndices.push_back(indices);

  };

  virtual void draw(Mat4 & mProject, Mat4 & mModelView){
    this->displayShader().bind();
    Mat4 matrix = this->matrix();
    Mat4 mvp = mProject*mModelView*matrix;
    this->displayShader().setUniform("MVP", mvp);
    
    /* Draw 2 triangles starting at index 0 */
    glLineWidth(5.0);
    this->displayShader().drawIndexed(GL_LINES, 0, mIndices.cols());
  }

private:
  std::vector<nanogui::MatrixXf> mImPositions;
  std::vector<nanogui::MatrixXu> mImIndices;
};



BufferObject * makeCube(){  
  using namespace nanogui;

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
  //cube->initPickingShader();
  
  return cube;
};
/*
BufferObject * makeLineCube(){  
  using namespace nanogui;

  BufferObject * cube = new BufferObject();
  cube->buildLineBuffer(12, 8,
		   [](nanogui::MatrixXu & indices, 
		      nanogui::MatrixXf & positions)->void {
		     indices.col(0)  << 0, 1;
		     indices.col(1)  << 1, 2;
		     indices.col(2)  << 2, 3;
		     indices.col(3)  << 3, 0;
		     
		     indices.col(4)  << 0, 5;
		     indices.col(5)  << 1, 6;
		     indices.col(6)  << 2, 7;
		     indices.col(7)  << 3, 4;

		     indices.col(8)   << 4, 5;
		     indices.col(9)   << 5, 6;
		     indices.col(10)  << 6, 7;
		     indices.col(11)  << 7, 4;

		     positions.col(0) <<  1, 1, 1;
		     positions.col(1) << -1, 1, 1;
		     positions.col(2) << -1,-1, 1;
		     positions.col(3) <<  1,-1, 1;

		     positions.col(4) <<  1,-1,-1;
		     positions.col(5) <<  1, 1,-1;
		     positions.col(6) << -1, 1,-1;
		     positions.col(7) << -1,-1,-1;
		   });
  //cube->initDisplayShader();
  cube->initLineShader();
  
  return cube;
};
*/

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
  //cone->initPickingShader();
  
  return cone;
};
#endif
