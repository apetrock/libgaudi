#ifndef __BUFFEROBJECT__
#define __BUFFEROBJECT__

#include <GaudiGraphics/viewer.hpp>
#include <GaudiMath/typedefs.h>
#include <nanogui/glutil.h>

#include <iostream>
#include <string>

namespace gg {
using namespace GaudiMath;
std::string realToString(float f) {
  std::ostringstream convert; // stream used for the conversion
  convert << f;         // insert the textual representation of 'Number' in the
  return convert.str(); // set 'Result' to the contents of the
}

class Drawable;
using DrawablePtr = std::shared_ptr<Drawable>;

class Drawable {
  // Really simple hierarchy, all classes derive from drawable, and thats really
  // the only method they have, so that we can collect all drawables
public:
  Drawable() {}
  ~Drawable() {}

  static DrawablePtr create() { return std::make_shared<Drawable>(); }

  virtual void draw(Mat4 &mProject, Mat4 &mModelView) {}
  virtual void init() {}

  virtual bool intersectBbox(Vec4 r0, Vec4 r1, Real &tnear, Real &tfar) {
    return true;
  };

  virtual bool rayBoxIntersect(Vec4 r_o, Vec4 r_d, Vec4 boxmin, Vec4 boxmax,
                               Real &tnear, Real &tfar) {

    Vec4 invR = Vec4(1.0 / r_d[0], 1.0 / r_d[1], 1.0 / r_d[2], 0.0);

    Vec4 tbot, ttop;
    for (int i = 0; i < 4; i++) {
      tbot[i] = invR[i] * (boxmin[i] - r_o[i]);
      ttop[i] = invR[i] * (boxmax[i] - r_o[i]);
    }
    // re-order intersections to find smallest and largest on each axis
    Vec4 tmin, tmax;
    for (int i = 0; i < 4; i++) {
      tmin[i] = std::min(ttop[i], tbot[i]);
      tmax[i] = std::max(ttop[i], tbot[i]);
    }

    // find the largest tmin and the smallest tmax
    Real largest_tmin =
        std::max(std::max(tmin[0], tmin[1]), std::max(tmin[0], tmin[2]));
    Real smallest_tmax =
        std::min(std::min(tmax[0], tmax[1]), std::min(tmax[0], tmax[2]));

    tnear = largest_tmin;
    tfar = smallest_tmax;

    return (smallest_tmax > largest_tmin);
  }

  bool isVisible = true;
};

class BufferObject;
using BufferObjectPtr = std::shared_ptr<BufferObject>;

class BufferObject : public Drawable {
public:
  static BufferObjectPtr create() { return std::make_shared<BufferObject>(); }

  BufferObject(){};
  ~BufferObject() {
    // free everything here

    mDispShader->free();
  };

  virtual void init() { this->initBuffer(); }

  virtual void initBuffer() {

    mDispShader = new nanogui::GLShader;

    mIndices = nanogui::MatrixXu(3, 0);
    mPositions = nanogui::MatrixXf(3, 0);

    mMatrix.setIdentity();
    mRot.setIdentity();
    mColor = Vec3(0.75, 0.75, 1.0);
    mCen = Vec4(0.0, 0.0, 0.0, 1.0); // add the homogeneous 1
    mScale = Vec4(1.0, 1.0, 1.0, 1.0);
    mOffset = Vec4(0.0, 0.0, 0.0, 1.0);
    selectionGroup = 0;
    isVisible = true;
    this->initDisplayShader();
  };

  virtual void allocateVerts(int nInd, int nVerts) {
    mIndices.resize(3, nInd);
    mPositions.resize(3, nVerts);
    mColors.resize(3, nVerts);
  };

  void fillBuffer(std::function<void(BufferObject &)> buildFunc) {
    buildFunc(*this);
    calcBbox();
    computeNormals();
  };

  virtual bool intersectBbox(Vec4 r0, Vec4 r1, Real &tnear, Real &tfar,
                             Vec4 &hNear, Vec4 &hFar) {
    // the model may have a local transform, so rotate the rays
    // so that they are in model space

    Mat4 iMatrix = mMatrix.inverse();
    Vec4 r0p = iMatrix * r0; // p for prime
    Vec4 r1p = iMatrix * r1;

    Real tpnear, tpfar;
    bool hit = rayBoxIntersect(r0p, r1p, bbmin, bbmax, tpnear, tpfar);
    Vec4 hpNear = (r0p + tpnear * r1p);
    Vec4 hpFar = (r0p + tpfar * r1p);

    hNear = mMatrix * hpNear;
    hFar = mMatrix * hpFar;
    // we need to reclaim our ts in the new coordinates
    tnear = (hNear - r0).norm() / r1.norm(); // these should be the same
    tfar = (hFar - r0).norm() / r1.norm();   // as the old ones
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

  void calcBbox() {
    bbmin = Vec4(99999, 99999, 99999, 1);
    bbmax = Vec4(-99999, -99999, -99999, 1);
    for (int i = 0; i < mPositions.cols(); i++) {
      Vec3 v = mPositions.col(i);
      bbmin(Eigen::seq(0, 2), 0) =
          bbmin(Eigen::seq(0, 2), 0).array().min(v.array());
      bbmax(Eigen::seq(0, 2), 0) =
          bbmax(Eigen::seq(0, 2), 0).array().min(v.array());
    }
  }

  void computeNormals() {
    // this isn't really that necessary, now, normals are computed
    // per fragment on the shader
    mVertNormals = nanogui::MatrixXf(mPositions.rows(), mPositions.cols());
    mVertNormals.setZero();

    for (int i = 0; i < mIndices.cols(); i++) {
      /*
      std::cout << i << " " << mIndices.cols() << " - " << mIndices(0, i) << " "
                << mIndices(1, i) << " " << mIndices(2, i) << " "
                << mPositions.rows() << " " << mPositions.cols() << std::endl;
      */
      Vec3 v0 = mPositions.col(mIndices(0, i));
      Vec3 v1 = mPositions.col(mIndices(1, i));
      Vec3 v2 = mPositions.col(mIndices(2, i));

      Vec3 N = (v1 - v0).cross(v2 - v0);
      N.normalize();
      mVertNormals.col(mIndices(0, i)) += N;
      mVertNormals.col(mIndices(1, i)) += N;
      mVertNormals.col(mIndices(2, i)) += N;
    }

    for (int i = 0; i < mVertNormals.cols(); i++) {
      mVertNormals.col(i).normalize();
    }
  }

  virtual void updateShaderAttributes() {
    //mDispShader->uploadAttrib("vertexNormal_modelspace", mVertNormals);
    mDispShader->uploadIndices(mIndices);
    mDispShader->uploadAttrib("vertexPosition_modelspace", mPositions);
    mDispShader->uploadAttrib("vertexColor_modelspace", mColors);
  }

  virtual void initDisplayShader() {
    // std::cout << " system: ";
    // system("less ../../src/shaders/standard.vs");
    mId = rand();
    int Number = mId;           // number to be converted to a string
    std::string Result;         // string which will contain the result
    std::ostringstream convert; // stream used for the conversion
    convert << "a_standard_shader: "
            << Number; // insert the textual representation of 'Number' in the
                       // characters in the stream
    Result = convert.str();

    mDispShader->initFromFiles(
        /* An identifying name */
        convert.str(), "src/shaders/standard.vs", "src/shaders/standard.fs");
  };

  virtual void updateModel() {
    // first rotation
    mMatrix = mRot;
    Vec3 cen;
    // then translate into world coordinates

    // std::cout << mRot << std::endl;
    Vec4 nCen = mCen + mOffset;

    nCen[0] /= mScale[0];
    nCen[1] /= mScale[1];
    nCen[2] /= mScale[2];
    nCen = mRot.transpose() * nCen;
    mMatrix = nanogui::scale(xyz(mScale));
    mMatrix = nanogui::translate(xyz(nCen));

    // std::cout << mMatrix << std::endl;
  }

  void bind() { mDispShader->bind(); }

  virtual void setScale(Vec4 t) {
    mScale = t;
    // mMatrix = nanogui::scale(mMatrix, xyz(t));
    updateModel();
  }

  virtual void translate(Vec4 t) {
    t(0) /= mScale[0];
    t(1) /= mScale[1];
    t(2) /= mScale[2];

    // t = mRot.transpose()*t;

    mCen(0) += t(0);
    mCen(1) += t(1);
    mCen(2) += t(2);
    // mMatrix = nanogui::translate(mMatrix, xyz(t));
    updateModel();
  }

  virtual void setOffset(Vec4 nOff) {
    mOffset = nOff;
    updateModel();
  }

  virtual void setCenter(Vec4 nCen) {
    mCen = nCen;
    updateModel();
  }

  virtual void applyRotation(Mat4 r) {
    mRot = r * mRot;
    updateModel();
  }

  virtual void setRotation(Mat4 r) {
    mRot = r;
    updateModel();
  }

  virtual void draw(Mat4 &mProject, Mat4 &mModelView) {
    bind();
    Mat4 matrix = this->matrix();
    Mat4 mvp = mProject * mModelView * matrix;
    this->displayShader().setUniform("MVP", mvp);
    this->displayShader().setUniform("V", mModelView);
    this->displayShader().setUniform("M", this->matrix());
    this->displayShader().setUniform("LightPosition_worldspace",
                                     Vec3(3, 3., 5.));
    updateShaderAttributes();
    /* Draw 2 triangles starting at index 0 */
    //
    // this->displayShader().drawIndexed(GL_LINES, 0, mIndices.cols());
    this->displayShader().drawIndexed(GL_TRIANGLES, 0, mIndices.cols());
  }

  // getters and setters.  Some variables can be gotten here, but only if they
  // don'te require special handling
  virtual int ID() const { return mId; }
  virtual int &ID() { return mId; }
  virtual Mat4 matrix() const { return mMatrix; }
  virtual Mat4 &matrix() { return mMatrix; }

  virtual Mat4 rmatrix() const { return mRot; }
  virtual Mat4 &rmatrix() { return mRot; }

  virtual Vec4 center() const { return mCen; }
  virtual Vec4 scale() const { return mScale; }
  virtual Vec4 offset() const { return mOffset; }

  nanogui::GLShader &displayShader() { return *mDispShader; }
  nanogui::MatrixXu &indices() { return mIndices; }
  nanogui::MatrixXf &positions() { return mPositions; }
  nanogui::MatrixXf &colors() { return mColors; }

  Vec4 bbmin;
  Vec4 bbmax;
  int selectionGroup;

protected:
  int mId;     // this is the color we use for picking
  Vec3 mColor; // this is the color we use for picking
  nanogui::MatrixXu mIndices;
  nanogui::MatrixXf mPositions;
  nanogui::MatrixXf mColors;
  nanogui::MatrixXf mVertNormals;

  nanogui::GLShader *mDispShader;

  bool isHovered;
  bool isSelected;
  Vec4 mScale;
  Vec4 mCen;
  Mat4 mRot;
  Vec4 mOffset;
  Mat4 mMatrix;
};

class PointBuffer;
using PointBufferPtr = std::shared_ptr<PointBuffer>;

class PointBuffer : public BufferObject {
public:
  static PointBufferPtr create() { return std::make_shared<PointBuffer>(); }

  PointBuffer(){};
  ~PointBuffer(){
      // free everything here
  };

  virtual void allocateVerts(int nInd, int nVerts) {
    mIndices.resize(1, nInd);
    mPositions.resize(3, nVerts);
    mColors.resize(3, nVerts);
  };

  void initDisplayShader() {
    // std::cout << " system: ";
    // system("less src/shaders/standard.vs");
    mDispShader->init(
        /* An identifying name */
        "a_point_shader",

        /* Vertex shader */

        "#version 330\n"
        "layout(location = 0) in vec3 vertexPosition_modelspace;\n"
        "layout(location = 1) in vec2 vertexUV;\n"
        "layout(location = 3) in vec3 vertexColor_modelspace;\n"
        "out vec2 UV;\n"
        "out vec3 Color_cameraspace;\n"
        "uniform mat4 MVP;\n"

        "void main() {\n"

        "    Color_cameraspace = vertexColor_modelspace;\n"
        "    UV = vertexUV;\n"
        "    gl_Position = MVP * vec4(vertexPosition_modelspace, 1.0);\n"
        "}",

        /* Fragment shader */
        "#version 330\n"
        "in vec2 UV;\n"
        "in vec3 Position_worldspace;\n"
        "in vec3 Color_cameraspace;\n"

        "out vec4 color;\n"
        "void main() {\n"
        "    color = vec4(Color_cameraspace, 1.0);\n"
        "}");

    mDispShader->bind();
    mDispShader->uploadIndices(mIndices);
    mDispShader->uploadAttrib("vertexPosition_modelspace", mPositions);
  };

  void fillBuffer(std::function<void(BufferObject &)> buildFunc) {
    buildFunc(*this);
    calcBbox();
    updateShaderAttributes();
  };

  virtual void updateShaderAttributes() {
    mDispShader->uploadIndices(mIndices);
    mDispShader->uploadAttrib("vertexPosition_modelspace", mPositions);
    mDispShader->uploadAttrib("vertexColor_modelspace", mColors);
  }

  virtual void draw(Mat4 &mProject, Mat4 &mModelView) {
    bind();
    Mat4 matrix = this->matrix();
    Mat4 mvp = mProject * mModelView * matrix;
    this->displayShader().setUniform("MVP", mvp);
    updateShaderAttributes();

    glPointSize(1.5);
    this->displayShader().drawIndexed(GL_POINTS, 0, mIndices.cols());
  }
}; // PointBuffer;

class LineBuffer;
using LineBufferPtr = std::shared_ptr<LineBuffer>;

class LineBuffer : public PointBuffer {
public:
  static LineBufferPtr create() { return std::make_shared<LineBuffer>(); }

  LineBuffer(){};
  ~LineBuffer(){
      // free everything here
  };

  virtual void allocateVerts(int nInd, int nVerts) {
    mIndices.resize(2, nInd);
    mPositions.resize(3, nVerts);
    mColors.resize(3, nVerts);
  };

  void fillBuffer(std::function<void(BufferObject &)> buildFunc) {
    buildFunc(*this);
    calcBbox();
    updateShaderAttributes();
  };

  virtual void draw(Mat4 &mProject, Mat4 &mModelView) {
    bind();
    Mat4 matrix = this->matrix();
    Mat4 mvp = mProject * mModelView * matrix;
    this->displayShader().setUniform("MVP", mvp);
    updateShaderAttributes();

    this->displayShader().drawIndexed(GL_LINES, 0, mIndices.cols());
  }

}; // LineBuffer;

class DebugBuffer;
using DebugBufferPtr = std::shared_ptr<DebugBuffer>;

class DebugBuffer : public Drawable {
public:
  static DebugBufferPtr create() { return std::make_shared<DebugBuffer>(); }

  DebugBuffer(){};
  ~DebugBuffer(){
      // free everything here
  };

  virtual void init() {
    if (!_lines) {
      _lines = gg::LineBuffer::create();
      _lines->initBuffer();
    }
  }

  void renderLines() {
    using namespace nanogui;

    int numVerts = 0, numIndices = 0;
    for (int i = 0; i < mLineIndices.size(); i++) {
      numIndices += mLineIndices[i].cols();
    }

    for (int i = 0; i < mLinePositions.size(); i++) {
      numVerts += mLinePositions[i].cols();
    }
    _lines->fillBuffer([&](gg::BufferObject &o) -> void {
      o.allocateVerts(numIndices, numVerts);

      auto &indices = o.indices();
      auto &positions = o.positions();
      auto &colors = o.colors();

      int iI = 0;
      int iP = 0;
      for (int i = 0; i < mLineIndices.size(); i++) {
        for (int j = 0; j < mLineIndices[i].cols(); j++) {

          indices.col(iI + j)
              << (Vec2i(iP, iP).cast<unsigned int>() + mLineIndices[i].col(j));
        }
        iP += mLinePositions[i].cols();
        iI += mLineIndices[i].cols();
      }

      iP = 0;
      for (int i = 0; i < mLinePositions.size(); i++) {
        for (int j = 0; j < mLinePositions[i].cols(); j++) {
          positions.col(iP + j) << mLinePositions[i].col(j);
          colors.col(iP + j) = mLineColors[i](Eigen::seq(0, 2), 0);
        }
        iP += mLinePositions[i].cols();
      }
    });
  };

  void pushBox(Vec4 cen4, Vec4 h, Vec4 col) {
    nanogui::MatrixXu indices = nanogui::MatrixXu(2, 12);
    nanogui::MatrixXf positions = nanogui::MatrixXf(3, 8);
    Vec3 cen = cen4(Eigen::seq(0, 2), 0);
    indices.col(0) << 0, 1;
    indices.col(1) << 1, 2;
    indices.col(2) << 2, 3;
    indices.col(3) << 3, 0;

    indices.col(4) << 0, 5;
    indices.col(5) << 1, 6;
    indices.col(6) << 2, 7;
    indices.col(7) << 3, 4;

    indices.col(8) << 4, 5;
    indices.col(9) << 5, 6;
    indices.col(10) << 6, 7;
    indices.col(11) << 7, 4;

    positions.col(0) << cen + Vec3(h[0], h[1], h[2]);
    positions.col(1) << cen + Vec3(-h[0], h[1], h[2]);
    positions.col(2) << cen + Vec3(-h[0], -h[1], h[2]);
    positions.col(3) << cen + Vec3(h[0], -h[1], h[2]);

    positions.col(4) << cen + Vec3(h[0], -h[1], -h[2]);
    positions.col(5) << cen + Vec3(h[0], h[1], -h[2]);
    positions.col(6) << cen + Vec3(-h[0], h[1], -h[2]);
    positions.col(7) << cen + Vec3(-h[0], -h[1], -h[2]);

    mLinePositions.push_back(positions);
    mLineIndices.push_back(indices);
    mLineColors.push_back(col);
  };

  void pushLine(Vec4 c0, Vec4 c1, Vec4 col) {
    nanogui::MatrixXu indices = nanogui::MatrixXu(2, 1);
    nanogui::MatrixXf positions = nanogui::MatrixXf(3, 2);

    indices.col(0) << 0, 1;

    positions.col(0) << c0(Eigen::seq(0, 2), 0);
    positions.col(1) << c1(Eigen::seq(0, 2), 0);
    mLinePositions.push_back(positions);
    mLineIndices.push_back(indices);
    mLineColors.push_back(col);
  };

  virtual void draw(Mat4 &mProject, Mat4 &mModelView) {
    _lines->draw(mProject, mModelView);
  }

  void clear() {
    //_lines->clear();
    mLineColors.clear();
    mLinePositions.clear();
    mLineIndices.clear();
  }

private:
  LineBufferPtr _lines;

  std::vector<Vec4> mLineColors;
  std::vector<nanogui::MatrixXf> mLinePositions;
  std::vector<nanogui::MatrixXu> mLineIndices;
};

BufferObject *makeCube() {
  using namespace nanogui;

  BufferObject *cube = new BufferObject();
  cube->initBuffer();
  cube->fillBuffer([](BufferObject &buf) -> void {
    buf.allocateVerts(12, 8);
    auto indices = buf.indices();
    auto positions = buf.positions();

    indices.col(0) << 0, 1, 2;
    indices.col(1) << 2, 3, 0;
    indices.col(2) << 0, 3, 4;
    indices.col(3) << 4, 5, 0;
    indices.col(4) << 0, 5, 6;
    indices.col(5) << 6, 1, 0;

    indices.col(6) << 1, 6, 7;
    indices.col(7) << 7, 2, 1;

    indices.col(8) << 7, 4, 3;
    indices.col(9) << 3, 2, 7;

    indices.col(10) << 4, 7, 6;
    indices.col(11) << 6, 5, 4;

    positions.col(0) << 1, 1, 1;
    positions.col(1) << -1, 1, 1;
    positions.col(2) << -1, -1, 1;
    positions.col(3) << 1, -1, 1;

    positions.col(4) << 1, -1, -1;
    positions.col(5) << 1, 1, -1;
    positions.col(6) << -1, 1, -1;
    positions.col(7) << -1, -1, -1;
  });
  cube->initDisplayShader();
  // cube->initPickingShader();

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
/*
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
*/
} // namespace gg
#endif
