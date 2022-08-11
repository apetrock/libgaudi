#ifndef __APP_VIEWER__
#define __APP_VIEWER__

#include <GaudiMath/typedefs.hpp>
#include <nanogui/glutil.h>

//#include <nanogui/glutil.h>
#include <nanogui/nanogui.h>

#include <iostream>
#include <memory>
#include <string>

namespace gg {
using namespace GaudiMath;

class Viewer;
using ViewerPtr = std::shared_ptr<Viewer>;

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

class Scene {
public:
  static ScenePtr create() { return std::make_shared<Scene>(); }

  Scene() {
    // init();
  }
  virtual void onAnimate(int frame){};
  virtual void save(int frame){};
  virtual void onDraw(gg::Viewer &viewer){};

  virtual void _save() { this->save(_frame); };
  virtual void _onAnimate() {
    _frame++;
    this->onAnimate(_frame);
  };

  int _frame = 0;
};

class Viewer {
public:
  static ViewerPtr create(Vec2i sz) { return std::make_shared<Viewer>(sz); }

  Viewer() {
    // init();
  }

  Viewer(Vec2i sz) { init(sz); }

  ~Viewer() {}

  Mat4 makeProjectionMatrix(Real fovY, Real aspectRatio, Real zNear,
                            Real zFar) {
    std::cout << "aspect ration: " << aspectRatio << std::endl;
    Real yScale = tan(M_PI_2 - fovY / 2.0);
    Real xScale = yScale / aspectRatio;
    Mat4 M;
    M.col(0) << xScale, 0, 0, 0;
    M.col(1) << 0, yScale, 0, 0;
    M.col(2) << 0, 0, -(zFar + zNear) / (zFar - zNear), -1;
    M.col(3) << 0, 0, -2 * zNear * zFar / (zFar - zNear), 0;
    return M;
    // return nanogui::frustum(-.1,.1,-.1,.1,.1,20);
  }

  void init(Vec2i size) {

    mProject = this->makeProjectionMatrix(
        60.0 * M_PI / 180.0, double(size[0]) / double(size[1]), 0.1, 20);
    mModelView.setIdentity();
    mModelViewOld.setIdentity();

    mModelRotOld.setIdentity();
    mModelRotNew.setIdentity();

    mDragging = false;

    ball = new nanogui::Arcball();
    ball->setSize(size);
    mPosition = Vec4(0, 0, 4, 0);
  }

  void updateFrame() {
    // cache the old rotation
    mModelRotOld = mModelRotNew;
    ball->setState(Quat::Identity()); // reset ball to identity
  }

  void updatePosition() {
    // This is my camera model, its not a seperat object, but doesn't
    // necessarily warrant one.  All the fancy stuff was implemented in the
    // arcball class
    Mat4 r = ball->matrix();
    // accumulate rotations... this likely could be done with the
    // quaternions in the arcball, but this works after some fiddling
    mModelRotNew = r * mModelRotOld;
    // set the modelView equal to the new rotation
    mModelView = mModelRotNew;
    mModelView.col(3) << Vec4(0, 0, -3, 1);
    // rotate the world, then take a step back
    // get the position from the inverse  of the matrix
    mPosition = mModelRotNew.transpose() * Vec4(0, 0, 4, 1);
  }

  virtual bool onMouseButton(const Eigen::Vector2i &p, int button, bool down,
                             int modifiers) {
    updateFrame();
    ball->button(p, down);
    mDragging = false;

    return true;
  }

  virtual bool onMouseMotion(const Eigen::Vector2i &p,
                             const Eigen::Vector2i &rel, int button,
                             int modifiers) {

    ball->motion(p);
    updatePosition();
    return true;
  }

  Mat4 &getModelView() { return mModelView; }
  Mat4 &getProjection() { return mProject; }

protected:
  // camera variables, this could be in its own class
  Vec4 mPosition;
  Mat4 mProject;
  Mat4 mModelView;
  Mat4 mModelRotNew;
  Mat4 mModelRotOld;
  Mat4 mModelViewOld;

  // variables for selection and dragging, I maintain two selection groups
  // one for widgets and one for objects
  bool mDragging;
  Vec4 mObjCenCache;
  Vec4 mDragStart;

  nanogui::Arcball *ball;
}; // viewer

class FrameGrabber;
using FrameGrabberPtr = std::shared_ptr<FrameGrabber>;

class FrameGrabber {
public:
  static FrameGrabberPtr create(int width, int height) {
    return std::make_shared<FrameGrabber>(width, height);
  }

  FrameGrabber(int width, int height) : _width(width), _height(height) {

    // start ffmpeg telling it to expect raw rgba 720p-60hz frames
    // -i - tells it to read frames from stdin
    const char *cmd = "ffmpeg -r 60 -f rawvideo -pix_fmt rgba -s 1280x720 -i - "
                      "-threads 0 -preset fast -y -pix_fmt yuv420p -crf 21 -vf "
                      "vflip output.mp4";

    // open pipe to ffmpeg's stdin in binary write mode
    ffmpeg = popen(cmd, "w");
    _buffer = new int[_width * _height];
  }

  ~FrameGrabber() { pclose(ffmpeg); }

  void onFrame() {
    std::cout << "grabbing frame" << std::endl;
    glReadPixels(0, 0, _width, _height, GL_RGBA, GL_UNSIGNED_BYTE,
                 this->_buffer);
    fwrite(_buffer, sizeof(int) * _width * _height, 1, this->ffmpeg);
  }

private:
  int _width, _height;
  int *_buffer;
  FILE *ffmpeg;
};

class SimpleApp;
using SimpleAppPtr = std::shared_ptr<SimpleApp>;

class SimpleApp : public nanogui::Screen {
public:
  static SimpleAppPtr create(int w = 1280, int h = 720) {
    return std::make_shared<SimpleApp>(w, h);
  }

  typedef double Real;

  SimpleApp(int w = 1280, int h = 720)
      : _width(w),
        _height(h), nanogui::Screen(Eigen::Vector2i(w, h), "App Simple") {
    using namespace nanogui;

    // now for GUI
    // Window *window = new Window(this, "coordinates");
    // window->setPosition(Vector2i(15, 15));
    // window->setLayout(new GroupLayout());

    performLayout(mNVGContext);

    _viewer = gg::Viewer::create(mSize);

    std::cout << "size: " << mSize << std::endl;

    glEnable(GL_DEPTH_TEST);
    // Accept fragment if it closer to the camera than the former one
    glDepthFunc(GL_LESS);
    // Cull triangles which normal is not towards the camera
    glEnable(GL_CULL_FACE);
    _frameGrabber = FrameGrabber::create(_width, _height);
  }

  ~SimpleApp() {}

  void setScene(ScenePtr scene) { mScene = scene; }
  ///////////////////////////
  //  events
  ///////////////////////////

  virtual bool mouseButtonEvent(const Eigen::Vector2i &p, int button, bool down,
                                int modifiers) {
    _viewer->onMouseButton(p, button, down, modifiers);
    Screen::mouseButtonEvent(p, button, down, modifiers);
    return true;
  }

  virtual bool mouseMotionEvent(const Eigen::Vector2i &p,
                                const Eigen::Vector2i &rel, int button,
                                int modifiers) {

    _viewer->onMouseMotion(p, rel, button, modifiers);
    Screen::mouseMotionEvent(p, rel, button, modifiers);
    return true;
  }

  virtual bool keyboardEvent(int key, int scancode, int action, int modifiers) {
    if (Screen::keyboardEvent(key, scancode, action, modifiers))
      return true;

    if (key == GLFW_KEY_A && action == GLFW_PRESS) {
      _animate = !_animate;
    }

    if (key == GLFW_KEY_S && action == GLFW_PRESS) {
      this->animate();
      this->drawContents();
    }

    if (key == GLFW_KEY_O && action == GLFW_PRESS) {
      if (mScene)
        mScene->_save();
    }

    if (key == GLFW_KEY_Q && action == GLFW_PRESS) {
      nanogui::leave();
    }
    return false;
  }

  ///////////////////////////
  //  Draw
  ///////////////////////////

  virtual void draw(NVGcontext *ctx) {
    // if (false)
    Screen::draw(ctx);
  }

  virtual void animate() {
    if (mScene) {
      mScene->_onAnimate();
    }
  }

  virtual void drawContents() {
    using namespace nanogui;
    glfwGetTime();

    if (this->_animate) {
      this->animate();
    }

    glEnable(GL_DEPTH_TEST);

    if (mScene)
      mScene->onDraw(*_viewer);

    glDisable(GL_DEPTH_TEST);

    if (this->_animate) {
      _frameGrabber->onFrame();
    }

    _frame++;
  }

private:
  bool _animate = false;
  unsigned int _frame = 0;
  ScenePtr mScene;
  gg::ViewerPtr _viewer;

  int _width, _height;
  FrameGrabberPtr _frameGrabber;

}; // Simple App

} // namespace gg
#endif
