#ifndef __APP_VIEWER__
#define __APP_VIEWER__

#include <GaudiMath/typedefs.hpp>
#include <nanogui/glutil.h>

//#include <nanogui/glutil.h>
#include <nanogui/nanogui.h>

#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

inline void GLAPIENTRY MessageCallback(GLenum source, GLenum type, GLuint id,
                                       GLenum severity, GLsizei length,
                                       const GLchar *message,
                                       const void *userParam) {
  fprintf(stderr,
          "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
          (type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : ""), type, severity,
          message);
}

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
  static ViewerPtr create(Vec2i sz, double d = 4.0) {
    return std::make_shared<Viewer>(sz, d);
  }

  Viewer() {
    // init();
  }

  Viewer(Vec2i sz, double d = 4.0) : mDist(d) { init(sz); }

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

    mPosition = Vec4(0, 0, mDist, 1);
    this->updatePosition();
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

    mModelView.col(3) << Vec4(0, 0, -mDist, 1);
    // rotate the world, then take a step back
    // get the position from the inverse  of the matrix
    mPosition = mModelRotNew.transpose() * Vec4(0, 0, mDist, 1);
  }

  void rotate_ball() {
    Vec3 rotation(0.0, M_PI / 360.0, 0.0);
    double angle = rotation.norm();
    Vec3 axis = rotation.normalized();
    Eigen::Quaternionf q(Eigen::AngleAxisf(angle, axis));
    ball->state() *= q;
  }

  virtual bool onMouseButton(const Eigen::Vector2i &p, int button, bool down,
                             int modifiers) {
    if (down) {
      mode = button;
      pLast = p;

    } else
      mode = -1;

    updateFrame();
    ball->button(p, down);
    mDragging = false;

    return true;
  }

  virtual bool onMouseMotion(const Eigen::Vector2i &p,
                             const Eigen::Vector2i &rel, int button,
                             int modifiers) {
    if (mode == 2) {
      int delt = p[1] - pLast[1];
      mDist += 0.01 * float(delt);
    } else
      ball->motion(p);
    updatePosition();
    pLast = p;
    return true;
  }

  void onFrame() {
    rotate_ball();
    // updatePosition();
  }

  Mat4 &getModelView() { return mModelView; }
  Mat4 &getProjection() { return mProject; }
  Vec4 &getPosition() { return mPosition; }
  Vec3 getPosition3() {
    return Vec3(mPosition[0], //
                mPosition[1], //
                mPosition[2]);
  }

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
  double mDist = 3.0;
  bool mDragging;
  Eigen::Vector2i pLast;
  Vec4 mObjCenCache;
  Vec4 mDragStart;
  int mode = -1;
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
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glReadPixels(0, 0, _width, _height, GL_RGBA, GL_UNSIGNED_BYTE,
                 this->_buffer);
    fwrite(_buffer, sizeof(int) * _width * _height, 1, this->ffmpeg);
  }

private:
  int _width, _height;
  int *_buffer;
  FILE *ffmpeg;
};

class RenderingEffect;
using RenderingEffectPtr = std::shared_ptr<RenderingEffect>;

class RenderingEffect {
public:
  virtual ~RenderingEffect() = default;

  RenderingEffect(int w = 1280, int h = 720) : _width(w), _height(h) {}

  virtual void initFbo(){};
  virtual void initShader(){};
  virtual unsigned int getFBO() { return 0; };
  virtual void render(ViewerPtr view){};

  void setIndinces() {
    Eigen::MatrixXi indices(3, 2);
    indices.col(0) << 0, 1, 2;
    indices.col(1) << 2, 1, 3;
    Eigen::MatrixXf positions(3, 4);
    positions.col(0) << -1, -1, 0;
    positions.col(1) << 1, -1, 0;
    positions.col(2) << -1, 1, 0;
    positions.col(3) << 1, 1, 0;
    // bind the shader and upload vertex positions and indices
    mShader.bind();
    mShader.uploadIndices(indices);
    mShader.uploadAttrib("aPos", positions);
  }

  void draw() {
    setIndinces();
    mShader.drawIndexed(GL_TRIANGLES, 0, 2);
  }

  nanogui::GLShader mShader;

  Real _width, _height;
};

class HdrEffect;
using HdrEffectPtr = std::shared_ptr<HdrEffect>;

class HdrEffect : public RenderingEffect {
public:
  static HdrEffectPtr create(int w = 1280, int h = 720) {
    return std::make_shared<HdrEffect>(w, h);
  }

  HdrEffect(int w = 1280, int h = 720) : RenderingEffect(w, h) {
    initFbo();
    initShader();
    std::cout << "HDR effect created" << std::endl;
    std::cout << _width << " " << _height << std::endl;
  }

  ~HdrEffect() {}

  virtual void initFbo() override {

    // configure floating point framebuffer
    // ------------------------------------

    glGenFramebuffers(1, &mFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, mFBO);
    // create floating point color buffer

    glGenTextures(1, &mColorBuffer);
    glBindTexture(GL_TEXTURE_2D, mColorBuffer);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, _width, _height, 0, GL_RGBA,
                 GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                           mColorBuffer, 0);
    // create depth buffer (renderbuffer)

    glGenRenderbuffers(1, &mRboDepth);
    glBindRenderbuffer(GL_RENDERBUFFER, mRboDepth);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, _width, _height);
    // attach buffers

    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, mRboDepth);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      std::cout << "Framebuffer not complete!" << std::endl;
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

  virtual void initShader() override {
    this->mShader.bind();
    this->mShader.init("hdr_shader", get_shader("sqr_vert"),
                       get_shader("hdr_sqr_frag"));
  }

  virtual unsigned int getFBO() override { return mFBO; }

  virtual void render(ViewerPtr view) override {

    mShader.bind();

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, mColorBuffer);
    this->mShader.setUniform("hdrBuffer", 0);
    this->mShader.setUniform("hdr", false);
    this->mShader.setUniform("exposure", 1.0);
    this->draw();
  }

  unsigned int mFBO;
  unsigned int mColorBuffer;
  unsigned int mRboDepth;
};

class DeferredShadingEffect;
using DeferredShadingEffectPtr = std::shared_ptr<DeferredShadingEffect>;

class DeferredShadingEffect : public RenderingEffect {
public:
  static DeferredShadingEffectPtr create(int w = 1280, int h = 720) {
    return std::make_shared<DeferredShadingEffect>(w, h);
  }

  DeferredShadingEffect(int w = 1280, int h = 720) : RenderingEffect(w, h) {
    initFbo();
    initShader();
    std::cout << "Deferred Shading effect created" << std::endl;
    std::cout << _width << " " << _height << std::endl;
  }

  ~DeferredShadingEffect() {}

  virtual void initFbo() override {

    // configure floating point framebuffer
    // ------------------------------------

    // configure g-buffer framebuffer
    // ------------------------------
    glGenFramebuffers(1, &gBuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, gBuffer);

    // position color buffer
    glGenTextures(1, &gPosition);
    glBindTexture(GL_TEXTURE_2D, gPosition);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, _width, _height, 0, GL_RGBA,
                 GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                           gPosition, 0);
    // normal color buffer
    glGenTextures(1, &gNormal);
    glBindTexture(GL_TEXTURE_2D, gNormal);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, _width, _height, 0, GL_RGBA,
                 GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D,
                           gNormal, 0);
    // color + specular color buffer
    glGenTextures(1, &gAlbedoSpec);
    glBindTexture(GL_TEXTURE_2D, gAlbedoSpec);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, _width, _height, 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D,
                           gAlbedoSpec, 0);
    // tell OpenGL which color attachments we'll use (of this framebuffer) for
    // rendering
    unsigned int attachments[3] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1,
                                   GL_COLOR_ATTACHMENT2};
    glDrawBuffers(3, attachments);
    // create and attach depth buffer (renderbuffer)

    glGenRenderbuffers(1, &mRboDepth);
    glBindRenderbuffer(GL_RENDERBUFFER, mRboDepth);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, _width, _height);

    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, mRboDepth);
    // finally check if framebuffer is complete
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      std::cout << "Framebuffer not complete!" << std::endl;
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

  virtual void initShader() override {
    this->mShader.bind();
    this->mShader.init("deffered_shader", get_shader("sqr_vert"),
                       get_shader("def_sqr_frag"));
    // this->mShader.init("hdr_shader", get_shader("hdr_sqr_vert"),
    //                    get_shader("hdr_sqr_frag"));
  }

  virtual unsigned int getFBO() override { return gBuffer; }

  virtual unsigned int getPositionTexture() { return gPosition; }
  virtual unsigned int getNormalTexture() { return gNormal; }
  virtual void setSsaoTexture(unsigned int tex) { ssaoTexture = tex; }

  virtual void render(ViewerPtr view) override {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    mShader.bind();

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, gPosition);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, gNormal);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, gAlbedoSpec);
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, ssaoTexture);

    std::cout << "gs: " << gBuffer << " " << gPosition << " " << gNormal << " "
              << gAlbedoSpec << std::endl;
    // this->mShader.setUniform("hdrBuffer", 2);
    this->mShader.setUniform("gPosition", 0);
    this->mShader.setUniform("gNormal", 1);
    this->mShader.setUniform("gAlbedoSpec", 2);
    this->mShader.setUniform("ssao", 3);

    this->mShader.setUniform("viewPos", view->getPosition3());
    this->draw();

    glBindFramebuffer(GL_READ_FRAMEBUFFER, gBuffer);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    glBlitFramebuffer(0, 0, _width, _height, 0, 0, _width, _height,
                      GL_DEPTH_BUFFER_BIT, GL_NEAREST);
  }

  unsigned int gBuffer;
  unsigned int gPosition, gNormal, gAlbedoSpec, ssaoTexture;
  unsigned int mRboDepth;
};

class SsaoShadingEffect;
using SsaoShadingEffectPtr = std::shared_ptr<SsaoShadingEffect>;

class SsaoShadingEffect : public RenderingEffect {
public:
  static SsaoShadingEffectPtr create(int w = 1280, int h = 720) {
    return std::make_shared<SsaoShadingEffect>(w, h);
  }

  SsaoShadingEffect(int w = 1280, int h = 720) : RenderingEffect(w, h) {
    initFbo();
    initKernel();
    initNoise();
    initShader();
    std::cout << "Deferred Shading effect created" << std::endl;
    std::cout << _width << " " << _height << std::endl;
  }

  ~SsaoShadingEffect() {}
  void initKernel() {
    // generate sample kernel
    // ----------------------
    std::uniform_real_distribution<GLfloat> randomFloats(
        0.0, 1.0); // generates random floats between 0.0 and 1.0
    std::default_random_engine generator;

    for (unsigned int i = 0; i < 64; ++i) {
      Vec3 sample(randomFloats(generator) * 2.0 - 1.0, //
                  randomFloats(generator) * 2.0 - 1.0, //
                  randomFloats(generator));
      sample.normalize();

      sample *= randomFloats(generator);
      Real scale = Real(i) / 64.0;

      // scale samples s.t. they're more aligned to center of kernel
      Real f = 0.1;
      scale = (1.0 - f) * 1.0 + f * scale * scale;
      sample *= scale;
      ssaoKernel.push_back(sample);
    }
  }

  void initNoise() {
    // generate noise texture
    // ----------------------
    std::uniform_real_distribution<GLfloat> randomFloats(
        0.0, 1.0); // generates random floats between 0.0 and 1.0
    std::default_random_engine generator;

    std::vector<Vec3> ssaoNoise;
    for (unsigned int i = 0; i < 16; i++) {
      Vec3 noise(randomFloats(generator) * 2.0 - 1.0,
                 randomFloats(generator) * 2.0 - 1.0,
                 0.0f); // rotate around z-axis (in tangent space)
      ssaoNoise.push_back(noise);
    }

    glGenTextures(1, &noiseTexture);
    glBindTexture(GL_TEXTURE_2D, noiseTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, 4, 4, 0, GL_RGB, GL_FLOAT,
                 &ssaoNoise[0]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  }

  virtual void initFbo() override {

    // also create framebuffer to hold SSAO processing stage
    // -----------------------------------------------------

    glGenFramebuffers(1, &ssaoFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, ssaoFBO);
    // SSAO color buffer
    glGenTextures(1, &ssaoColorBuffer);
    glBindTexture(GL_TEXTURE_2D, ssaoColorBuffer);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, _width, _height, 0, GL_RED, GL_FLOAT,
                 NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                           ssaoColorBuffer, 0);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      std::cout << "SSAO Framebuffer not complete!" << std::endl;

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

  virtual void initShader() override {
    this->mShader.bind();
    this->mShader.init("deffered_shader", get_shader("sqr_vert"),
                       get_shader("ssao_sqr_frag"));
    // this->mShader.init("hdr_shader", get_shader("hdr_sqr_vert"),
    //                    get_shader("hdr_sqr_frag"));
  }

  virtual unsigned int getFBO() override { return ssaoFBO; }

  virtual void setPositionTexture(unsigned int tex) { gPosition = tex; }
  virtual void setNormalTexture(unsigned int tex) { gNormal = tex; }
  virtual unsigned int getSsaoTexture() { return ssaoColorBuffer; }

  virtual void render(ViewerPtr view) override {

    // 2. generate SSAO texture
    // ------------------------
    glBindFramebuffer(GL_FRAMEBUFFER, ssaoFBO);
    glClear(GL_COLOR_BUFFER_BIT);

    mShader.bind();
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, gPosition);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, gNormal);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, noiseTexture);

    for (unsigned int i = 0; i < 64; ++i)
      this->mShader.setUniform("samples[" + std::to_string(i) + "]",
                               ssaoKernel[i]);
    std::cout << view->getProjection() << std::endl;
    this->mShader.setUniform("projection", view->getProjection());
    this->mShader.setUniform("gPosition", 0);
    this->mShader.setUniform("gNormal", 1);
    this->mShader.setUniform("texNoise", 2);
    this->draw();
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }
  unsigned int gPosition, gNormal;
  unsigned int noiseTexture;
  std::vector<Vec3> ssaoKernel;
  unsigned int ssaoFBO;
  unsigned int ssaoColorBuffer;

  // nanogui::GLShader mBlurShader;
};

class SimpleApp;
using SimpleAppPtr = std::shared_ptr<SimpleApp>;

class SimpleApp : public nanogui::Screen {
public:
  static SimpleAppPtr create(int w = 1280, int h = 720) {
    return std::make_shared<SimpleApp>(w, h);
  }

  typedef double Real;

  SimpleApp(int w = 1280, int h = 720, double d = 4.0, bool grab = true)
      : _width(w),
        _height(h), nanogui::Screen(Eigen::Vector2i(w, h), "App Simple") {
    using namespace nanogui;

    // now for GUI
    // Window *window = new Window(this, "coordinates");
    // window->setPosition(Vector2i(15, 15));
    // window->setLayout(new GroupLayout());

    performLayout(mNVGContext);

    _viewer = gg::Viewer::create(mSize, d);

    std::cout << "size: " << mSize << std::endl;
    // During init, enable debug output
    glEnable(GL_DEBUG_OUTPUT);
    glDebugMessageCallback(MessageCallback, 0);

    glEnable(GL_DEPTH_TEST);
    // Accept fragment if it closer to the camera than the former one
    glDepthFunc(GL_LESS);
    // Cull triangles which normal is not towards the camera
    glEnable(GL_CULL_FACE);
    if (grab)
      _frameGrabber = FrameGrabber::create(_width, _height);

    _renderingEffects.resize(2);
    _renderingEffects[0] = SsaoShadingEffect::create(_width, _height);
    _renderingEffects[1] = DeferredShadingEffect::create(_width, _height);
    SsaoShadingEffectPtr ssao =
        std::dynamic_pointer_cast<SsaoShadingEffect>(_renderingEffects[0]);
    DeferredShadingEffectPtr deffered =
        std::dynamic_pointer_cast<DeferredShadingEffect>(_renderingEffects[1]);
    ssao->setPositionTexture(deffered->getPositionTexture());
    ssao->setNormalTexture(deffered->getNormalTexture());
    deffered->setSsaoTexture(ssao->getSsaoTexture());
    //_renderingEffect = HdrEffect::create(_width, _height);
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
    _viewer->onFrame();
  }

  virtual void drawContents() {
    using namespace nanogui;
    glfwGetTime();

    if (this->_animate) {
      this->animate();
    }

    glEnable(GL_DEPTH_TEST);
    glBindFramebuffer(GL_FRAMEBUFFER, _renderingEffects[1]->getFBO());
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    if (mScene)
      mScene->onDraw(*_viewer);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    for (auto &effect : _renderingEffects)
      effect->render(_viewer);

    glDisable(GL_DEPTH_TEST);

    if (this->_animate && _frameGrabber) {
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
  std::vector<RenderingEffectPtr> _renderingEffects;

}; // Simple App

} // namespace gg
#endif
