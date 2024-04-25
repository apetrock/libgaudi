#include "GaudiMath/typedefs.hpp"
#include <algorithm>
#include <cmath>
#include <exception>
#include <iterator>
#include <vector>
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

#include <complex>
#include <iostream>
#include <random>
#include <string>
#include <bit>

#include "GaudiGraphics/buffers.hpp"
#include "GaudiGraphics/geometry_logger.h"
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"
// #include "gaudi/asawa/asawa.h"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/logger.hpp"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;
using uint = unsigned int;
using uint2 = std::array<uint, 2>;

// Expands a 10-bit integer into 30 bits
// by inserting 2 zeros after each bit.
inline uint expandBits(uint v)
{
  v = (v * 0x00010001u) & 0xFF0000FFu;
  v = (v * 0x00000101u) & 0x0F00F00Fu;
  v = (v * 0x00000011u) & 0xC30C30C3u;
  v = (v * 0x00000005u) & 0x49249249u;
  return v;
}

inline uint scale(float x, float c)
{
  return min(max(x * c, 0.0f), c - 1.0f);
}

// dump binary to a stream to use in a cout ie << dump_binary(i) << endl;
std::string dump_binary(uint i)
{
  std::string s;
  for (int j = 31; j >= 0; j--)
  {
    s += ((i & (1 << j)) ? '1' : '0');
  }
  return s;
}

// Calculates a 30-bit Morton code for the
// given 3D point located within the unit cube [0,1].
inline uint morton3D(float x, float y, float z)
{
  x = scale(x, 1024.0f);
  y = scale(y, 1024.0f);
  z = scale(z, 1024.0f);
  uint xx = expandBits((uint)x);
  uint yy = expandBits((uint)y);
  uint zz = expandBits((uint)z);
  return xx * 4 + yy * 2 + zz;
}

inline int clz(uint i, uint j, const std::vector<int> &ids, const std::vector<int> &hash)
{

  if (j < 0)
    return -1;
  if (j > ids.size() - 1)
    return -1;
  int code_i = hash[ids[i]];
  int code_j = hash[ids[j]];
  // std::cout << "     " << code_i << " " << dump_binary(code_i) << std::endl;
  // std::cout << "     " << code_j << " "<< dump_binary(code_j) << std::endl;
  // std::cout << "     " << __builtin_clz(code_i ^ code_j) << std::endl;
  return __builtin_clz(code_i ^ code_j);
}

uint find_split(int start, int end, const std::vector<int> &ids, const std::vector<int> &hash)
{

  // uint first_code = hash[ids[start]];
  // uint last_code = hash[ids[end]];
  // if (first_code == last_code)
  //   return (first_code + last_code) >> 1;

  // uint common_prefix = std::countl_zero(first^last);
  int common_prefix_dist = clz(start, end, ids, hash);
  int last_prefix_dist = clz(end - 1, end, ids, hash);
  uint split = start;
#if 0
  for (int i = start; i < end; i++){
    uint curr = hash[ids[i]];
    uint new_prefix = __builtin_clz(first ^ curr);
    std::cout << "new_prefix: " << new_prefix << " " << common_prefix << std::endl;
    if (new_prefix > common_prefix)
    {
      split = i;
      //common_prefix = new_prefix;
    }
  }
#endif

#if 1
  uint step = end - start;
  while (step > 1)
  {
    step = (step + 1) >> 1; // exponential decrease
    int new_split = split + step;
    if (new_split < end)
    {

      int new_prefix_dist = clz(start, new_split, ids, hash);

      if (new_prefix_dist > common_prefix_dist)
      {
        split = new_split;
        last_prefix_dist = new_prefix_dist;
      }
    }
  }
#endif

  return split;
}

uint2 find_range(const int &i, const std::vector<int> &ids, const std::vector<int> &hash)
{
  uint N = ids.size();

  int dir = va::sgn(clz(i, i + 1, ids, hash) - clz(i, i - 1, ids, hash));
  int sig_min = clz(i, i - dir, ids, hash);

  int lmax = 2;
  while (clz(i, i + lmax * dir, ids, hash) > sig_min)
    lmax *= 2;

  uint l = 0;
  uint t = lmax;

  while (t >= 1)
  {
    t = t >> 1;
    if (clz(i, i + (l + t) * dir, ids, hash) > sig_min)
      l += t;
  }
  uint j = i + l * dir;

  return dir < 0 ? uint2{j, uint(i)} : uint2{uint(i), j};
}

const uint UNULL = uint(-1);

struct radix_tree_node
{
  uint start = UNULL;
  uint end = UNULL;
  uint split = UNULL; // split is the index of the last element in the left child
  uint parent = UNULL;
};

void dump_nodes(const std::vector<radix_tree_node> &nodes)
{
  int i = 0;
  for (auto &n : nodes)
  {
    std::cout << "node[" << i++ << "]: " << n.start << " " << n.end << " " << n.split << " " << n.parent << std::endl;
  }
}

std::vector<radix_tree_node> build_tree(const std::vector<int> &ids, const std::vector<int> &hash)
{
  std::vector<radix_tree_node> nodes(ids.size() + ids.size() - 1);
  uint leaf_start = ids.size() - 1;
  for (int i = 0; i < ids.size(); i++)
  {
    nodes[leaf_start + i].start = leaf_start + i;
    nodes[leaf_start + i].end = leaf_start + i + 1;
    nodes[leaf_start + i].split = UNULL;
    nodes[leaf_start + i].parent = UNULL;
  }

  for (int i = 0; i < ids.size() - 1; i++)
  {
    uint2 range = find_range(i, ids, hash);
    uint split = find_split(range[0], range[1], ids, hash);

    nodes[i].start = range[0];
    nodes[i].end = range[1];
    nodes[i].split = split;
    if (split == range[0])
      nodes[leaf_start + range[0]].parent = i;
    else
      nodes[split].parent = i;

    if (split + 1 == range[1])
      nodes[leaf_start + range[1]].parent = i;
    else
      nodes[split + 1].parent = i;
  }
#if 0 // dump nodes
  dump_nodes(nodes);
#endif
  return nodes;
}

void test_tree(const std::vector<radix_tree_node> &nodes, const std::vector<int> &ids, const std::vector<int> &hash)
{
  std::vector<bool> visited(ids.size(), false);
  std::stack<uint> stack;
  stack.push(0);
  uint leaf_start = ids.size() - 1;
  bool no_asserts = true;
  while (stack.size() > 0)
  {
    uint i = stack.top();
    stack.pop();
    std::cout << "i: " << i << std::endl;
    // cout current node
    std::cout << "node["
              << i << "]: "
              << nodes[i].start << " "
              << nodes[i].end << " "
              << nodes[i].split << " "
              << nodes[i].parent << " " << std::endl;

    if (nodes[i].split + 0 == nodes[i].start ||
        nodes[i].split + 1 == nodes[i].end)
    {
      uint j0 = nodes[i].start;
      uint j1 = nodes[i].end;
      visited[leaf_start + j0] = true;
      visited[leaf_start + j1] = true;
      assert(nodes[leaf_start + j0].parent == i);
      assert(nodes[leaf_start + j1].parent == i);
      //
    }
    else
    {
      uint start_i = nodes[i].start;
      uint end_i = nodes[i].end;

      uint j0 = nodes[i].split;
      uint j1 = nodes[i].split + 1;
      // assert parents
      assert(nodes[j0].parent == i);
      assert(nodes[j1].parent == i);
      // assert rangings
      assert(nodes[j0].start >= nodes[i].start);
      assert(nodes[j0].end <= nodes[i].end);
      assert(nodes[j1].start >= nodes[i].start);
      assert(nodes[j1].end <= nodes[i].end);
      stack.push(j0);
      stack.push(j1);
    }
  }
  bool visited_all = false;
  for (int i = 0; i < visited.size(); i++)
  {
    if (!visited[i])
    {
      visited_all = false;
      break;
    }
  }
}

void unit_test_tree()
{
  std::vector<int> hash = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  // std::vector<int> hash = {0, 1, 2, 3, 4, 5, 6, 7};

  std::vector<int> ids = hash;
  std::vector<radix_tree_node> nodes = build_tree(ids, hash);
  test_tree(nodes, ids, hash);
}

namespace gaudi
{
  using extents_3 = std::array<vec3, 2>;

  template <int STRIDE>
  ext::extents_t calc_extents(index_t i, const std::vector<index_t> &indices,
                              const std::vector<vec3> &vertices)
  {
    ext::extents_t extents = ext::init();
    extents[0] = vertices[indices[STRIDE * i + 0]];
    extents[1] = extents[0];
    for (int k = 0; k < STRIDE; k++)
    {
      vec3 p = vertices[indices[STRIDE * i + k]];
      extents = ext::expand(extents, p);
    }

    return extents;
  }

  ext::extents_t pyramid(const ext::extents_t &a, const ext::extents_t &b)
  {
    vec3 mn = va::min(a[0], b[0]); // assume atomic min/max on gpu
    vec3 mx = va::max(a[1], b[1]);
    return {mn, mx};
  }

  // function object prototype that takes T A and B and returns T, to generically perform pyramid ops

  template <typename T>
  std::vector<T> build_pyramid(const std::vector<T> &data, int leaf_start, const std::vector<radix_tree_node> &nodes,
                               const std::function<T()> &init,
                               const std::function<T(const T &, const T &)> &op)
  {
    // THE GPU is probably going to need a max tree depth, but we can just do a max recursion rate of 32?
    const int max_depth = 32;
    std::vector<T> pyramid(nodes.size(), init());

    for (int i = 0; i < data.size(); i++)
    {
      // insert leaf then propagate up the pyramid
      int leaf_start = data.size() - 1;
      int ii = leaf_start + i;

      T data_i = data[i];

      pyramid[ii] = T(data_i);
      int parent = nodes[ii].parent;

      int j = 0;
      while (j < max_depth && parent != UNULL)
      {

        pyramid[parent] = op(pyramid[parent], data_i);
        data_i = T(pyramid[parent]);
        parent = nodes[parent].parent;

        j++;
      }
    }
    return pyramid;
  }

  void dump_cube_list()
  {
    for (int i = 0; i < 8; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        int neighbor = i ^ (1 << j); // Flip j-th bit of i
        if (neighbor <= i)
          continue; // Avoid drawing the same line twice
        std::cout << ", " << i << ", " << neighbor;
      }
    }
    std::cout << std::endl;
  }

  const int cube_edges[24] = {0, 1, 0, 2, 0, 4, 1, 3, 1, 5, 2, 3, 2, 6, 3, 7, 4, 5, 4, 6, 5, 7, 6, 7};
  const int cube_verts[8][3] = {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}};
  inline std::array<vec3, 8> create_cube_vertices(const ext::extents_t &extents)
  {
    vec3 xf0 = extents[0];
    vec3 xf1 = extents[1];
    return {vec3(xf0[0], xf0[1], xf0[2]),
            vec3(xf0[0], xf0[1], xf1[2]),
            vec3(xf0[0], xf1[1], xf0[2]),
            vec3(xf0[0], xf1[1], xf1[2]),
            vec3(xf1[0], xf0[1], xf0[2]),
            vec3(xf1[0], xf0[1], xf1[2]),
            vec3(xf1[0], xf1[1], xf0[2]),
            vec3(xf1[0], xf1[1], xf1[2])};
  }

  void log_extents(const ext::extents_t &extents, vec4 color = vec4(1.0, 0.0, 0.0, 1.0))
  {

    std::array<vec3, 8> vertices = create_cube_vertices(extents);
    for (int i = 0; i < 12; i++)
    {
      int v0 = cube_edges[2 * i + 0];
      int v1 = cube_edges[2 * i + 1];
      // logger::line(vertices[v0], vertices[v1], vec4(0.0, 0.0, 1.0, 1.0));
      //  alternatively using the cube_verts array
      auto [i0, j0, k0] = cube_verts[v0];
      auto [i1, j1, k1] = cube_verts[v1];
      logger::line(vec3(extents[i0][0], extents[j0][1], extents[k0][2]),
                   vec3(extents[i1][0], extents[j1][1], extents[k1][2]), color);
      // logger::line(vec3(cube_verts[v0][0], cube_verts[v0][1], cube_verts[v0][2]), vec3(cube_verts[v1][0], cube_verts[v1][1], cube_verts[v1][2]), vec4(0.0, 0.0, 1.0, 1.0));
    }
  }

  void log_extents_to_parent(int id,
                             int leaf_start,
                             const std::vector<radix_tree_node> &nodes, const std::vector<ext::extents_t> &extents)
  {

    int parent = nodes[leaf_start + id].parent;
    int j = 0;
    const int max_depth = 32;
    log_extents(extents[leaf_start + id], vec4(1.0, 1.0, 0.0, 1.0));

    for (int i = nodes[parent].start; i < nodes[parent].end; i++)
    {
      std::cout << id << "  i: " << i << " p: " << parent << std::endl;
      log_extents(extents[leaf_start + i], vec4(0.0, 1.0, 0.0, 1.0));
    }

    while (j < max_depth && parent != UNULL)
    {
      log_extents(extents[parent]);
      parent = nodes[parent].parent;
      j++;
    }
  }

  namespace duchamp
  {
    using namespace asawa;

    class aabb_build
    {
    public:
      typedef std::shared_ptr<aabb_build> ptr;

      static ptr create() { return std::make_shared<aabb_build>(); }

      aabb_build() { load_shell(); };

      void load_shell()
      {
        //__M = load_cube();
        __M = shell::load_bunny();
        //__M = shell::load_crab();

        shell::triangulate(*__M);

        std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
        asawa::center(x, 2.0);

        _eps = asawa::shell::avg_length(*__M, x);
      }

      void step_dynamics(int frame)
      {
        unit_test_tree();
        asawa::shell::shell &M = *__M;

        std::vector<vec3> &x = asawa::get_vec_data(M, 0);
        std::vector<vec3> cens = asawa::shell::face_centers(M, x);

        std::vector<index_t> face_vert_ids = M.get_face_vert_ids();

        // calc bounding box of cens
        vec3 min = cens[0];
        vec3 max = cens[0];
        for (int i = 0; i < cens.size(); i++)
        {
          min = va::min(min, cens[i]);
          max = va::max(max, cens[i]);
        }

        std::vector<int> hash(cens.size(), 0);
        std::vector<int> indices(cens.size(), 0);
        for (int i = 0; i < cens.size(); i++)
        {
          indices[i] = i;
          vec3 c = cens[i] - min;
          vec3 n = max - min;
          c = va::div(c, n);
          hash[i] = morton3D(c[0], c[1], c[2]);
        }

        std::sort(indices.begin(), indices.end(), [&](int a, int b)
                  { return hash[a] < hash[b]; });

        std::vector<radix_tree_node> nodes = build_tree(indices, hash);
        int r0 = nodes[nodes[nodes[0].split].split].start;
        int r1 = nodes[nodes[nodes[0].split].split].end;
        int leaf_start = indices.size() - 1;

        std::vector<extents_3> extents(face_vert_ids.size() / 3);
        for (int i = 0; i < face_vert_ids.size() / 3; i++)
        {
          extents[i] = calc_extents<3>(indices[i], face_vert_ids, x);
          //log_extents(extents[i]);
        }

        std::vector<extents_3> ext = build_pyramid<extents_3>(
            extents, leaf_start, nodes,
            []()
            { return ext::init(); },
            [](const extents_3 &a, const extents_3 &b)
            { return pyramid(a, b); });

        int test_id = leaf_start + 1;
        log_extents_to_parent(1000, leaf_start, nodes, ext);
        log_extents_to_parent(1, leaf_start, nodes, ext);
        log_extents_to_parent(2, leaf_start, nodes, ext);
        
        /*
        for (int i = 0; i < (indices.size() - 1); i++)
        {
          vec3 xf0 = cens[indices[i + 0]];
          vec3 xf1 = cens[indices[i + 1]];
          if (i > r0 && i < r1)
            logger::line(xf0, xf1, vec4(0.0, 1.0, 0.0, 1.0));
          else
            logger::line(xf0, xf1, vec4(1.0, 0.0, 0.0, 1.0));
        }
        */
      }

      void step(int frame)
      {
        std::cout << "frame: " << frame << std::endl;
        std::cout << "  -surface" << std::endl;
        std::cout << "    -corners: " << __M->corner_count() << std::endl;
        std::cout << "    -verts: " << __M->vert_count() << std::endl;
        std::cout << "    -faces: " << __M->face_count() << std::endl;

        step_dynamics(frame);
      }

      real _eps;
      shell::shell::ptr __M;
    };

  } // namespace duchamp
} // namespace gaudi

using namespace GaudiMath;

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

class Scene : public gg::Scene
{

public:
  static ScenePtr create() { return std::make_shared<Scene>(); }

  Scene() : gg::Scene() { initScene(); }

  void initScene()
  {
    //_experiment = duchamp::mean_shift_experiment<growth>::create();

    _objs.resize(2);

    _objs[0] = gg::BufferObject::create();
    _objs[0]->init();
    mSceneObjects.push_back(_objs[0]);

    __surf = gaudi::duchamp::aabb_build::create();
    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);

    colors = {
        gg::colorRGB(0.5, 0.5, 0.5, 1.0),
        gg::colorRGB(0.0, 1.0, 1.0, 1.0),
    };
  }
  virtual void onAnimate(int frame)
  {

    __surf->step(frame);

    gg::fillBuffer_ref(*__surf->__M, _objs[0], colors[0]);

    //   std::cout << "rendering debug" << std::endl;
    //   asawa::test();

    gg::geometry_logger::render();
  }

  virtual void onDraw(gg::Viewer &viewer)
  {

    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable
                  {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });

    gg::geometry_logger::clear();
  }

private:
  // gaudi::duchamp::fast_summation_test::ptr __surf;
  gaudi::duchamp::aabb_build::ptr __surf;

  std::vector<gg::DrawablePtr> mSceneObjects;
  std::vector<gg::BufferObjectPtr> _objs;
  vector<gg::colorRGB> colors;
};

std::string GetCurrentWorkingDir(void)
{
  char buff[FILENAME_MAX];
  GetCurrentDir(buff, FILENAME_MAX);
  std::string current_working_dir(buff);
  return current_working_dir;
}

class App;
using AppPtr = std::shared_ptr<App>;

class App : public gg::SimpleApp
{
public:
  static AppPtr create(std::string file) { return std::make_shared<App>(file); }

  typedef double Real;

  App(std::string file) : gg::SimpleApp(1280, 720, 4.0, true, "vortex_test_")
  {
    this->setScene(scene = Scene::create());
    this->initUI();
  }

  void initUI()
  {
    using namespace nanogui;
    int w = 256;
    performLayout();
    // window->center();
  }

  ~App() {}

  ScenePtr scene;
};

int main(int argc, char *argv[])
{
  try
  {
    cout << "You have entered " << argc << " arguments:"
         << "\n";

    for (int i = 0; i < argc; ++i)
      cout << argv[i] << "\n";

    nanogui::init();

    AppPtr app = App::create(std::string(argv[0]));

    // app->setScene(Scene::create());
    app->drawAll();
    app->setVisible(true);
    nanogui::mainloop();
    // delete app;
    nanogui::shutdown();
  }
  catch (const std::runtime_error &e)
  {
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
