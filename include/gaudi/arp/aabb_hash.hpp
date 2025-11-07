#include "gaudi/arp/hasn_tree.hpp"


// S = stride
template <int S> struct aabb_tree {

  public:
    typedef std::shared_ptr<aabb_tree<S>> ptr;
    typedef aabb_node<S> node;
  
    static ptr create(const std::vector<index_t> &indices,
                      const std::vector<vec3> &vertices, int lvl = 8) {
      return std::make_shared<aabb_tree<S>>(indices, vertices, lvl);
    }
  
    aabb_tree() {}
    ~aabb_tree() {}
  
    aabb_tree(const aabb_tree &other) { *this = other; }
  
    aabb_tree(const std::vector<index_t> &indices,
              const std::vector<vec3> &vertices, int lvl = 8)
        : __x(vertices), __indices(indices) {
      this->build(__indices, __x, lvl);
    }
  
    aabb_tree &operator=(const aabb_tree &rhs) {
      if (this != &rhs) {
        //update_thse
        nodes = rhs.nodes;
        leafNodes = rhs.leafNodes;
        permutation = rhs.permutation;
      }
      return *this;
    }
  
    void build(const std::vector<index_t> &indices,
               const std::vector<vec3> &vertices, int maxLevel) {

    }
  
    const std::vector<index_t> &indices() const { return this->__indices; }
    const std::vector<vec3> &verts() const { return this->__x; }
    const vec3 &vert(index_t i) const { return __x[__indices[i]]; }
  
    void debug() {
    }
  
    void debug_half() {
    }

    const std::vector<vec3> &data,
    const std::vector<index_t> &indices,
    const std::vector<radix_tree_node> &internal_nodes,
    const std::vector<radix_tree_node> &leaf_nodes,
    const TreeResult<ext::extents_t> &bvh_result
  };