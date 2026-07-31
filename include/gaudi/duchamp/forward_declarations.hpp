#ifndef __GAUDI_FORWARD_DECLARATIONS__
#define __GAUDI_FORWARD_DECLARATIONS__

#include <memory>

// Forward declarations for commonly used types
// This reduces compilation dependencies when only pointers are needed

namespace hepworth {
    namespace block {
        class projection_solver;
        class projection_constraint;
        class vec3_block;
        class quat_block;
        class sim_block;
        class pinned;
        class stretch_shear;
        class bend;
        class twist;
        class rod_collision;
    }
}

namespace asawa {
    namespace rod {
        class rod;
        class dynamic;
        typedef std::shared_ptr<rod> rod_ptr;
        typedef std::shared_ptr<dynamic> dynamic_ptr;
    }
}

namespace sdf {
    class sdf_base;
    class sdf_sphere;
    class sdf_multi_sphere;
    typedef std::shared_ptr<sdf_base> ptr;
}

namespace gaudi {
    namespace geometry_logger {
        // Forward declarations for geometry logging
    }
}

#endif 