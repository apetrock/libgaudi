#ifndef __HENSON_RIGID_BODY__
#define __HENSON_RIGID_BODY__

//named after Jim Henson, creator of the Muppets
//basis for linkage and rigid body framework
//using projection based dynamics
#include "gaudi/common.h"
#include "gaudi/asawa/graph/graph.hpp"

namespace gaudi {
namespace henson {

enum class joint_type { fixed, hinge, slider, ball, universal, planar, free, custom };

//joint constraints project the relative position of two rigid bodies
// -given current orientation of body A and B, compute appropriate relative position
//rigid_body constraints project the position/orientation of a rigid body
// -given current position of all joints on body, compute appropriate position and orientation



class node_base{ {
  //a quat and a vec3, thats about all we need to
  //describe a rigid body, maybe turn that into a
  //biquaternion later
  public:
    DEFINE_CREATE_FUNC(node_base)
    node_base() {}
    virtual quat R(){return _R;}
    virtual vec3 p(){return _p;}
    virtual ~node_base() {}
}

class port_base{
    //a joint is a connection between two rigid bodies
    //it has a relative rotation and translation
    public:
      DEFINE_CREATE_FUNC(port_base)
      port_base() {}
      virtual ~port_base() {}
      virtual quat R(){return _R;}
      virtual vec3 p(){return _p;}
      class joint_type _type;
  }




}
} // namespace gaudi::henson