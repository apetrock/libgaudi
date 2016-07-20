/*
 *  manifold_singleton.cpp
 *  Phase Vocoder
 *
 *  Created by John Delaney on 12/29/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "manifold_singleton.h"

namespace m2 {
  bool ID::instance_flag = false;
  ID* ID::global_instance = NULL;
	
  ID& ID::get_instance(){
    static ID* id;
    if (!id->initialized()) {
      id = new ID();
      id->initialized() = true;
    }
    return *id;
  }
	
  int ID::new_face_id(){
    ID& id = ID::get_instance();
    int out = id.get_next_face();
    return out;
  }
	
  int ID::new_edge_id(){
    ID& id = ID::get_instance();
    int out = id.get_next_edge();
    return out;
  }
	
  int ID::new_vertex_id(){
    ID& id = ID::get_instance();
    int out = id.get_next_vertex();
    return out;
  }
	
  int ID::new_face_vertex_id(){
    ID& id = ID::get_instance();
    int out = id.get_next_face_vertex();
    return out;
  }
	
}//end m2
