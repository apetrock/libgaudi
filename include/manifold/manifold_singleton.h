/*
 *  manifold_singleton.h
 *  Phase Vocoder
 *
 *  Created by John Delaney on 12/29/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>

#ifndef __MANIFOLD_SINGLETON__
#define __MANIFOLD_SINGLETON__
namespace m2 {
  class ID{
  public:
		
    static int		 new_face_id();	
    static int		 new_vertex_id();	
    static int		 new_edge_id();
    static int		 new_face_vertex_id();
    static ID& get_instance();
    bool& initialized(){return instance_flag;}
    bool  initialized() const {return instance_flag;}
		
  private:
    ID(){
      face_id = 0;
      face_vertex_id = 0;
      edge_id = 0;
      vertex_id = 0;
      global_instance = this;
    }
		
    ID(const ID&);
    ID & operator=(const ID &);
		
    int	get_next_face(){
      int out_ = face_id;
      face_id ++;
      return out_;
    }
		
    int	get_next_edge(){
      int out_ = edge_id;
      edge_id ++;
      return out_;
    }
		
    int	get_next_vertex(){
      int out_ = vertex_id;
      vertex_id ++;
      return out_;
    }
		
    int	get_next_face_vertex(){
      int out_ = face_vertex_id;
      face_vertex_id ++;
      return out_;
    }
		
    static ID* global_instance;
    static bool instance_flag;
    int face_id;
    int face_vertex_id;
    int vertex_id;
    int edge_id;
  };
}//end singleton

#endif
