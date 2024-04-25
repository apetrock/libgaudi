#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cxxabi.h>
#include <execinfo.h>
#include <iostream>
#include <memory.h>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>
#include <zlib.h>

#ifndef __ASAWA_GRAPH__
#define __ASAWA_GRAPH__

namespace gaudi {
namespace asawa {

class datum;

namespace graph {

  class graph{
    //a joint network is a set of joints
    //that connect rigid bodies
    public:
      DEFINE_CREATE_FUNC(graph)
      graph() {}
      virtual ~graph() {}
      virtual void add_body(quat R, vec3 p, node_base::ptr node){
        _nodes.push_back(node);
        _R.push_back(R);
        _p.push_back(p);
        _port_start.push_back(-1);
      }

      index_t other(index_t i) {
        return (i % 2 == 0) ? i + 1 : i - 1;
      }

      virtual void add_node(){
        _port_start.push_back(-1);
      }

      virtual void insert_port(const index_t & nid, const index_t & pid){
        _port_start[nid] = pid;
        
        _node[pid] = node0;
        indext_t pAc =  _port_start[nid];
        indext_t pAn =  _port_next[pAc];
        indext_t pAp =  _port_prev[pAc];
        _port_start[nid] = pid;
        _port_next[pAc] = pid;
        _port_prev[pAn] = pid;
        _port_next[pid] = pAn;
        _port_prev[pid] = pAc;

      }
            
      virtual void connect_nodes(const index_t & node0, const index_t & node1){

        index_t pid0 = _port_next.size();
        index_t pid1 = _port_next.size()+1;
        _port_next.resize(Nc+2);
        _port_prev.resize(Nc+2);
        _node.resize(Nc+2);

        insert_port(node0, pid0);
        insert_port(node1, pid1);
      }

      datum_ptr &get_datum(index_t i) { return __data[i]; }
      const datum_ptr &const_get_datum(index_t i) const { return __data[i]; }

      std::vector<datum_ptr> &get_data() { return __data; }

      std::vector<index_t> _port_start; //index of first port for each node
      std::vector<index_t> _port_next; //next port on the node for given port
      std::vector<index_t> _port_prev; //previous port on the node for given port
      std::vector<index_t> _node; //2*port count, index of node A/B for each port
      std::vector<datum_ptr> __data; //data for nodes/edges
  }
}
}}