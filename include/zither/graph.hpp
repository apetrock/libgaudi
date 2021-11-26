#include <iostream>
#include <memory>
#include <vector>
#include <map>
#include <bitset>
#include <string>

#ifndef __NET_GRAPH__
#define __NET_GRAPH__

namespace zither {
  template <class T, int S>
  class Buffer{
    int mRead;
    int mWrite;
    
    std::array<T,S> mArray;
  };


  class Process {
  public:
    Process(){};
    virtual std::vector<std::string> & getPorts() {return mPorts;}
    virtual void process(){localProcess();}
    virtual void localProcess() = 0;
    std::vector<std::string> mPorts;
    
  };


  using ProcPtr = std::shared_ptr<zither::Process>;


class SinProc : public Process {
public:
  SinProc() : 
    Process(),
    freq(440.0),
    phase(0.0)
  {
    mPorts.push_back(std::string("frequency"));
    mPorts.push_back(std::string("phase"));
    mPorts.push_back(std::string("output"));
  }
  virtual void localProcess(){std::cout << "sin proc" << std::endl;};
  double freq, phase;
};

class OutProc : public Process {
public:
  OutProc() :
  Process()
  {
    mPorts.push_back("left_input"); 
    mPorts.push_back("right_input");
  }
  virtual void localProcess(){std::cout << "out proc" << std::endl;};
  
};

  class edge;
  class node;
  class network_graph;

  using GraphType = zither::network_graph;  
  using EdgeType = zither::edge;
  using NodeType = zither::node;
  
  using GraphPtr = std::shared_ptr<GraphType>;  
  using EdgePtr =  std::shared_ptr<EdgeType>;
  using NodePtr = std::shared_ptr<NodeType>;

  class node{
  public:
  
    node(){
      
      /*
        mProcess = std::make_shared<Process>();
        mPorts = mProcess->getPorts();
        int i = 0;
        for(auto port : mPorts){
        mPortIds[port] = i;
        i++;
        }
      */
    }

    node(std::shared_ptr<Process> proc) : mData(proc){
      
      /*
        mProcess = std::make_shared<Process>();
        mPorts = mProcess->getPorts();
        int i = 0;
        for(auto port : mPorts){
        mPortIds[port] = i;
        i++;
        }
      */
    }
    
    void addEdge(EdgePtr e){
      mEdges.push_back(e);
    }
    
    void setId(int id){mId = id;};
    int id() const {return mId;};
    
    int portId(std::string port) {return mPortIds[port];};
    
    bool visited() const {return bits[0];} 
    void setVisited(bool v) {bits[0] = v;} 
    
    std::shared_ptr<Process> data() {return mData;} 
    
    void setData(std::shared_ptr<Process> data) {
      mData = data;
    } 

    std::vector<EdgePtr> & getEdges(){return mEdges;}
    
  protected:
    int mId;
    std::vector<EdgePtr>       mEdges;
    std::map<std::string, int> mPortIds;
    std::vector<std::string>   mPorts;
    std::bitset<32> bits;

    std::shared_ptr<Process> mData;
  };

  
  template <typename PROC>
  class ProcNode : public node {
  public:
    ProcNode(){
      
      std::shared_ptr<PROC> proc = std::make_shared<PROC>();
      std::vector<std::string> & ports = proc->getPorts();
      for(auto p : ports){
        mPortIds[p] = mPorts.size();
        mPorts.push_back(p);  
      }
      //node(std::dynamic_pointer_cast<Process>(proc));
      setData(std::dynamic_pointer_cast<Process>(proc));
    }

    std::shared_ptr<PROC> getProc(){
      return std::dynamic_pointer_cast<std::shared_ptr<PROC>>(data());
    }

    std::shared_ptr<Process> getBaseProc(){
      return data();
    }
    
  protected:
    std::shared_ptr<PROC>  mProcess;
  };
 
  template <typename PROC>
  std::shared_ptr<node> 
  getBaseNode(std::shared_ptr<ProcNode<PROC>>pn){
    return std::dynamic_pointer_cast<node>(pn);
  }



  class edge{
  public:
    
    edge(){};
    edge(size_t n0, size_t p0,
         size_t n1, size_t p1):
      node0(n0),
      node1(n1),
      port0(p0),
      port1(p1)
    {}

    bool outgoing(int id){
      return id == node0;
    }
    
    bool incoming(int id){
      return id == node1;
    }
    
    bool other(int id){
      return (id == node0) ? node1 : node0;
    }

    size_t node0, node1;
    size_t port0, port1;
  };



  class network_graph{
  public:

    network_graph(){}
    
    NodePtr addNode(NodePtr n){
      
      int id = mNodes.size();
      n->setId(id);
      mNodes.push_back(n);
      return n;
    }

    template <typename PROC>
    std::shared_ptr<ProcNode<PROC>> addNode(){
      std::shared_ptr<ProcNode<PROC>> pn = 
        std::make_shared<ProcNode<PROC>>();
      std::shared_ptr<node> n = getBaseNode(pn);
      this->addNode(n);
      return pn;
    }

    EdgePtr connectNode(size_t n0, size_t p0, 
                        size_t n1, size_t p1){
      
      EdgePtr e = std::make_shared<EdgeType>(n0, p0,
                                             n1, p1);
      NodePtr node0 = mNodes[n0];
      NodePtr node1 = mNodes[n1];
      node0->addEdge(e);
      node1->addEdge(e);
      return e;
    }
    
    EdgePtr connectNode(NodePtr n0, std::string p0, 
                        NodePtr n1, std::string p1){
      size_t nId0 = n0->id();
      size_t nId1 = n1->id();
      size_t pId0 = n0->portId(p0);
      size_t pId1 = n1->portId(p1);
      return this->connectNode(nId0, pId0,
                       nId1, pId1);
      
    }
    
    std::vector<NodePtr> & nodes(){return mNodes;}
    std::vector<EdgePtr> & edges(){return mEdges;}
  private:
    std::vector<NodePtr> mNodes;
    std::vector<EdgePtr> mEdges;
  };
}
#endif
