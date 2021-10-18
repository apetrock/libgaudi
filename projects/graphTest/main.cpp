#include <iostream>
#include <string>
#include <memory>
#include <stack>

#include "zither/graph.hpp"

namespace zither {


  
  /*
    L ← Empty list that will contain the sorted nodes
    while there are unmarked nodes do
    select an unmarked node n
    visit(n) 
    function visit(node n)
    if n has a temporary mark then stop (not a DAG)
    if n is not marked (i.e. has not been visited yet) then
    mark n temporarily
    for each node m with an edge from n to m do
    visit(m)
    mark n permanently
    unmark n temporarily
    add n to head of L
  */

  std::vector<NodePtr> topoSort(GraphPtr g){
    std::vector<NodePtr> sorted;
    std::vector<NodePtr> & nodes = g->nodes();
    std::stack<NodePtr> stack;
    
    for(auto n : nodes){
      n->data()->process();
      stack.push(n);
    }
    while(stack.size() > 0){
      NodePtr cn = stack.top();
      stack.pop();
      if(cn->visited()) continue;
      for(auto e : cn->getEdges()){
        if(e->incoming(cn->id()))
          stack.push(nodes[e->other(cn->id())]);
      }
      cn->setVisited(true);
      sorted.push_back(cn);
    }
    return sorted;
  }
}

zither::GraphPtr graph;

int main(int argc, char *argv[]) {
  graph = std::make_shared<zither::GraphType>();
  std::shared_ptr<zither::ProcNode<zither::OutProc>> n0 = 
    graph->addNode<zither::OutProc>();
  std::shared_ptr<zither::ProcNode<zither::SinProc>> n1 = 
    graph->addNode<zither::SinProc>();
  std::shared_ptr<zither::ProcNode<zither::SinProc>> n2 = 
    graph->addNode<zither::SinProc>();

  graph->connectNode(getBaseNode(n1), "output", getBaseNode(n0), "left_input");
  graph->connectNode(getBaseNode(n2), "output", getBaseNode(n0), "right_input");
  std::vector<zither::NodePtr> sortedNodes = topoSort(graph);

  
  for(auto n : sortedNodes){
    std::shared_ptr<zither::Process> p = 
      std::dynamic_pointer_cast<zither::Process>(n->data());
    p->process();
  }

return 0;
}
