//----------------------------------------------------------------
// HeapItem.h
// Simple class with which to build the heap demonstration.
//
// Author: Dr. Rick Coleman
//----------------------------------------------------------------
#ifndef HEAPITEM_H
#define HEAPITEM_H
#include <iostream>

using namespace std;
template  <typename ITEM_T>
class HeapItem {
private:
  ITEM_T  m_dData;                            // Dummy data value
   
public:
  int index;
  bool inHeap;
  //-----------------------------------
  // Default constructor
  //-----------------------------------
  HeapItem()
  {
    //m_dData = 0.0;
  }

  //-----------------------------------
  // Constructor
  //-----------------------------------
  HeapItem(ITEM_T data)
  {
    m_dData = data;
  }

  //-----------------------------------
  // Destructor
  //-----------------------------------
  ~HeapItem()
  {
  }

  //-----------------------------------
  // Return data item
  //-----------------------------------
  
  ITEM_T getData()
  {
    return m_dData;
  }
  
};

template <typename ITEM_T, typename KEY_T, typename COMPARE, typename ASSIGN>
class Heap
{
public:
  typedef HeapItem<ITEM_T> heap_t;
  std::vector<heap_t> m_Elements;
  COMPARE    * m_comp;
  ASSIGN     * m_assign;
  int          m_iNumElements;              // Number of elements in the heap
  int          m_iHeapLength;               // Size of the array
  

  Heap(int size, COMPARE* comp, ASSIGN *assign)
  {
    // Create heap of given size
    m_comp = comp;
    m_assign = assign;
    m_Elements.reserve(size);
  }

  int size(){
    return m_Elements.size();
  };

  int heap(ITEM_T object) {
    auto * newItem = new HeapItem<ITEM_T>(object);
    newItem->inHeap = true;
    m_Elements.push_back(*newItem);
    int loc = m_Elements.size()-1;

    (*m_assign)(newItem->getData(),loc);
    this->up(newItem->index = loc);
      
    return m_Elements.size();
  };

  
  heap_t getItem(int i){
    return m_Elements[i];
  }

  heap_t pop() {
    
    heap_t  removed = m_Elements[0];
    heap_t & object = m_Elements.back();

    removed.inHeap = false;

    m_Elements.pop_back();
    if (m_Elements.size()) {
      (*m_assign)(object.getData(),0);
      m_Elements[object.index = 0] = object;
      this->down(0);
    }
    return removed;    
  };


  int remove(int i) {
    auto removed = m_Elements[i];    
    removed.inHeap = false;

    auto object = m_Elements.back();
    m_Elements.pop_back();
    if (i != m_Elements.size()) {
      
      (*m_assign)(object.getData(),i);
      m_Elements[object.index = i] = object;
      
      if((*m_comp)(object.getData(), removed.getData()))
	this->down(i);
      else
	this->up(i);
    }
    return i;
  };

  int remove(heap_t removed) {
    auto i = removed.index;
    this->remove(i);
  };

  int up(int i) {
    auto object = m_Elements[i];
    while (i > 0) {
      int up = ((i + 1) >> 1) - 1;
      
      auto parent = m_Elements[up];
      if ((*m_comp)(object.getData(), parent.getData())) break;
            
      (*m_assign)(parent.getData(),i);
      (*m_assign)(object.getData(),up);

      m_Elements[parent.index = i] = parent;
      m_Elements[object.index = i = up] = object;
      
    }
  }

  int down(int i) {
    auto object = m_Elements[i];
    while (true) {
      int right = (i + 1) << 1;
      int left = right - 1;
      int down = i;
      auto child = m_Elements[down];
      if (left < m_Elements.size() && 
	  !(*m_comp)(m_Elements[left].getData(), child.getData())) 
	child = m_Elements[down = left];
      
      if (right < m_Elements.size() && 
	  !(*m_comp)(m_Elements[right].getData(), child.getData())) 
	child = m_Elements[down = right];

      if (down == i) break;
      
      (*m_assign)(child.getData(),i);
      (*m_assign)(object.getData(),down);
      
      m_Elements[child.index = i] = child;
      m_Elements[object.index = i = down] = object;
    }
  }

};


#endif


