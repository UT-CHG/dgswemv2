#ifndef CLASS_MESH_H
#define CLASS_MESH_H
#include<memory>
#include<type_traits>
#include<vector>

#include "class_element.h"
#include "class_edge.h"

class MESH {
private: 
  std::vector<ELEMENT*> elements;
  std::vector<EDGE*> edges;
public:
  MESH();

  uint number_of_elements() { return elements.size(); }
  uint number_of_edges() { return edges.size(); }

  const ELEMENT* get_elt(int i) { return elements.at(i); }
  const EDGE*    get_edg(int i) { return edges.at(i); }



};

#endif
