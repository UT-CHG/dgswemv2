#ifndef SHAPE_FACTORY_H
#define SHAPE_FACTORY_H

#include "util.hxx"

enum class element_type {
  element_tri
};

enum class edge_type {

};

template<class... Args>
static ELEMENT* create_elt(element_type elt_typ, Args&&... args )
{
  switch (elt_typ) {
  case element_type::element_tri:
    return create<ELEMENT_TRI>(std::forward<Args>(args)...);
  default:
    return nullptr;
  }
}

template<class... Args>
static EDGE* create_edg(edge_type edg_typ, Args&&... args )
{
  switch (edg_typ) {
  default:
    return nullptr;
  }
}

#endif