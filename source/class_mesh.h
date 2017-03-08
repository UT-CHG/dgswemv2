#ifndef CLASS_MESH_H
#define CLASS_MESH_H
#include<vector>

#include "class_element.h"
#include "class_edge.h"
#include "mesh/shape_factory.h"
#include "mesh/singleton_factory.h"

struct elt_sig
{
  int ID;
  element_enumerator elt_typ;
  basis_enumerator basis_typ;
  basis_geom_enumerator basis_geom;
  double* x_nodal_coordinates; 
  double* y_nodal_coordinates;
  //... what ever else
};

struct edg_sig
{
  int ID;
  edge_enumerator edg_typ;
};

struct MESH_DATA {
  std::vector<elt_sig> elts;
  std::vector<edg_sig> edgs;
};


enum class MESH_FORMAT {
  adcirc_mesh
};

class MESH
{
private:
  SINGLETON_FACTORY<basis> bases;
  SINGLETON_FACTORY<basis_geometry> geometric_bases;
  SINGLETON_FACTORY<integration_rule> integration_rules;

  std::vector<std::unique_ptr<ELEMENT> > elements;
  std::vector<std::unique_ptr<EDGE> > edges;
public:
  MESH(std::string mesh_file_id, std::string params_file_id);

  uint number_of_elements() const { return elements.size(); }
  uint number_of_edges() const { return edges.size(); }

  const ELEMENT* get_elt(int i) const { return elements.at(i); }
  const EDGE*    get_edg(int i) const { return edges.at(i); }
};

#endif
