#ifndef MESH_METADATA_HPP
#define MESH_METADATA_HPP

#include "../general_definitions.hpp"
#include "../problem/SWE/swe_definitions.hpp"
#include "adcirc_format.hpp"

// ID: corresponds to the element ID
// coordinates: correspond to the vertices of the element
//   moving counter clockwise around the triangle
// boundary type: coresponds to the type of the element edge (Interface, Tidal boundary, etc.)
// neighbor_IDs: corresponds to the neighbor ID if it exists. If it does not exist the default ID is set to
//   DEFAULT_ID as set in general_definitions.hpp
struct ElementMetaData
{
  ElementMetaData(uint n_faces)
    : coordinates(n_faces), boundary_type(n_faces), neighbor_IDs(n_faces)
  {}

  std::vector<Point<2>> coordinates;
  std::vector<SWE::BoundaryConditions> boundary_type;
  std::vector<uint> neighbor_IDs;
};


class MeshMetaData
{
public:
  MeshMetaData(const AdcircFormat& mesh_file);

  using const_iterator = std::vector<ElementMetaData>::const_iterator;

  inline
  const_iterator cbegin()
  { return _meta.cbegin(); }

  inline
  const_iterator cend()
  { return _meta.cend(); }

private:
  std::unordered_map<uint,ElementMetaData> _meta;
};


#endif
