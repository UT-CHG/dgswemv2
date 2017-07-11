#include "mesh_metadata.hpp"

MeshMetaData::MeshMetaData(const AdcircFormat& mesh_file)
{
 for ( const auto& elt : mesh_file.elements ) {

   if ( elt.first < 0 ) { throw std::logic_error( "ERROR: Element ID is negative; not supported\n" ); }

   uint ID = elt.first;
   _meta.insert({ID, ElementMetaData(3) });

   for ( uint i = 0; i < 3; ++i ) {
     std::array<double,3> tmp = mesh_file.nodes.at(elt.second[i+1]);
     _meta.at(ID).nodal_coordinates[i][0] = tmp[0];
     _meta.at(ID).nodal_coordinates[i][1] = tmp[1];
   }
  }

  {//make all edges
    using eltID_faceID = std::pair<int,int>;

    std::unordered_map<std::uint64_t,std::pair<eltID_faceID,eltID_faceID> > edge_dictionary;
    for ( const auto& elt : mesh_file.elements ) {
      std::vector<int> node { elt.second[1], elt.second[2], elt.second[3] };

      for ( int k = 0; k < 3; ++k ) {
	std::uint64_t curr_key =  (static_cast<uint>(std::min(node[k],node[(k+1)%3])))<<32
	  | static_cast<uint>(std::max(node[k],node[(k+1)%3]));

	if (edge_dictionary.count(curr_key)) {//if already one element has the edge.

	  edge_dictionary.at(curr_key).second = std::make_pair(elt.first,k);

	} else {

	  std::pair<eltID_faceID,eltID_faceID> edge_info {{elt.first,k},{-1,0}};
	  edge_dictionary.insert({curr_key, edge_info});

	}
      }
    }

    for ( const auto& edge : edge_dictionary ) {
      //check if there are two elements associated with this edge
      if ( edge.second.second.first != -1 ) {
	uint eltA_id = edge.second.first.first;
	uint eltB_id = edge.second.second.first;
	uint fidA = edge.second.first.second;
	uint fidB = edge.second.second.second;

	_meta.at(eltA_id).neighbor_ID[fidA] = eltB_id;
	_meta.at(eltA_id).boundary_type[fidA] = SWE::BoundaryConditions::internal;

	_meta.at(eltB_id).neighbor_ID[fidB] = eltA_id;
	_meta.at(eltB_id).boundary_type[fidB] = SWE::BoundaryConditions::internal;

      } else {
	//treat boundary conditions
	int elt_id = edge.second.first.first;
	int fid = edge.second.first.second;
/*
	std::array<int,2> nodes { mesh_file.elements[elt_id][fid],
	    mesh_file.elements[elt_id][ (fid+1) %3] };
/*
	_meta.at(elt_id).boundary_type[fid] = mesh_file.get_ibtype( nodes );
  */    }
    }
  }
}