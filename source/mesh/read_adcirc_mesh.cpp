#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

/*! \file */

/*! /fn void read_adcirc_mesh(const std::string& fort14)
 *   /brief Constructs mesh from adcirc mesh file
 *   /param fort14 Name of the mesh file
 *
 *   The mesh is constructed by determining the edge connectivity
 *   by enumerating the edges using a bit shift to define a single
 *   long integer for which the first 32 bits correspond to the
 *   lower of the two nodes attached to that edge, and the other
 *   node corresponds to the remaining 32 bits. We then use a hash
 *   table to determine which edges are neighbors.
 */
void read_adcirc_mesh(const std::string& fort14) {
  std::ifstream ifs(fort14);

  if (!ifs) {
    std::cerr << "Fatal Error: Mesh named" << fort14 << " not found \n";
    exit(1);
  }

  int  number_nodes, number_eles;

  //read in mesh name
  ifs.ignore(1000,'\n');

  ifs >> number_eles;
  ifs >> number_nodes;

  std::cout << "It has "<< number_eles << " elements,\n";
  std::cout << "and " << number_nodes << " nodes.\n";

  //read node locations
  for (int j = 0; j < number_nodes+1; ++j)
    ifs.ignore(1000,'\n');


  //read in mesh connectivity information
  //store edges in a hash table
  std::unordered_map<uint64_t,std::pair<int,int> > edgedict;
  std::vector<uint64_t> node(3);
  int temp;

  for (int j = 0; j < number_eles; ++j)
    {

      //process a line
      //file is formatted as element
      //[element number, number of nodes, node names]
      ifs >> temp;
      ifs >> temp;

      ifs >> node[ 0 ];
      ifs >> node[ 1 ];
      ifs >> node[ 2 ];

      ifs.ignore(1000,'\n');

      for(int k = 0; k <3; ++k)
	{
	  uint64_t curr_key =  (std::min(node[k],node[(k+1)%3]))<<32
	    | std::max(node[k],node[(k+1)%3]);

	  if (edgedict.count(curr_key)) //if already one element has the edge.
	    {
	      auto& edge_neigh = edgedict.at(curr_key);
	      edge_neigh.second = j;

	      edgedict.at(curr_key) = edge_neigh;
	    }
	  else
	    {
	      std::pair<int,int> edge_neigh {j,-1};
	      edgedict.insert({curr_key, edge_neigh});
	    }
	}
    }
}
