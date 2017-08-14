#include "preprocessor/ADCIRC_reader/adcirc_format.hpp"
#include "preprocessor/mesh_metadata.hpp"

#include "numa_configuration.hpp"

#include <cassert>
#include <vector>

std::vector<std::vector<MeshMetaData> > partition(const MeshMetaData& mesh_meta,
                                                  const int num_partitions,
                                                  const int num_nodes,
                                                  const NumaConfiguration& numa_config);

int main(int argc, char **argv) {

    std::cout << "?????????????????????????????????????????????????????????????????????\n";
    std::cout << "?  Mesh Preprocessor\n";
    std::cout << "?????????????????????????????????????????????????????????????????????\n\n";

    if ( argc < 3 ) {
        std::cout << "\nUsage:\n";
        std::cout << "  path/to/partitioner <input_file_name> <number of partitions>\n";
        std::cout << "                      <number of nodes> <NUMA configuration>(optional)\n";

        return 0;
    }

    std::cout << "Mesh Partitioner Configuration\n";
    std::string input_mesh_str(argv[1]);
    std::cout << "  Mesh Name: " << input_mesh_str << '\n';
    int num_partitions = atoi(argv[2]);
    std::cout << "  Number of partitions: " << num_partitions << '\n';
    int num_nodes = atoi(argv[3]);
    std::cout << "  Number of compute nodes: " << num_nodes << '\n';

    NumaConfiguration numa_config;
    if (argc > 4 ) {
        std::string numa_str(argv[4]);
        numa_config = NumaConfiguration(numa_str);
        std::cout << "  NUMA configuration: " << numa_str << "n\n";
    } else {
        numa_config = NumaConfiguration("default");
        std::cout << "  NUMA configuration: default\n\n";
    }

    AdcircFormat mesh(input_mesh_str);
    mesh.write_to("test_mesh");
//    try {
        MeshMetaData mesh_meta(mesh);
        std::vector<std::vector<MeshMetaData>> submeshes = partition(mesh_meta, num_partitions,
                                                                     num_nodes, numa_config);
/*    for ( uint n = 0; n < submeshes.size(); ++n ) {
        for ( uint m = 0; n < submeshes[n].size(); ++m ) {
            std::string outname = input_mesh_str + "_" + std::to_string(static_cast<long long>(n))
                + "_" + std::to_string(static_cast<long long>(m));

            AdcircFormat tmp = submeshes[n].write_AdcircFormat();
            tmp.write_to(outname);
        }
        }*/

        //  } catch (const std::exception& e ) {
        //std::cerr << "Exception thrown in constructing MeshMetaData: " << e.what() << '\n';
        //return 1;
        //}
}
