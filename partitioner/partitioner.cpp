#include "preprocessor/input_parameters.hpp"
#include "preprocessor/ADCIRC_reader/adcirc_format.hpp"
#include "preprocessor/mesh_metadata.hpp"

#include "numa_configuration.hpp"

#include <cassert>
#include <chrono>
#include <vector>

std::vector<std::vector<MeshMetaData>> partition(const MeshMetaData& mesh_meta,
                                                 const int num_partitions,
                                                 const int num_nodes,
                                                 const NumaConfiguration& numa_config);

void write_distributed_edge_metadata(const std::string& file_name,
                                     const InputParameters& input,
                                     const MeshMetaData& mesh_meta,
                                     const std::vector<std::vector<MeshMetaData>>& submeshes);

int main(int argc, char** argv) {

    std::cout << "?????????????????????????????????????????????????????????????"
                 "????????\n";
    std::cout << "?  Mesh Preprocessor\n";
    std::cout << "?????????????????????????????????????????????????????????????"
                 "????????\n\n";

    if (argc < 4 || argc > 5) {
        std::cout << "Usage:\n";
        std::cout << "  path/to/partitioner <input_file_name> <number of "
                     "partitions>\n";
        std::cout << "                      <number of nodes> <NUMA "
                     "configuration>(optional)\n";

        return 0;
    }


    std::cout << "Mesh Partitioner Configuration\n";
    InputParameters input(argv[1]);
    std::cout << "  Input File: " << argv[1] << '\n';
    std::string path_to_input(argv[1]);
    path_to_input = path_to_input.substr(0,path_to_input.find_last_of("/\\")+1);
    std::string input_mesh_str(path_to_input + input.mesh_file_name);
    std::cout << "  Mesh Name: " << input_mesh_str << '\n';
    int num_partitions = atoi(argv[2]);
    std::cout << "  Number of partitions: " << num_partitions << '\n';
    int num_nodes = atoi(argv[3]);
    std::cout << "  Number of compute nodes: " << num_nodes << '\n';

    NumaConfiguration numa_config;
    if (argc == 5) {
        std::string numa_str(argv[4]);
        numa_config = NumaConfiguration(numa_str);
        std::cout << "  NUMA configuration: " << numa_str << "n\n";
    } else {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
        numa_config = NumaConfiguration("default");
        std::cout << "  NUMA configuration: default (1 NUMA domain per node)\n\n";
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    AdcircFormat mesh(input_mesh_str);
    MeshMetaData mesh_meta(mesh);
    std::vector<std::vector<MeshMetaData>> submeshes = partition(mesh_meta, num_partitions, num_nodes, numa_config);
    for (uint n = 0; n < submeshes.size(); ++n) {
        for (uint m = 0; m < submeshes[n].size(); ++m) {
            std::string outname = input_mesh_str;
            outname = outname.substr(0,outname.find_last_of(".")-1);
            
            outname += "_" + std::to_string(static_cast<long long>(n)) + "_" +
                       std::to_string(static_cast<long long>(m)) + ".meta";
	    
            submeshes[n][m].WriteTo(outname);
        }
    }

    write_distributed_edge_metadata(input_mesh_str, input, mesh_meta, submeshes);

    // finish out by writing updated output file
    std::string updated_input_filename = std::string(argv[1]);
    updated_input_filename.erase(updated_input_filename.size() - 3);
    updated_input_filename += "_parallelized.15";
    input.mesh_format = "Meta";
    input.WriteTo(updated_input_filename);

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "\nTime Elapsed (in us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
              << std::endl;
}
