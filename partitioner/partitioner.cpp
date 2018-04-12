#include "preprocessor/input_parameters.hpp"
#include "preprocessor/ADCIRC_reader/adcirc_format.hpp"
#include "preprocessor/mesh_metadata.hpp"

#include "problem/swe_partitioner_inputs.hpp"
#include "problem/default_partitioner_inputs.hpp"

#include <cassert>
#include <chrono>
#include <sstream>
#include <vector>

std::vector<std::vector<MeshMetaData>> partition(const MeshMetaData&                                 mesh_meta,
                                                 const std::unordered_map<int, std::vector<double>>& problem_weights,
                                                 const int                                           num_partitions,
                                                 const int                                           num_nodes,
                                                 const int                                           ranks_per_locality,
                                                 const bool                                          rank_balanced);

void write_distributed_edge_metadata(const std::string&                            file_name,
                                     const InputParameters<>&                      input,
                                     const MeshMetaData&                           mesh_meta,
                                     const std::vector<std::vector<MeshMetaData>>& submeshes);

int main(int argc, char** argv) {
    std::cout << "?????????????????????????????????????????????????????????????"
                 "????????\n";
    std::cout << "?  Mesh Preprocessor\n";
    std::cout << "?????????????????????????????????????????????????????????????"
                 "????????\n\n";

    if (argc < 4 || argc > 6) {
        std::cout << "Usage:\n";
        std::cout << "  path/to/partitioner <input_file_name> <number of "
                     "partitions>\n";
        std::cout << "                      <number of nodes> <ranks per locality>(optional) <rank "
                     "balanced>(optional)\n";

        return 0;
    }

    std::cout << "Mesh Partitioner Configuration\n";
    InputParameters<> input(argv[1]);
    std::cout << "  Input File: " << argv[1] << '\n';
    std::string input_mesh_str(input.mesh_input.mesh_file_name);
    std::cout << "  Mesh Path: " << input_mesh_str << '\n';
    int num_partitions = atoi(argv[2]);
    std::cout << "  Number of partitions: " << num_partitions << '\n';
    int num_nodes = atoi(argv[3]);
    std::cout << "  Number of compute nodes: " << num_nodes << '\n';

    int ranks_per_locality{1};
    if (argc >= 5) {
        ranks_per_locality = atoi(argv[4]);
        std::cout << "  Number of ranks per locality: " << ranks_per_locality << '\n';
    };

    bool rank_balanced{false};
    if (argc == 6) {
        std::stringstream ss(argv[5]);
        ss >> std::boolalpha >> rank_balanced;
    }
    std::cout << std::boolalpha << "  Rank balanced: " << rank_balanced << "\n\n";

    auto t1 = std::chrono::high_resolution_clock::now();

    input.read_mesh();
    MeshMetaData& mesh_meta = input.mesh_input.mesh_data;

    // initialize problem specific inputs
    std::unique_ptr<ProblemPartitionerInputs> problem_inputs;
    if (input.problem_input.node["name"]) {
        if (input.problem_input.node["name"].as<std::string>() == "swe") {
            SWE::Inputs swe_inputs(input.problem_input.node);
            problem_inputs = std::make_unique<SWE::PartitionerInputs>(mesh_meta, swe_inputs);
        } else {
            problem_inputs = std::make_unique<DefaultPartitionerInputs>(mesh_meta);
        }
    } else {
        problem_inputs = std::make_unique<DefaultPartitionerInputs>(mesh_meta);
    }

    std::vector<std::vector<MeshMetaData>> submeshes = partition(
        mesh_meta, problem_inputs->GetWeights(), num_partitions, num_nodes, ranks_per_locality, rank_balanced);

    for (uint n = 0; n < submeshes.size(); ++n) {
        for (uint m = 0; m < submeshes[n].size(); ++m) {
            std::string outname = input_mesh_str;
            outname             = outname.substr(0, outname.find_last_of("."));

            outname += "_" + std::to_string(static_cast<long long>(n)) + "_" +
                       std::to_string(static_cast<long long>(m)) + ".meta";

            submeshes[n][m].write_to(outname);
        }
    }

    write_distributed_edge_metadata(input_mesh_str, input, mesh_meta, submeshes);

    problem_inputs->PartitionAuxiliaryFiles();

    // finish out by writing updated output file
    std::string updated_input_filename = std::string(argv[1]);
    updated_input_filename.erase(updated_input_filename.size() - 3);
    updated_input_filename += "_parallelized.15";
    input.mesh_input.mesh_file_name =
        input.mesh_input.mesh_file_name.substr(0, input.mesh_input.mesh_file_name.find_last_of(".")) + ".meta";
    input.mesh_input.mesh_format = "Meta";
    input.write_to(updated_input_filename);

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Time Elapsed (in us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
              << std::endl;
}
