#include "preprocessor/ADCIRC_reader/adcirc_format.hpp"
#include "preprocessor/mesh_metadata.hpp"

const static auto is_equal = [](const MeshMetaData& meshA, const MeshMetaData& meshB) -> bool {
    bool is_the_same{true};

    if (meshA.mesh_name.compare(meshB.mesh_name) != 0) {
        std::cerr << "Strings not equal\n";
        is_the_same = false;
    }

    if (meshA.elements != meshB.elements) {
        std::cerr << "Elements not equal\n";
        is_the_same = false;
    }

    if (meshA.nodes != meshB.nodes) {
        std::cerr << "Nodes not equal\n";
        is_the_same = false;
    }

    return is_the_same;
};

int main(int argc, char** argv) {
    bool error_found{false};
    for (int i = 1; i < argc; ++i) {
        AdcircFormat mesh1(argv[1]);
        MeshMetaData meshA(mesh1);

        std::string out_name = argv[1];
        out_name             = out_name + ".meta.out";

        meshA.write_to(out_name);

        MeshMetaData meshB(out_name);

        if (!is_equal(meshA, meshB)) {
            error_found = true;
            std::cerr << "Error: in reading a writing mesh: " << out_name << '\n' << "       for MeshMeta format.\n";
        }
    }

    return error_found;
}
