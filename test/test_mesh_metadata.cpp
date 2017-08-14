#include "preprocessor/ADCIRC_reader/adcirc_format.hpp"
#include "preprocessor/mesh_metadata.hpp"

const static auto is_equal = [](const MeshMetaData& meshA,
                                const MeshMetaData& meshB)-> bool {

    bool is_the_same {true};

    if ( meshA._mesh_name.compare(meshB._mesh_name) != 0 ) {
        std::cerr << "Strings not equal\n";
        is_the_same = false;
    }

    if ( meshA._elements != meshB._elements ) {
        std::cerr << "Elements not equal\n";
        is_the_same = false;
    }

    if ( meshA._nodes != meshB._nodes ) {
        std::cerr << "Nodes not equal\n";
        is_the_same = false;
    }

    return is_the_same;
};

int main(int argc, char** argv) {

    AdcircFormat mesh1(argv[1]);
    MeshMetaData meshA(mesh1);

    std::string out_name = argv[1];
    out_name = out_name + ".meta.out";

    meshA.WriteTo(out_name);

    MeshMetaData meshB(out_name);

    if (!is_equal(meshA, meshB)) {
        return 1;
    }

    return 0;
}
