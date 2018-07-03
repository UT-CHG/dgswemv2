#include "preprocessor/ADCIRC_reader/adcirc_format.hpp"

#include <iostream>

const static auto compare = [](AdcircFormat meshA, AdcircFormat meshB) -> bool {
    if ((meshA.name == meshB.name) && (meshA.nodes == meshB.nodes) && (meshA.elements == meshB.elements) &&
        (meshA.NOPE == meshB.NOPE) && (meshA.NETA == meshB.NETA) && (meshA.NBDV == meshB.NBDV) &&
        (meshA.NBOU == meshB.NBOU) && (meshA.NVEL == meshB.NVEL) && (meshA.IBTYPE == meshB.IBTYPE) &&
        (meshA.NBVV == meshB.NBVV) && (meshA.BARINTH == meshB.BARINTH) && (meshA.BARINCFSB == meshB.BARINCFSB) &&
        (meshA.BARINCFSP == meshB.BARINCFSP) && (meshA.NGEN == meshB.NGEN) && (meshA.NGEN == meshB.NGEN) &&
        (meshA.NNGN == meshB.NNGN) && (meshA.NBGN == meshB.NBGN)) {
        return true;
    }

    return false;
};

int main(int argc, char** argv) {
    bool error_found{false};

    for (int i = 1; i < argc; ++i) {
        AdcircFormat mesh1(argv[i]);

        std::string out_name = argv[i];
        out_name             = out_name + ".out";

        mesh1.write_to(out_name.c_str());

        AdcircFormat mesh2(out_name.c_str());

        if (!compare(mesh1, mesh2)) {
            std::cerr << "Error: in writing and reading mesh " << out_name << '\n';
            error_found = true;
        }
    }

    return error_found;
}
