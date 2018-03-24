#include "preprocessor/ADCIRC_reader/adcirc_format.hpp"

#include <iostream>

const static auto compare = [](AdcircFormat meshA, AdcircFormat meshB) -> bool {
    if ((meshA.name == meshB.name) && (meshA.nodes == meshB.nodes) && (meshA.elements == meshB.elements) &&
        (meshA.NOPE == meshB.NOPE) && (meshA.NETA == meshB.NETA) && (meshA.NBDV == meshB.NBDV) &&
        (meshA.IBTYPE == meshB.IBTYPE) && (meshA.NBOU == meshB.NBOU) && (meshA.NVEL == meshB.NVEL) &&
        (meshA.IBTYPE == meshB.IBTYPE) && (meshA.NBVV == meshB.NBVV)) {
        return true;
    }

    return false;
};

int main(int argc, char** argv) {

    AdcircFormat mesh1(argv[1]);

    std::string out_name = argv[1];
    out_name             = out_name + ".out";

    mesh1.write_to(out_name.c_str());

    AdcircFormat mesh2(out_name.c_str());

    if (!compare(mesh1, mesh2)) {
        return 1;
    }

    return 0;
}
