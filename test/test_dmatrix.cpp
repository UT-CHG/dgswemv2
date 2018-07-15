#include "general_definitions.hpp"

int main(int argc, char** argv) {
    DynMatrix<double> mm;
    mm = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};

    DynMatrix<DynMatrix<double>> mmm;

    mmm.resize(2, 2);

    mmm(0, 0) = mm;
    mmm(1, 0) = mm;
    mmm(0, 1) = mm;
    mmm(1, 1) = mm;

    std::cout << mmm << std::endl;
}