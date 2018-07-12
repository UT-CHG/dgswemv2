#include "general_definitions.hpp"

int main(int argc, char** argv) {
    DMatrix<double> mm;

    // mm.resize(3, 3);

    // blaze::submatrix(mm, 0, 0, 3, 3) = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};

    mm = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};

    std::cout << mm << std::endl;
}