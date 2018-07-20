#include "general_definitions.hpp"

int main(int argc, char* argv[]) {
    DynVector<double> vec{0, 1, 2, 3, 4, 5, 6, 7, 8};

    auto vector = transpose(vec);

    DynMatrix<double> mat(3, 3, vec.data());
    DynMatrix<double> mmat(3, 3, vector.data());

    std::cout << vector << std::endl;
    std::cout << vec << std::endl;
    std::cout << mat << std::endl;
    std::cout << mmat << std::endl;

    StatVector<double, 9> vect{0, 1, 2, 3, 4, 5, 6, 7, 8};
    StatMatrix<double, 3, 3> matt = reshape<double, 3, 3, 9>(vect);

    std::cout << vec << std::endl;

    return 0;
}