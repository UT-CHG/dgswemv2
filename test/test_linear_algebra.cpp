#include "general_definitions.hpp"

int main() {
#ifdef USE_BLAZE  
    DynMatrix<double> A(3,3); 
    DynMatrix<double> B(3,3);

    A = {{1,2,3},{4,5,6},{7,8,9}};
    B = {{1,2,3},{4,5,6},{7,8,9}};

    DynMatrix<double> C = cwise_division(A,B);

    std::cout << A << B << C << "works\n";
#endif

#ifdef USE_EIGEN  
    DynVector<double> a(9);
    a << 1,2,3,4,5,6,7,8,9;

    DynMatrix<double> A = reshape<double, 3>(a);

    std::cout << a << '\n' << A << "\nworks\n";
#endif
}
