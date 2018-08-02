#include <hpx/runtime/serialization/serialize.hpp>
#include "utilities/linear_algebra/use_blaze.hpp"
#include "utilities/almost_equal.hpp"

#include <iostream>

template <typename Vector>
bool test_serialize_vector(Vector& i_vec) {
    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << i_vec;

    hpx::serialization::input_archive i_archive(buffer);
    Vector o_vec;
    i_archive >> o_vec;

    bool error_found = (o_vec.size() != i_vec.size()) || (o_vec.capacity() != i_vec.capacity());

    error_found |= !Utilities::almost_equal(0., blaze::norm(i_vec - o_vec));

    return error_found;
}

int main() {
    bool error_found{false};

    {  // serialize vector
        StatVector<double, 3> a{1, 2, 3};
        error_found |= test_serialize_vector(a);

        DynVector<double> b{4, 5, 6};
        error_found |= test_serialize_vector(b);

        SparseVector<double> c{1, 0, 3, 0, 5};
        error_found |= test_serialize_vector(c);
    }
    return error_found;
}