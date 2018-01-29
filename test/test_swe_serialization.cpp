#include "general_definitions.hpp"
#include "problem/SWE/data_structure/swe_data.hpp"

#include "utilities/almost_equal.hpp"

using Utilities::almost_equal;

bool double_vectors_are_same(const std::vector<double>& a, const std::vector<double>& b) {
    bool are_same{ a.size() == b.size() };

    for ( uint i = 0; i < std::min(a.size(), b.size()); ++i ) {
        are_same = are_same && ( almost_equal(a[i],b[i]) );
    }
    return are_same;
}

bool test_swe_data_boundary() {
    const std::size_t ngp{4};
    SWE::Boundary o_bdry(ngp);
    for ( uint i = 0; i < ngp; ++ i ) {
        o_bdry.ze_at_gp[i] = i;
        o_bdry.qx_at_gp[i] = i + 1*ngp;
        o_bdry.qy_at_gp[i] = i + 2*ngp;

        o_bdry.ze_numerical_flux_at_gp[i] = i + 3*ngp;
        o_bdry.qx_numerical_flux_at_gp[i] = i + 4*ngp;
        o_bdry.qy_numerical_flux_at_gp[i] = i + 5*ngp;

        o_bdry.bath_at_gp[i] = i + 6*ngp;
    }

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_bdry;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::Boundary i_bdry;
    i_archive >> i_bdry;

    bool error_found = !( double_vectors_are_same(o_bdry.ze_at_gp,i_bdry.ze_at_gp)
                          && double_vectors_are_same(o_bdry.qx_at_gp,i_bdry.qx_at_gp)
                          && double_vectors_are_same(o_bdry.qy_at_gp,i_bdry.qy_at_gp)
                          && double_vectors_are_same(o_bdry.ze_numerical_flux_at_gp,i_bdry.ze_numerical_flux_at_gp)
                          && double_vectors_are_same(o_bdry.qx_numerical_flux_at_gp,i_bdry.qx_numerical_flux_at_gp)
                          && double_vectors_are_same(o_bdry.qy_numerical_flux_at_gp,i_bdry.qy_numerical_flux_at_gp)
                          && double_vectors_are_same(o_bdry.bath_at_gp,i_bdry.bath_at_gp)
        );

    return error_found;
}

int main() {

    bool error_found{false};

    error_found |= test_swe_data_boundary();

    if ( error_found ) {
        return 1;
    }

    return 0;
}