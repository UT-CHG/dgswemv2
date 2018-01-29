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
    for ( uint i = 0; i < ngp; ++i ) {
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

bool test_swe_data_state() {
    const std::size_t ndof{5};
    SWE::State o_state(ndof);
    for ( uint i = 0; i < ndof; ++i ) {
        o_state.ze[i] = i;
        o_state.qx[i] = i + 1*ndof;
        o_state.qy[i] = i + 2*ndof;
        o_state.bath[i] = i + 3*ndof;

        o_state.rhs_ze[i] = i + 4*ndof;
        o_state.rhs_qx[i] = i + 5*ndof;
        o_state.rhs_qy[i] = i + 6*ndof;
    }

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_state;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::State i_state;
    i_archive >> i_state;

    bool error_found = !( double_vectors_are_same(o_state.ze,i_state.ze) &&
                          double_vectors_are_same(o_state.qx,i_state.qx) &&
                          double_vectors_are_same(o_state.qy,i_state.qy) &&
                          double_vectors_are_same(o_state.bath,i_state.bath) &&
                          double_vectors_are_same(o_state.rhs_ze,i_state.rhs_ze) &&
                          double_vectors_are_same(o_state.rhs_qx,i_state.rhs_qx) &&
                          double_vectors_are_same(o_state.rhs_qy,i_state.rhs_qy)
        );

    return error_found;
}

int main() {

    bool error_found{false};

    if ( test_swe_data_boundary() ) {
        std::cerr << "Error: Serialization of SWE::Boundary produces incorrect output\n";
        error_found  = true;
    }

    if ( test_swe_data_state() ) {
        std::cerr << "Error: Serialization of SWE::State produces incorrect output\n";
        error_found = true;
    }

    if ( error_found ) {
        return 1;
    }

    return 0;
}