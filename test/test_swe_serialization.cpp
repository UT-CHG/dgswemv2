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

SWE::Boundary build_boundary(const std::size_t ngp, double& counter) {
    SWE::Boundary bdry(ngp);
    for ( uint i = 0; i < ngp; ++i ) {
        bdry.ze_at_gp[i] = counter++;
        bdry.qx_at_gp[i] = counter++;
        bdry.qy_at_gp[i] = counter++;

        bdry.ze_numerical_flux_at_gp[i] = counter++;
        bdry.qx_numerical_flux_at_gp[i] = counter++;
        bdry.qy_numerical_flux_at_gp[i] = counter++;

        bdry.bath_at_gp[i] = counter++;
    }
    return bdry;
}

bool check_equality(const SWE::Boundary& o_bdry, const SWE::Boundary& i_bdry) {
return  double_vectors_are_same(o_bdry.ze_at_gp,i_bdry.ze_at_gp)
    && double_vectors_are_same(o_bdry.qx_at_gp,i_bdry.qx_at_gp)
    && double_vectors_are_same(o_bdry.qy_at_gp,i_bdry.qy_at_gp)
    && double_vectors_are_same(o_bdry.ze_numerical_flux_at_gp,i_bdry.ze_numerical_flux_at_gp)
    && double_vectors_are_same(o_bdry.qx_numerical_flux_at_gp,i_bdry.qx_numerical_flux_at_gp)
    && double_vectors_are_same(o_bdry.qy_numerical_flux_at_gp,i_bdry.qy_numerical_flux_at_gp)
    && double_vectors_are_same(o_bdry.bath_at_gp,i_bdry.bath_at_gp);
}

bool test_swe_data_boundary() {
    const std::size_t ngp{4};
    double counter{0};
    SWE::Boundary o_bdry = build_boundary(ngp,counter);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_bdry;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::Boundary i_bdry;
    i_archive >> i_bdry;

    return check_equality(o_bdry, i_bdry);
}

SWE::State build_state(const std::size_t ndof, double& counter) {
    SWE::State state(ndof);
    for ( uint i = 0; i < ndof; ++i ) {
        state.ze[i] = counter++;
        state.qx[i] = counter++;
        state.qy[i] = counter++;
        state.bath[i] = counter++;

        state.rhs_ze[i] = counter++;
        state.rhs_qx[i] = counter++;
        state.rhs_qy[i] = counter++;
    }
    return state;
}

bool check_equality(const SWE::State& o_state, const SWE::State& i_state) {
    return  double_vectors_are_same(o_state.ze,i_state.ze) &&
        double_vectors_are_same(o_state.qx,i_state.qx) &&
        double_vectors_are_same(o_state.qy,i_state.qy) &&
        double_vectors_are_same(o_state.bath,i_state.bath) &&
        double_vectors_are_same(o_state.rhs_ze,i_state.rhs_ze) &&
        double_vectors_are_same(o_state.rhs_qx,i_state.rhs_qx) &&
        double_vectors_are_same(o_state.rhs_qy,i_state.rhs_qy);
}

bool test_swe_data_state() {
    const std::size_t ndof{5};
    double counter{0};
    SWE::State o_state = build_state(ndof,counter);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_state;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::State i_state;
    i_archive >> i_state;

    return check_equality(o_state, i_state);
}

SWE::Internal build_internal(const std::size_t ngp, double& counter) {
    SWE::Internal _int(ngp);

    for ( uint i = 0; i < ngp; ++i ) {
        _int.ze_flux_at_gp[0][i] = counter++;
        _int.ze_flux_at_gp[1][i] = counter++;
        _int.qx_flux_at_gp[0][i] = counter++;
        _int.qx_flux_at_gp[1][i] = counter++;
        _int.qy_flux_at_gp[0][i] = counter++;
        _int.qy_flux_at_gp[1][i] = counter++;

        _int.ze_source_term_at_gp[i] = counter++;
        _int.qx_source_term_at_gp[i] = counter++;
        _int.qy_source_term_at_gp[i] = counter++;

        _int.ze_at_gp[i] = counter++;
        _int.qx_at_gp[i] = counter++;
        _int.qy_at_gp[i] = counter++;
        _int.bath_at_gp[i] = counter++;
        _int.h_at_gp[i] = counter++;

        _int.bath_deriv_wrt_x_at_gp[i] = counter++;
        _int.bath_deriv_wrt_y_at_gp[i] = counter++;
    }
    return _int;
}

bool check_equality(const SWE::Internal& o_int, const SWE::Internal& i_int) {
return  double_vectors_are_same(o_int.ze_flux_at_gp[0],i_int.ze_flux_at_gp[0]) &&
    double_vectors_are_same(o_int.ze_flux_at_gp[1],i_int.ze_flux_at_gp[1]) &&
    double_vectors_are_same(o_int.qx_flux_at_gp[0],i_int.qx_flux_at_gp[0]) &&
    double_vectors_are_same(o_int.qx_flux_at_gp[1],i_int.qx_flux_at_gp[1]) &&
    double_vectors_are_same(o_int.qy_flux_at_gp[0],i_int.qy_flux_at_gp[0]) &&
    double_vectors_are_same(o_int.qy_flux_at_gp[1],i_int.qy_flux_at_gp[1]) &&
    double_vectors_are_same(o_int.ze_source_term_at_gp, i_int.ze_source_term_at_gp) &&
    double_vectors_are_same(o_int.qx_source_term_at_gp, i_int.qx_source_term_at_gp) &&
    double_vectors_are_same(o_int.qy_source_term_at_gp, i_int.qy_source_term_at_gp) &&
    double_vectors_are_same(o_int.bath_at_gp,i_int.bath_at_gp) &&
    double_vectors_are_same(o_int.h_at_gp, i_int.h_at_gp) &&
    double_vectors_are_same(o_int.bath_deriv_wrt_x_at_gp, i_int.bath_deriv_wrt_x_at_gp) &&
    double_vectors_are_same(o_int.bath_deriv_wrt_y_at_gp, i_int.bath_deriv_wrt_y_at_gp);
}

bool test_swe_data_internal() {
    const std::size_t ngp{4};
    double counter{0};
    SWE::Internal o_int = build_internal(ngp,counter);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_int;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::Internal i_int;
    i_archive >> i_int;

    return check_equality(o_int,i_int);
}

SWE::SlopeLimit build_slope_limit(const std::size_t nbound, double& counter) {
    SWE::SlopeLimit sl(nbound);
    for ( uint i = 0; i < nbound; ++i ) {
        for ( uint dim = 0; dim < 2; ++dim ) {
            sl.surface_normal[i][dim] = counter++;
            sl.midpts_coord[i][dim] = counter++;
            sl.baryctr_coord_neigh[i][dim] = counter++;
        }

        sl.alpha_1[i] = counter++;
        sl.alpha_2[i] = counter++;
        sl.r_sq[i] = counter++;

        sl.ze_at_vrtx[i] = counter++;
        sl.qx_at_vrtx[i] = counter++;
        sl.qy_at_vrtx[i] = counter++;

        sl.ze_at_midpts[i] = counter++;
        sl.qx_at_midpts[i] = counter++;
        sl.qy_at_midpts[i] = counter++;
        sl.bath_at_midpts[i] = counter++;

        sl.ze_at_baryctr_neigh[i] = counter++;
        sl.qx_at_baryctr_neigh[i] = counter++;
        sl.qy_at_baryctr_neigh[i] = counter++;
        sl.bath_at_baryctr_neigh[i] = counter++;
    }


    for ( uint j = 0; j < 3; ++j ) {
        sl.w_midpt_char[j] = counter++;
        sl.delta_char[j] = counter++;

        for ( uint c = 0; c < 3; ++c ) {
            sl.w_baryctr_char[c][j] = counter++;
            sl.delta[c][j] = counter++;
            sl.L[c][j] = counter++;
            sl.R[c][j] = counter++;
        }
    }

    sl.baryctr_coord[0] = counter++;
    sl.baryctr_coord[1] = counter++;

    sl.ze_at_baryctr = counter++;
    sl.qx_at_baryctr = counter++;
    sl.qy_at_baryctr = counter++;
    sl.bath_at_baryctr = counter++;

    return sl;
}

bool check_equality( const SWE::SlopeLimit& o_sl, const SWE::SlopeLimit& i_sl ) {
    bool is_equal{true};

    for ( uint i = 0; i < o_sl.midpts_coord.size(); ++i ) {
        is_equal &= double_vectors_are_same(o_sl.surface_normal[i], i_sl.surface_normal[i]);

        is_equal &= almost_equal(o_sl.midpts_coord[i][0], i_sl.midpts_coord[i][0]);
        is_equal &= almost_equal(o_sl.midpts_coord[i][1], i_sl.midpts_coord[i][1]);

        is_equal &= almost_equal(o_sl.baryctr_coord_neigh[i][0], i_sl.baryctr_coord_neigh[i][0]);
        is_equal &= almost_equal(o_sl.baryctr_coord_neigh[i][1], i_sl.baryctr_coord_neigh[i][1]);
    }

    is_equal &= double_vectors_are_same(o_sl.w_midpt_char, i_sl.w_midpt_char);
    is_equal &= double_vectors_are_same(o_sl.delta_char, i_sl.delta_char);

    for ( uint c = 0; c < 3; ++c ) {
        is_equal &= double_vectors_are_same(o_sl.w_baryctr_char[c], i_sl.w_baryctr_char[c]);
        is_equal &= double_vectors_are_same(o_sl.delta[c], i_sl.delta[c]);
        is_equal &= double_vectors_are_same(o_sl.L[c], i_sl.L[c]);
        is_equal &= double_vectors_are_same(o_sl.R[c], i_sl.R[c]);
    }

    is_equal &= double_vectors_are_same(o_sl.alpha_1, i_sl.alpha_1);
    is_equal &= double_vectors_are_same(o_sl.alpha_2, i_sl.alpha_2);
    is_equal &= double_vectors_are_same(o_sl.r_sq, i_sl.r_sq);


    is_equal &= double_vectors_are_same(o_sl.ze_at_vrtx, i_sl.ze_at_vrtx);
    is_equal &= double_vectors_are_same(o_sl.qx_at_vrtx, i_sl.qx_at_vrtx);
    is_equal &= double_vectors_are_same(o_sl.qy_at_vrtx, i_sl.qy_at_vrtx);

    is_equal &= double_vectors_are_same(o_sl.ze_at_midpts, i_sl.ze_at_midpts);
    is_equal &= double_vectors_are_same(o_sl.qx_at_midpts, i_sl.qx_at_midpts);
    is_equal &= double_vectors_are_same(o_sl.qy_at_midpts, i_sl.qy_at_midpts);
    is_equal &= double_vectors_are_same(o_sl.bath_at_midpts, i_sl.bath_at_midpts);

    is_equal &= double_vectors_are_same(o_sl.ze_at_baryctr_neigh, i_sl.ze_at_baryctr_neigh);
    is_equal &= double_vectors_are_same(o_sl.qx_at_baryctr_neigh, i_sl.qx_at_baryctr_neigh);
    is_equal &= double_vectors_are_same(o_sl.qy_at_baryctr_neigh, i_sl.qy_at_baryctr_neigh);
    is_equal &= double_vectors_are_same(o_sl.bath_at_baryctr_neigh, i_sl.bath_at_baryctr_neigh);

    is_equal &= almost_equal(o_sl.baryctr_coord[0],i_sl.baryctr_coord[0]);
    is_equal &= almost_equal(o_sl.baryctr_coord[1],i_sl.baryctr_coord[1]);

    is_equal &= almost_equal(o_sl.ze_at_baryctr, i_sl.ze_at_baryctr);
    is_equal &= almost_equal(o_sl.qx_at_baryctr, i_sl.qx_at_baryctr);
    is_equal &= almost_equal(o_sl.qy_at_baryctr, i_sl.qy_at_baryctr);
    is_equal &= almost_equal(o_sl.bath_at_baryctr, i_sl.bath_at_baryctr);

    return is_equal;
}

bool test_swe_data_slope_limit() {
    const std::size_t nbound{10};
    double counter{1};

    SWE::SlopeLimit o_sl = build_slope_limit(nbound,counter);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_sl;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::SlopeLimit i_sl;
    i_archive >> i_sl;

    return check_equality(o_sl, i_sl);
}

SWE::WetDry build_wet_dry(const std::size_t nvrtx, double& counter) {
    SWE::WetDry wd(nvrtx);

    for ( uint i = 0; i < nvrtx; ++i ) {
        wd.ze_at_vrtx[i] = counter++;
        wd.qx_at_vrtx[i] = counter++;
        wd.qy_at_vrtx[i] = counter++;
        wd.bath_at_vrtx[i] = counter++;
        wd.h_at_vrtx[i] = counter++;
        wd.h_at_vrtx_temp[i] = counter++;
    }
    wd.wet = true;
    wd.bath_min = counter++;
    wd.water_volume = counter++;

    return wd;
}

bool check_equality( const SWE::WetDry& o_wd, const SWE::WetDry& i_wd ) {
    return double_vectors_are_same(o_wd.ze_at_vrtx,i_wd.ze_at_vrtx) &&
        double_vectors_are_same(o_wd.qx_at_vrtx,i_wd.qx_at_vrtx) &&
        double_vectors_are_same(o_wd.qy_at_vrtx,i_wd.qy_at_vrtx) &&
        double_vectors_are_same(o_wd.bath_at_vrtx,i_wd.bath_at_vrtx) &&
        double_vectors_are_same(o_wd.h_at_vrtx,i_wd.h_at_vrtx) &&
        double_vectors_are_same(o_wd.h_at_vrtx_temp, i_wd.h_at_vrtx_temp) &&
        almost_equal(o_wd.bath_min, i_wd.bath_min) &&
        almost_equal(o_wd.water_volume, i_wd.water_volume) &&
        (o_wd.wet == i_wd.wet);
}

bool test_swe_data_wet_dry() {
    const std::size_t nvrtx{3};
    double counter{0};
    SWE::WetDry o_wd = build_wet_dry(nvrtx,counter);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_wd;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::WetDry i_wd;
    i_archive >> i_wd;


    return check_equality(o_wd,i_wd);
}

SWE::Data build_data() {
    SWE::Data data;

    uint nvrtx{3};
    uint nbound{3};
    uint ndof {14};
    uint nstates{3};
    uint ngp_internal{12};

    data.set_nvrtx(nvrtx);
    data.set_nbound(nbound);
    data.set_ndof(ndof);
    data.set_ngp_internal(ngp_internal);



    for ( uint i = 0; i < 3; ++i ) {
        data.set_ngp_boundary(i,8+i);
    }

    double counter{0};
    data.wet_dry_state = build_wet_dry(nvrtx,counter);
    data.slope_limit_state = build_slope_limit(nbound,counter);
    data.state.push_back(build_state(ndof,counter));
    data.resize(nstates);

    for ( uint s = 0; s < nstates; ++s ) {
        data.state.at(s) = build_state(ndof,counter);
    }


    data.internal = build_internal( data.get_ngp_internal(), counter );

    for ( uint b = 0; b < nbound; ++b ) {
        data.boundary.push_back(build_boundary( data.get_ngp_boundary(b), counter ));
    }

    return data;
}

bool check_equality(const SWE::Data& o_data, const SWE::Data& i_data) {
    bool is_equal{true};
    for ( uint i =0; i < o_data.state.size(); ++i ) {
        is_equal &= check_equality( o_data.state.at(i), i_data.state.at(i) );
    }

    for ( uint j = 0; j < o_data.boundary.size(); ++j ) {
        is_equal &= check_equality( o_data.boundary.at(j), i_data.boundary.at(j) );
    }

    return check_equality(o_data.wet_dry_state, i_data.wet_dry_state)
        && check_equality(o_data.slope_limit_state, i_data.slope_limit_state)
        && check_equality(o_data.internal, i_data.internal)
        && is_equal;
}

bool test_swe_data() {
    SWE::Data o_data = build_data();

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_data;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::Data i_data;
    i_archive >> i_data;

    return check_equality(o_data, i_data);
}

int main() {

    bool error_found{false};

    if ( !test_swe_data_boundary() ) {
        std::cerr << "Error: Serialization of SWE::Boundary produces incorrect output\n";
        error_found  = true;
    }

    if ( !test_swe_data_state() ) {
        std::cerr << "Error: Serialization of SWE::State produces incorrect output\n";
        error_found = true;
    }

    if ( !test_swe_data_internal() ) {
        std::cerr << "Error: Serialization of SWE::Internal produces incorrect output\n";
        error_found = true;
    }

    if ( !test_swe_data_slope_limit() ) {
        std::cerr << "Error:Serialization of SWE::SlopeLimit produces incorrect output\n";
        error_found = true;
    }

    if ( !test_swe_data_wet_dry() ) {
        std::cerr << "Error: Serialization of SWE::WetDry produces incorrect output\n";
        error_found = true;
    }

    if ( !test_swe_data() ) {
        std::cerr << "Error: Serialization of SWE::Data produces incorrect output\n";
        error_found = true;
    }

    if ( error_found ) {
        return 1;
    }

    return 0;
}