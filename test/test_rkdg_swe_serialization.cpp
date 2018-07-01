#include "general_definitions.hpp"
#include "problem/SWE/discretization_RKDG/data_structure/rkdg_swe_data.hpp"

#include "utilities/almost_equal.hpp"

using Utilities::almost_equal;

bool double_vectors_are_same(const std::vector<double>& a, const std::vector<double>& b) {
    bool are_same{ a.size() == b.size() };

    for ( uint i = 0; i < std::min(a.size(), b.size()); ++i ) {
        are_same = are_same && ( almost_equal(a[i],b[i]) );
    }
    return are_same;
}

SWE::RKDG::Boundary build_boundary(const std::size_t ngp, double& counter) {
    SWE::RKDG::Boundary bdry(ngp);
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

bool check_equality(const SWE::RKDG::Boundary& o_bdry, const SWE::RKDG::Boundary& i_bdry) {
return  double_vectors_are_same(o_bdry.ze_at_gp,i_bdry.ze_at_gp)
    && double_vectors_are_same(o_bdry.qx_at_gp,i_bdry.qx_at_gp)
    && double_vectors_are_same(o_bdry.qy_at_gp,i_bdry.qy_at_gp)
    && double_vectors_are_same(o_bdry.ze_numerical_flux_at_gp,i_bdry.ze_numerical_flux_at_gp)
    && double_vectors_are_same(o_bdry.qx_numerical_flux_at_gp,i_bdry.qx_numerical_flux_at_gp)
    && double_vectors_are_same(o_bdry.qy_numerical_flux_at_gp,i_bdry.qy_numerical_flux_at_gp)
    && double_vectors_are_same(o_bdry.bath_at_gp,i_bdry.bath_at_gp);
}

bool test_rkdg_swe_boundary() {
    const std::size_t ngp{4};
    double counter{0};
    SWE::RKDG::Boundary o_bdry = build_boundary(ngp,counter);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_bdry;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::RKDG::Boundary i_bdry;
    i_archive >> i_bdry;

    return check_equality(o_bdry, i_bdry);
}

SWE::RKDG::State build_state(const std::size_t ndof, double& counter) {
    SWE::RKDG::State state(ndof);
    for ( uint i = 0; i < ndof; ++i ) {
        state.ze[i] = counter++;
        state.qx[i] = counter++;
        state.qy[i] = counter++;
        state.bath[i] = counter++;

        state.rhs_ze[i] = counter++;
        state.rhs_qx[i] = counter++;
        state.rhs_qy[i] = counter++;

        state.solution_ze[i] = counter++;
        state.solution_qx[i] = counter++;
        state.solution_qy[i] = counter++;
    }
    return state;
}

bool check_equality(const SWE::RKDG::State& o_state, const SWE::RKDG::State& i_state) {
    return  double_vectors_are_same(o_state.ze,i_state.ze) &&
        double_vectors_are_same(o_state.qx,i_state.qx) &&
        double_vectors_are_same(o_state.qy,i_state.qy) &&
        double_vectors_are_same(o_state.bath,i_state.bath) &&
        double_vectors_are_same(o_state.rhs_ze,i_state.rhs_ze) &&
        double_vectors_are_same(o_state.rhs_qx,i_state.rhs_qx) &&
        double_vectors_are_same(o_state.rhs_qy,i_state.rhs_qy) &&
        double_vectors_are_same(o_state.solution_ze,i_state.solution_ze) &&
        double_vectors_are_same(o_state.solution_qx,i_state.solution_qx) &&
        double_vectors_are_same(o_state.solution_qy,i_state.solution_qy);
}

bool test_rkdg_swe_state() {
    const std::size_t ndof{5};
    double counter{0};
    SWE::RKDG::State o_state = build_state(ndof,counter);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_state;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::RKDG::State i_state;
    i_archive >> i_state;

    return check_equality(o_state, i_state);
}

SWE::RKDG::Internal build_internal(const std::size_t ngp, double& counter) {
    SWE::RKDG::Internal _int(ngp);

    for ( uint i = 0; i < ngp; ++i ) {
        _int.ze_flux_at_gp[0][i] = counter++;
        _int.ze_flux_at_gp[1][i] = counter++;
        _int.qx_flux_at_gp[0][i] = counter++;
        _int.qx_flux_at_gp[1][i] = counter++;
        _int.qy_flux_at_gp[0][i] = counter++;
        _int.qy_flux_at_gp[1][i] = counter++;

        _int.tau_s_at_gp[0][i] = counter++;
        _int.tau_s_at_gp[1][i] = counter++;
        _int.dp_atm_at_gp[0][i] = counter++;
        _int.dp_atm_at_gp[1][i] = counter++;
        _int.dtide_pot_at_gp[0][i] = counter++;
        _int.dtide_pot_at_gp[1][i] = counter++;

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

bool check_equality(const SWE::RKDG::Internal& o_int, const SWE::RKDG::Internal& i_int) {
return  double_vectors_are_same(o_int.ze_flux_at_gp[0],i_int.ze_flux_at_gp[0]) &&
    double_vectors_are_same(o_int.ze_flux_at_gp[1],i_int.ze_flux_at_gp[1]) &&
    double_vectors_are_same(o_int.qx_flux_at_gp[0],i_int.qx_flux_at_gp[0]) &&
    double_vectors_are_same(o_int.qx_flux_at_gp[1],i_int.qx_flux_at_gp[1]) &&
    double_vectors_are_same(o_int.qy_flux_at_gp[0],i_int.qy_flux_at_gp[0]) &&
    double_vectors_are_same(o_int.qy_flux_at_gp[1],i_int.qy_flux_at_gp[1]) &&
    double_vectors_are_same(o_int.tau_s_at_gp[0],i_int.tau_s_at_gp[0]) &&
    double_vectors_are_same(o_int.tau_s_at_gp[1],i_int.tau_s_at_gp[1]) &&
    double_vectors_are_same(o_int.dp_atm_at_gp[0],i_int.dp_atm_at_gp[0]) &&
    double_vectors_are_same(o_int.dp_atm_at_gp[1],i_int.dp_atm_at_gp[1]) &&
    double_vectors_are_same(o_int.dtide_pot_at_gp[0],i_int.dtide_pot_at_gp[0]) &&
    double_vectors_are_same(o_int.dtide_pot_at_gp[1],i_int.dtide_pot_at_gp[1]) &&
    double_vectors_are_same(o_int.ze_source_term_at_gp, i_int.ze_source_term_at_gp) &&
    double_vectors_are_same(o_int.qx_source_term_at_gp, i_int.qx_source_term_at_gp) &&
    double_vectors_are_same(o_int.qy_source_term_at_gp, i_int.qy_source_term_at_gp) &&
    double_vectors_are_same(o_int.bath_at_gp,i_int.bath_at_gp) &&
    double_vectors_are_same(o_int.h_at_gp, i_int.h_at_gp) &&
    double_vectors_are_same(o_int.bath_deriv_wrt_x_at_gp, i_int.bath_deriv_wrt_x_at_gp) &&
    double_vectors_are_same(o_int.bath_deriv_wrt_y_at_gp, i_int.bath_deriv_wrt_y_at_gp);
}

bool test_rkdg_swe_internal() {
    const std::size_t ngp{4};
    double counter{0};
    SWE::RKDG::Internal o_int = build_internal(ngp,counter);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_int;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::RKDG::Internal i_int;
    i_archive >> i_int;

    return check_equality(o_int,i_int);
}

SWE::RKDG::Spherical build_spherical(const uint nnode, const uint nbound, const uint ngp_internal,
                                     const std::vector<uint>& ngp_boundary, double& counter) {
    SWE::RKDG::Spherical sp(nnode,nbound,ngp_internal, ngp_boundary);

    for ( auto& x : sp.x_node) {
        x = counter++;
    }

    for ( auto& y : sp.y_node) {
        y = counter++;
    }

    for ( auto& gp_int : sp.sp_at_gp_internal ) {
        gp_int = counter++;
    }

    for ( auto& vec : sp.sp_at_gp_boundary ) {
        for ( auto& gp_bd : vec ) {
            gp_bd = counter ++;
        }
    }

    return sp;
}

bool check_equality(const SWE::RKDG::Spherical& o_sp, const SWE::RKDG::Spherical& i_sp) {
    bool are_same =   double_vectors_are_same(o_sp.x_node, i_sp.x_node) &&
        double_vectors_are_same(o_sp.y_node, i_sp.y_node) &&
        double_vectors_are_same(o_sp.sp_at_gp_internal, i_sp.sp_at_gp_internal);

    if ( o_sp.sp_at_gp_boundary.size() != i_sp.sp_at_gp_boundary.size() ) {
        return false;
    }

    for ( uint i = 0; i < o_sp.sp_at_gp_boundary.size(); ++i ) {
        are_same &= double_vectors_are_same(o_sp.sp_at_gp_boundary[i],
                                            i_sp.sp_at_gp_boundary[i]);
    }

    return are_same;
}

bool test_rkdg_swe_spherical() {
    uint nnode{3};
    uint nbound{3};
    uint ngp_internal{4};
    std::vector<uint> ngp_boundary{3,4,5};
    double counter{0};
    SWE::RKDG::Spherical o_sp = build_spherical(nnode,nbound,ngp_internal,ngp_boundary,counter);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_sp;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::RKDG::Spherical i_sp;
    i_archive >> i_sp;

    return check_equality(o_sp,i_sp);
}

SWE::RKDG::Source build_source(const uint nnode, double& counter) {
    SWE::RKDG::Source src(nnode);
    src.coriolis_f = 3.14;
    src.manning = true;
    src.g_manning_n_sq = -1.0;

    for ( uint n = 0; n < nnode; ++n ) {
        src.tau_s[0][n] = counter++;
        src.tau_s[1][n] = counter++;
        src.p_atm[n] = counter++;

        src.tide_pot[n] = counter++;
        src.manning_n[n] = counter++;
    }

    return src;
}

bool check_equality(const SWE::RKDG::Source& o_src, const SWE::RKDG::Source& i_src) {
    return almost_equal(o_src.coriolis_f, i_src.coriolis_f) &&
        (o_src.manning == i_src.manning) &&
        (o_src.coriolis_f == i_src.coriolis_f) &&
        almost_equal(o_src.coriolis_f, i_src.coriolis_f) &&
        double_vectors_are_same(o_src.tau_s[0],i_src.tau_s[0]) &&
        double_vectors_are_same(o_src.tau_s[1],i_src.tau_s[1]) &&
        double_vectors_are_same(o_src.p_atm,i_src.p_atm) &&
        double_vectors_are_same(o_src.tide_pot,i_src.tide_pot) &&
        double_vectors_are_same(o_src.manning_n,i_src.manning_n);
}

bool test_rkdg_swe_source() {
    uint nnode{4};
    double counter{0};
    SWE::RKDG::Source o_src = build_source(nnode,counter);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_src;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::RKDG::Source i_src;
    i_archive >> i_src;

    return check_equality(o_src,i_src);
}

SWE::RKDG::SlopeLimit build_slope_limit(const std::size_t nbound, double& counter) {
    SWE::RKDG::SlopeLimit sl(nbound);
    for ( uint i = 0; i < nbound; ++i ) {
        for ( uint dim = 0; dim < 2; ++dim ) {
            sl.surface_normal[i][dim] = counter++;
            sl.midpts_coord[i][dim] = counter++;
            sl.baryctr_coord_neigh[i][dim] = counter++;
        }

        sl.alpha_1[i] = counter++;
        sl.alpha_2[i] = counter++;
        sl.r_sq[i] = counter++;

        sl.ze_lin[i] = counter++;
        sl.qx_lin[i] = counter++;
        sl.qy_lin[i] = counter++;

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

        sl.wet_neigh[i] = (bool)( i % 2 );
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

bool check_equality( const SWE::RKDG::SlopeLimit& o_sl, const SWE::RKDG::SlopeLimit& i_sl ) {
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

    is_equal &= double_vectors_are_same(o_sl.ze_lin, i_sl.ze_lin);
    is_equal &= double_vectors_are_same(o_sl.qx_lin, i_sl.qx_lin);
    is_equal &= double_vectors_are_same(o_sl.qy_lin, i_sl.qy_lin);

    is_equal &= almost_equal(o_sl.ze_at_baryctr,i_sl.ze_at_baryctr);
    is_equal &= almost_equal(o_sl.qx_at_baryctr,i_sl.qx_at_baryctr);
    is_equal &= almost_equal(o_sl.qy_at_baryctr,i_sl.qy_at_baryctr);

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

    is_equal &= almost_equal(o_sl.baryctr_coord[0],i_sl.baryctr_coord[0]);
    is_equal &= almost_equal(o_sl.baryctr_coord[1],i_sl.baryctr_coord[1]);

    is_equal &= almost_equal(o_sl.ze_at_baryctr, i_sl.ze_at_baryctr);
    is_equal &= almost_equal(o_sl.qx_at_baryctr, i_sl.qx_at_baryctr);
    is_equal &= almost_equal(o_sl.qy_at_baryctr, i_sl.qy_at_baryctr);
    is_equal &= almost_equal(o_sl.bath_at_baryctr, i_sl.bath_at_baryctr);

    is_equal &= (o_sl.wet_neigh == i_sl.wet_neigh);

    return is_equal;
}

bool test_rkdg_swe_slope_limit() {
    const std::size_t nbound{10};
    double counter{1};

    SWE::RKDG::SlopeLimit o_sl = build_slope_limit(nbound,counter);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_sl;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::RKDG::SlopeLimit i_sl;
    i_archive >> i_sl;

    return check_equality(o_sl, i_sl);
}

SWE::RKDG::WetDry build_wet_dry(const std::size_t nvrtx, double& counter) {
    SWE::RKDG::WetDry wd(nvrtx);

    for ( uint i = 0; i < nvrtx; ++i ) {
        wd.ze_lin[i] = counter++;
        wd.qx_lin[i] = counter++;
        wd.qy_lin[i] = counter++;

        wd.ze_at_vrtx[i] = counter++;
        wd.qx_at_vrtx[i] = counter++;
        wd.qy_at_vrtx[i] = counter++;
        wd.bath_at_vrtx[i] = counter++;
        wd.h_at_vrtx[i] = counter++;
        wd.h_at_vrtx_temp[i] = counter++;
    }
    wd.wet = true;
    wd.went_completely_dry = false;
    wd.bath_min = counter++;

    return wd;
}

bool check_equality( const SWE::RKDG::WetDry& o_wd, const SWE::RKDG::WetDry& i_wd ) {
    return double_vectors_are_same(o_wd.ze_lin,i_wd.ze_lin) &&
        double_vectors_are_same(o_wd.qx_lin,i_wd.qx_lin) &&
        double_vectors_are_same(o_wd.qy_lin,i_wd.qy_lin) &&
        double_vectors_are_same(o_wd.ze_at_vrtx,i_wd.ze_at_vrtx) &&
        double_vectors_are_same(o_wd.qx_at_vrtx,i_wd.qx_at_vrtx) &&
        double_vectors_are_same(o_wd.qy_at_vrtx,i_wd.qy_at_vrtx) &&
        double_vectors_are_same(o_wd.bath_at_vrtx,i_wd.bath_at_vrtx) &&
        double_vectors_are_same(o_wd.h_at_vrtx,i_wd.h_at_vrtx) &&
        double_vectors_are_same(o_wd.h_at_vrtx_temp, i_wd.h_at_vrtx_temp) &&
        almost_equal(o_wd.bath_min, i_wd.bath_min) &&
        (o_wd.went_completely_dry == i_wd.went_completely_dry) &&
        (o_wd.wet == i_wd.wet);
}

bool test_rkdg_swe_wet_dry() {
    const std::size_t nvrtx{3};
    double counter{0};
    SWE::RKDG::WetDry o_wd = build_wet_dry(nvrtx,counter);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_wd;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::RKDG::WetDry i_wd;
    i_archive >> i_wd;


    return check_equality(o_wd,i_wd);
}

SWE::RKDG::Data build_data() {
    SWE::RKDG::Data data;

    uint nvrtx{3};
    uint nbound{3};
    uint ndof {14};
    uint nstates{3};
    uint ngp_internal{12};
    uint nnode{3};
    std::vector<uint> ngp_boundary{8,9,10};

    data.set_nvrtx(nvrtx);
    data.set_nbound(nbound);
    data.set_ndof(ndof);
    data.set_ngp_internal(ngp_internal);
    data.set_nnode(nnode);

    for ( uint i = 0; i < ngp_boundary.size(); ++i ) {
        data.set_ngp_boundary(i,ngp_boundary[i]);
    }

    double counter{0};
    data.state.push_back(build_state(ndof,counter));
    data.resize(nstates);

    for ( uint s = 0; s < nstates; ++s ) {
        data.state.at(s) = build_state(ndof,counter);
    }

    data.internal = build_internal( data.get_ngp_internal(), counter );

    for ( uint b = 0; b < nbound; ++b ) {
        data.boundary.push_back(build_boundary( data.get_ngp_boundary(b), counter ));
    }

    data.spherical_projection = build_spherical(nnode,nbound,ngp_internal,ngp_boundary,counter);
    data.source = build_source(nnode,counter);
    data.wet_dry_state = build_wet_dry(nvrtx,counter);
    data.slope_limit_state = build_slope_limit(nbound,counter);

    return data;
}

bool check_equality(const SWE::RKDG::Data& o_data, const SWE::RKDG::Data& i_data) {
    bool is_equal{true};
    for ( uint i =0; i < o_data.state.size(); ++i ) {
        is_equal &= check_equality( o_data.state.at(i), i_data.state.at(i) );
    }

    for ( uint j = 0; j < o_data.boundary.size(); ++j ) {
        is_equal &= check_equality( o_data.boundary.at(j), i_data.boundary.at(j) );
    }

    if ( o_data.get_nbound() != i_data.get_nbound() ) {
        return false;
    }

    for ( uint j = 0; j < o_data.get_nbound(); ++j ) {
        is_equal &= ( o_data.get_ngp_boundary(j) == i_data.get_ngp_boundary(j) );
    }

    return ( o_data.get_nnode() == i_data.get_nnode() )
        && ( o_data.get_nvrtx() == i_data.get_nvrtx() )
        && ( o_data.get_ndof() == i_data.get_ndof() )
        && ( o_data.get_ngp_internal() == i_data.get_ngp_internal() )
        && check_equality( o_data.spherical_projection, i_data.spherical_projection )
        && check_equality(o_data.source, i_data.source)
        && check_equality(o_data.wet_dry_state, i_data.wet_dry_state)
        && check_equality(o_data.slope_limit_state, i_data.slope_limit_state)
        && check_equality(o_data.internal, i_data.internal)
        && is_equal;
}

bool test_rkdg_swe_data() {
    SWE::RKDG::Data o_data = build_data();

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_data;

    hpx::serialization::input_archive i_archive(buffer);
    SWE::RKDG::Data i_data;
    i_archive >> i_data;

    return check_equality(o_data, i_data);
}

int main() {

    bool error_found{false};

    if ( !test_rkdg_swe_boundary() ) {
        std::cerr << "Error: Serialization of SWE::Boundary produces incorrect output\n";
        error_found  = true;
    }

    if ( !test_rkdg_swe_state() ) {
        std::cerr << "Error: Serialization of SWE::State produces incorrect output\n";
        error_found = true;
    }

    if ( !test_rkdg_swe_internal() ) {
        std::cerr << "Error: Serialization of SWE::Internal produces incorrect output\n";
        error_found = true;
    }

    if ( !test_rkdg_swe_spherical() ) {
        std::cerr << "Error: Serialization of SWE::RKDG::Spherical produces incorrect output\n";
        error_found = true;
    }

    if ( !test_rkdg_swe_source() ) {
        std::cerr << "Error: Serialization of SWE::RKDG::Source produces incorrect output\n";
        error_found = true;
    }

    if ( !test_rkdg_swe_slope_limit() ) {
        std::cerr << "Error:Serialization of SWE::SlopeLimit produces incorrect output\n";
        error_found = true;
    }

    if ( !test_rkdg_swe_wet_dry() ) {
        std::cerr << "Error: Serialization of SWE::WetDry produces incorrect output\n";
        error_found = true;
    }

    if ( !test_rkdg_swe_data() ) {
        std::cerr << "Error: Serialization of SWE::Data produces incorrect output\n";
        error_found = true;
    }

    if ( error_found ) {
        return 1;
    }

    return 0;
}