#ifndef TEST_ELEMENT_TRIANGLE_HPP
#define TEST_ELEMENT_TRIANGLE_HPP
#include "general_definitions.hpp"
#include "utilities/almost_equal.hpp"
#include "geometry/mesh_definitions.hpp"
#include "preprocessor/input_parameters.hpp"
#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"
#include "problem/SWE/problem_function_files/swe_true_solution_functions.hpp"
#include "problem/SWE/discretization_RKDG/rkdg_swe_problem.hpp"

using MasterType  = Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>;
using ShapeType   = Shape::StraightTriangle;
using ElementType = Geometry::Element<2, MasterType, ShapeType, SWE::RKDG::Accessor>;

using Utilities::almost_equal;

const std::vector<double> IntegrationPhi_true = {
    6.881941874331419e-01, -8.89156081756484e-02, 7.216878364870322e-02, 1.20281306081172e-02, 0.00000000000000e+00,
    4.811252243246882e-03, 0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00};

const std::vector<double> IntegrationDPhiDX_true = {
    0.00000000000000e+00,  0.00000000000000e+00,  1.376388374866284,     0.00000000000000e+00,  6.212068893253612e-01,
    4.330127018922193e-01, 0.00000000000000e+00,  4.541823898251592e-01, 1.732050807568877e-01, 8.09844822527753e-01,
    0.00000000000000e+00,  3.593266739736606e-01, 8.66025403784439e-02,  4.752989329117831e-01, 2.886751345948129e-01,
    0.00000000000000e+00,  2.976844657418943e-01, 4.948716593053934e-02, 3.477531724792759e-01, 1.649572197684645e-01,
    5.677806551742285e-01, 0.00000000000000e+00,  2.542695064773368e-01, 3.092947870658709e-02, 2.799779589772887e-01,
    1.030982623552903e-01, 3.839882003961648e-01, 2.165063509461096e-01, 0.00000000000000e+00,  2.219874007710344e-01,
    2.061965247105806e-02, 2.368352966075515e-01, 6.873217490352688e-02, 2.924290225727635e-01, 1.443375672974064e-01,
    4.382557445972097e-01, 0.00000000000000e+00,  1.970208683922866e-01, 1.443375672974064e-02, 2.063521709566122e-01,
    4.811252243246881e-02, 2.388944320048427e-01, 1.010362971081845e-01, 3.216406511354539e-01, 1.732050807568877e-01,
    0.00000000000000e+00,  1.771259393312325e-01, 1.049727762162956e-02, 1.833727084180284e-01, 3.499092540543186e-02,
    2.038020106931869e-01, 7.348094335140692e-02, 2.541597625891523e-01, 1.259673314595547e-01, 3.571900656194553e-01,
    0.00000000000000e+00,  1.6089440460229e-01,   7.87295821622217e-03,  1.652850866122932e-01, 2.624319405407389e-02,
    1.788406944880994e-01, 5.511070751355518e-02, 2.112510229573667e-01, 9.4475498594666e-02,   2.765124421822678e-01,
    1.443375672974064e-01};

const std::vector<double> IntegrationDPhiDY_true = {
    0.00000000000000e+00,  2.383974596215561,      0.00000000000000e+00,  -1.744016935856293,    4.166666666666668e-01,
    5.639957660359269e-01, 2.337820631441114,      -3.00000000000000e-01, 3.443590455313351e-01, 2.00000000000000e-01,
    -1.879743494847109,    3.50000000000000e-01,   1.913828550551446e-01, 1.00000000000000e-01,  3.851221159381804e-01,
    2.252991532071854,     -3.238095238095238e-01, 2.131929837166789e-01, 5.714285714285716e-02, 2.453297042212707e-01,
    1.428571428571429e-01, -1.937912020128888,     3.392857142857143e-01, 1.178190458792872e-01, 3.571428571428573e-02,
    1.83043385083546e-01,  8.92857142857143e-02,   2.855387752616669e-01, 2.210576982387224,     -3.293650793650794e-01,
    1.641331427947897e-01, 2.380952380952382e-02,  1.483567763383696e-01, 5.952380952380954e-02, 2.021367719068524e-01,
    1.111111111111111e-01, -1.970227867507653,     3.361111111111112e-01, 8.21387896549881e-02,  1.666666666666667e-02,
    1.26044828084989e-01,  4.166666666666668e-02,  1.569447932368951e-01, 7.777777777777779e-02, 2.272196374916589e-01,
    2.185128252576446,     -3.313131313131313e-01, 1.367871463596913e-01, 1.212121212121213e-02, 1.102706931333926e-01,
    3.030303030303032e-02, 1.293162666349232e-01,  5.656565656565658e-02, 1.71701644642061e-01,  9.09090909090909e-02,
    -1.990792497657777,    3.348484848484849e-01,  6.041767889566378e-02, 9.0909090909091e-03,   9.83877593333601e-02,
    2.272727272727273e-02, 1.1080228133985e-01,    4.242424242424243e-02, 1.378321850860735e-01, 6.818181818181821e-02,
    1.8880147989604e-01};

bool check_for_error(ElementType& triangle, DynMatrix<double>& f_vals) {
    // Check integrations
    bool error_found{false};

    if (!almost_equal(IntegrationPhi_true[0], triangle.Integration(f_vals)[0], 1.e+04)) {
        error_found = true;
        std::cerr << "Error found in Triangle element in Integration" << std::endl;
    }

    for (uint dof = 0; dof < 66; ++dof) {
        if (!almost_equal(IntegrationPhi_true[dof], triangle.IntegrationPhi(dof, f_vals)[0], 1.e+04)) {
            error_found = true;
            std::cerr << "Error found in Triangle element in IntegrationPhi" << std::endl;
        }
    }

    for (uint dof = 0; dof < 66; ++dof) {
        if (!almost_equal(
                IntegrationDPhiDX_true[dof], triangle.IntegrationDPhi(GlobalCoord::x, dof, f_vals)[0], 1.e+04)) {
            error_found = true;
            std::cerr << "Error found in Triangle element in IntegrationDPhi "
                         "in x direction"
                      << std::endl;
        }
    }
    // Add 7 more modes
    for (uint dof = 0; dof < 66; ++dof) {
        if (!almost_equal(
                IntegrationDPhiDY_true[dof], triangle.IntegrationDPhi(GlobalCoord::y, dof, f_vals)[0], 1.e+04)) {
            error_found = true;

            std::cerr << "Error found in Triangle element in IntegrationDPhi "
                         "in y direction"
                      << std::endl;
        }
    }

    // Check linears through integration
    // u(x,y) = 3 + 2*x - 2*sqrt(3)*y plane
    DynMatrix<double> u(1, 3);
    u(0, 0) = 2.0;
    u(0, 1) = 4.0;
    u(0, 2) = 0.0;

    DynMatrix<double> u_gp(1, triangle.data.get_ngp_internal());
    DynMatrix<double> du_dx_gp(1, triangle.data.get_ngp_internal());
    DynMatrix<double> du_dy_gp(1, triangle.data.get_ngp_internal());

    u_gp     = triangle.ComputeLinearUgp(u);
    du_dx_gp = triangle.ComputeLinearDUgp(0, u);
    du_dy_gp = triangle.ComputeLinearDUgp(1, u);

    if (!almost_equal(0.866025403784442, triangle.Integration(u_gp)[0], 1.e+04)) {
        error_found = true;

        std::cout << triangle.Integration(u_gp);

        std::cerr << "Error found in Triangle element in ComputeLinearUgp" << std::endl;
    }

    if (!almost_equal(std::sqrt(3.0) / 2.0, triangle.Integration(du_dx_gp)[0], 1.e+04)) {
        error_found = true;

        std::cerr << "Error found in Triangle element in ComputeLinearDUgp "
                     "in x direction"
                  << std::endl;
    }

    if (!almost_equal(-3.0 / 2.0, triangle.Integration(du_dy_gp)[0], 1.e+04)) {
        error_found = true;

        std::cerr << "Error found in Triangle element in ComputeLinearDUgp "
                     "in y direction"
                  << std::endl;
    }

    // Check nodals through the same integration
    u_gp     = triangle.ComputeNodalUgp(u);
    du_dx_gp = triangle.ComputeNodalDUgp(0, u);
    du_dy_gp = triangle.ComputeNodalDUgp(1, u);

    if (!almost_equal(0.866025403784442, triangle.Integration(u_gp)[0], 1.e+04)) {
        error_found = true;

        std::cerr << "Error found in Triangle element in ComputeNodalUgp" << std::endl;
    }

    if (!almost_equal(std::sqrt(3.0) / 2.0, triangle.Integration(du_dx_gp)[0], 1.e+04)) {
        error_found = true;

        std::cerr << "Error found in Triangle element in ComputeNodalDUgp "
                     "in x direction"
                  << std::endl;
    }

    if (!almost_equal(-3.0 / 2.0, triangle.Integration(du_dy_gp)[0], 1.e+04)) {
        error_found = true;

        std::cerr << "Error found in Triangle element in ComputeNodalDUgp "
                     "in y direction"
                  << std::endl;
    }

    // Check ComputeUgp and ApplyMinv
    DynMatrix<double> mod_vals(1, triangle.data.get_ndof());
    DynMatrix<double> gp_vals(1, triangle.data.get_ngp_internal());
    DynMatrix<double> solution_vals(1, triangle.data.get_ndof());

    for (uint dof = 0; dof < 66; ++dof) {
        set_constant(mod_vals, 0.0);
        row(mod_vals, 0)[dof] = 1.0;

        gp_vals       = triangle.ComputeUgp(mod_vals);
        solution_vals = triangle.ApplyMinv(mod_vals);

        for (uint doff = 0; doff < 66; ++doff) {
            if (dof == doff) {
                if (!almost_equal((1. / triangle.IntegrationPhi(doff, gp_vals)[0]), solution_vals(0, dof), 1.e+03)) {
                    error_found = true;

                    std::cerr << "Error found in Triangle element in ApplyMinv" << std::endl;
                }
            } else {
                if (!almost_equal(triangle.IntegrationPhi(doff, gp_vals)[0], 0.0)) {
                    error_found = true;

                    std::cerr << "Error found in Triangle element in ComputeUgp" << std::endl;
                }
            }
        }
    }

    // Check ComputeDUgp
    DynMatrix<double> gp_dvals(1, triangle.data.get_ngp_internal());
    set_constant(gp_vals, 1.0);

    for (uint dof = 0; dof < 66; ++dof) {
        set_constant(mod_vals, 0.0);
        row(mod_vals, 0)[dof] = 1.0;

        gp_dvals = triangle.ComputeDUgp(GlobalCoord::x, mod_vals);

        if (!almost_equal(triangle.IntegrationPhi(0, gp_dvals)[0],
                          triangle.IntegrationDPhi(GlobalCoord::x, dof, gp_vals)[0])) {
            error_found = true;

            std::cerr << "Error found in Triangle element in ComputeDUgp in x "
                         "direction"
                      << std::endl;
        }

        gp_dvals = triangle.ComputeDUgp(GlobalCoord::y, mod_vals);

        if (!almost_equal(triangle.IntegrationPhi(0, gp_dvals)[0],
                          triangle.IntegrationDPhi(GlobalCoord::y, dof, gp_vals)[0])) {
            error_found = true;

            std::cerr << "Error found in Triangle element in ComputeDUgp in y "
                         "direction"
                      << std::endl;
        }
    }

    // Check L2 projection
    DynMatrix<double> nodal_vals(1, 3);
    nodal_vals(0, 0) = 1.0;
    nodal_vals(0, 1) = 2.0;
    nodal_vals(0, 2) = 3.0;

    DynMatrix<double> modal_vals_true(1, 3);
    modal_vals_true(0, 0) = 2.0;
    modal_vals_true(0, 1) = 0.5;
    modal_vals_true(0, 2) = 0.5;

    DynMatrix<double> modal_vals_computed(1, triangle.data.get_ndof());

    modal_vals_computed = triangle.L2ProjectionNode(nodal_vals);

    for (uint i = 0; i < 3; ++i) {
        if (!almost_equal(modal_vals_computed(0, i), modal_vals_true(0, i))) {
            error_found = true;

            std::cerr << "Error found in Triangle element in L2 projection" << std::endl;
        }
    }

    for (uint i = 3; i < 66; ++i) {
        if (!almost_equal(modal_vals_computed(0, i), 0.0, 1.e+04)) {
            error_found = true;

            std::cerr << "Error found in Triangle element in L2 projection" << std::endl;
        }
    }

    // Check integrations PhiPhi PhiDPhi
    DynMatrix<double> unit_gp(1, triangle.data.get_ngp_internal());
    set_constant(unit_gp, 1.0);

    for (uint dof = 0; dof < 66; ++dof) {
        set_constant(mod_vals, 0.0);
        row(mod_vals, 0)[dof] = 1.0;

        gp_vals = triangle.ComputeUgp(mod_vals);

        for (uint doff = 0; doff < 66; ++doff) {
            if (!almost_equal(triangle.IntegrationPhi(doff, gp_vals)[0],
                              triangle.IntegrationPhiPhi(dof, doff, unit_gp)[0],
                              1.e+03)) {
                error_found = true;

                std::cerr << "Error found in Triangle element in IntegrationPhiPhi" << std::endl;
                std::cerr << std::setprecision(16) << "  " << triangle.IntegrationPhi(doff, gp_vals)[0] << ' '
                          << triangle.IntegrationPhiPhi(dof, doff, unit_gp)[0] << '\n';
            }
        }

        gp_dvals = triangle.ComputeDUgp(GlobalCoord::x, mod_vals);

        for (uint doff = 0; doff < 66; ++doff) {
            if (!almost_equal(triangle.IntegrationDPhi(GlobalCoord::x, doff, gp_vals)[0],
                              triangle.IntegrationPhiDPhi(dof, GlobalCoord::x, doff, unit_gp)[0],
                              1.e+03)) {
                error_found = true;

                std::cerr << "Error found in Triangle element in IntegrationPhiDPhi" << std::endl;
                std::cerr << std::setprecision(16) << "  " << triangle.IntegrationDPhi(GlobalCoord::x, doff, gp_vals)[0]
                          << ' ' << triangle.IntegrationPhiDPhi(dof, GlobalCoord::x, doff, unit_gp)[0] << '\n';
            }

            if (!almost_equal(triangle.IntegrationPhi(doff, gp_dvals)[0],
                              triangle.IntegrationPhiDPhi(doff, GlobalCoord::x, dof, unit_gp)[0],
                              1.e+05)) {
                error_found = true;

                std::cerr << "Error found in Triangle element in IntegrationPhiDPhi" << std::endl;
            }
        }

        gp_dvals = triangle.ComputeDUgp(GlobalCoord::y, mod_vals);

        for (uint doff = 0; doff < 66; ++doff) {
            if (!almost_equal(triangle.IntegrationDPhi(GlobalCoord::y, doff, gp_vals)[0],
                              triangle.IntegrationPhiDPhi(dof, GlobalCoord::y, doff, unit_gp)[0],
                              1.e+03)) {
                error_found = true;

                std::cerr << "Error found in Triangle element in IntegrationPhiDPhi" << std::endl;
            }

            if (!almost_equal(triangle.IntegrationPhi(doff, gp_dvals)[0],
                              triangle.IntegrationPhiDPhi(doff, GlobalCoord::y, dof, unit_gp)[0],
                              1.e+05)) {
                error_found = true;

                std::cerr << "Error found in Triangle element in IntegrationPhiDPhi" << std::endl;
            }
        }
    }

    return error_found;
}

#endif