#include "general_definitions.hpp"
#include "utilities/almost_equal.hpp"
#include "problem/SWE/swe_definitions.hpp"

#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

#include <iostream>

bool test_configuration(const int configuration,
                        const double ze_in,
                        const double ze_ex,
                        const double qx_in,
                        const double qx_ex,
                        const double qy_in,
                        const double qy_ex,
                        const double bath,
                        const double sp,
                        const StatVector<double, 2>& normal,
                        const double true_ze_flux,
                        const double true_qx_flux,
                        const double true_qy_flux) {
    bool error_found = false;

    // pass using doubles
    std::array<double,3> aux_in{bath, ze_in + bath, sp};
    std::array<double,2> norm{normal[0], normal[1]};

    std::array<double,3> F_hat;
    SWE::RKDG::LLF_flux( SWE::Global::g,
                         ze_in,
                         qx_in,
                         qy_in,
                         ze_ex,
                         qx_ex,
                         qy_ex,
                         aux_in,
                         norm,
                         F_hat);

    if (!Utilities::almost_equal(F_hat[0], true_ze_flux)) {
        std::cerr << "Error in configuration " << configuration << " in surface elevation flux\n";
        std::cerr << "Got: " << F_hat[0] << " Should be:  " << true_ze_flux << "\n";
        error_found = true;
    }

    if (!Utilities::almost_equal(F_hat[1], true_qx_flux)) {
        std::cerr << "Error in configuration " << configuration << " in x-momentum flux\n";
        std::cerr << "Got: " << F_hat[1] << " Should be:  " << true_qx_flux << "\n";
        error_found = true;
    }

    if (!Utilities::almost_equal(F_hat[2], true_qy_flux)) {
        std::cerr << "Error in configuration " << configuration << " in y-momentum flux\n";
        std::cerr << "Got: " << F_hat[2] << " Should be:  " << true_qy_flux << "\n";
        error_found = true;
    }

    return error_found;
}

bool test_configuration_vectorized(const DynRowVector<double>& ze_in,
                                   const DynRowVector<double>& ze_ex,
                                   const DynRowVector<double>& qx_in,
                                   const DynRowVector<double>& qx_ex,
                                   const DynRowVector<double>& qy_in,
                                   const DynRowVector<double>& qy_ex,
                                   const std::array<DynRowVector<double>, SWE::n_auxiliaries>& aux,
                                   const std::array<DynRowVector<double>, SWE::n_dimensions>& surface_normal,
                                   const DynRowVector<double>& true_ze_flux,
                                   const DynRowVector<double>& true_qx_flux,
                                   const DynRowVector<double>& true_qy_flux) {
    bool error_found = false;

    std::array<DynRowVector<double>,3> F_hat;
    F_hat.fill(DynRowVector<double>(ze_in.size()));
    SWE::RKDG::LLF_flux( SWE::Global::g,
                         ze_in,
                         qx_in,
                         qy_in,
                         ze_ex,
                         qx_ex,
                         qy_ex,
                         aux,
                         surface_normal,
                         F_hat);

    double ze_norm = norm( F_hat[SWE::Variables::ze] - true_ze_flux);
    if (!Utilities::almost_equal(0,ze_norm)) {
        std::cerr << "Error in vectorized configuration in surface elevation flux\n";
        error_found = true;
    }

    double qx_norm = norm( F_hat[SWE::Variables::qx] - true_qx_flux);
    if (!Utilities::almost_equal(0, qx_norm)) {
        std::cerr << "Error in vectorized configuration in x-momentum flux\n";
        error_found = true;
    }

    double qy_norm = norm( F_hat[SWE::Variables::qy] - true_qy_flux);
    if (!Utilities::almost_equal(0,qy_norm)) {
        std::cerr << "Error in vectorized_configuration in y-momentum flux\n";
        error_found = true;
    }

    return error_found;

}


int main() {
    using Utilities::almost_equal;

    if (!almost_equal(SWE::Global::g, 9.81)) {
        std::cerr << "Error: this test was designed for a static gravity of 9.81\n";
        return 1;
    }

    bool error_found = false;

    DynRowVector<double> ze_in(3), ze_ex(3), qx_in(3), qx_ex(3), qy_in(3), qy_ex(3);

    ze_in[0] = 1.0;
    ze_in[1] = 1.0;
    ze_in[2] = 1.0;

    ze_ex[0] = 1.5;
    ze_ex[1] = 1.5;
    ze_ex[2] = 1.5;

    qx_in[0] = 2.2;
    qx_in[1] = 0;
    qx_in[2] = 5.2;

    qx_ex[0] = 3.2;
    qx_ex[1] = 0;
    qx_ex[2] = 6.2;

    qy_in[0] = 0.1;
    qy_in[1] = 0;
    qy_in[2] = 7.1;

    qy_ex[0] = -0.2;
    qy_ex[1] = 0;
    qy_ex[2] = 6.2;


    std::array<DynRowVector<double>, SWE::n_auxiliaries> aux;
    aux.fill(DynRowVector<double>(3));
    set_constant(aux[SWE::Auxiliaries::bath], 0);
    set_constant(aux[SWE::Auxiliaries::sp], 1);

    std::array<DynRowVector<double>, SWE::n_dimensions> surface_normal;
    surface_normal.fill(DynRowVector<double>(3));
    surface_normal[GlobalCoord::x][0] = 1. / std::sqrt(2.);
    surface_normal[GlobalCoord::x][1] = 1. / std::sqrt(2.);
    surface_normal[GlobalCoord::x][2] = 1. / std::sqrt(2.);

    surface_normal[GlobalCoord::y][0] = 1. / std::sqrt(2.);
    surface_normal[GlobalCoord::y][1] = -1. / std::sqrt(2.);
    surface_normal[GlobalCoord::y][2] = 1. / std::sqrt(2.);

    DynRowVector<double> true_ze_flux(3), true_qx_flux(3), true_qy_flux(3);
    true_ze_flux[0] = 0.561276190610246;
    true_qx_flux[0] = 7.062691284925730;
    true_qy_flux[0] = 6.363512979114632;

    true_ze_flux[1] = -0.959003388940832;
    true_qx_flux[1] = 5.636082987795024;
    true_qy_flux[1] = -5.636082987795024;

    true_ze_flux[2] = 5.775392407336936;
    true_qx_flux[2] = 40.455394948714215;
    true_qy_flux[2] = 59.955967780083341;

    // we test a few configurations
    for (uint config_id = 0; config_id < 3; ++config_id) {  // Configuration 1

        StatVector<double, 2> normal{surface_normal[GlobalCoord::x][config_id],
                surface_normal[GlobalCoord::y][config_id]};

        if (test_configuration(config_id,
                               ze_in[config_id],
                               ze_ex[config_id],
                               qx_in[config_id],
                               qx_ex[config_id],
                               qy_in[config_id],
                               qy_ex[config_id],
                               aux[SWE::Auxiliaries::bath][config_id],
                               aux[SWE::Auxiliaries::sp][config_id],
                               normal,
                               true_ze_flux[config_id],
                               true_qx_flux[config_id],
                               true_qy_flux[config_id])) {
            error_found = true;
        }
    }

    // Test all configurations at once
    if ( test_configuration_vectorized(ze_in,
                                       ze_ex,
                                       qx_in,
                                       qx_ex,
                                       qy_in,
                                       qy_ex,
                                       aux,
                                       surface_normal,
                                       true_ze_flux,
                                       true_qx_flux,
                                       true_qy_flux)) {
        error_found  = true;
    }

    if (error_found) {
        return 1;
    }
    return 0;
}
