#include "utilities/almost_equal.hpp"

#include "dg_problem/SWE/numerical_fluxes/swe_numerical_fluxes.hpp"

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
                        std::vector<double>& normal,
                        const double true_ze_flux,
                        const double true_qx_flux,
                        const double true_qy_flux) {
    bool error_found = false;

    double ze_flux, qx_flux, qy_flux;

    SWE::LLF_flux(
        SWE::Global::g, ze_in, ze_ex, qx_in, qx_ex, qy_in, qy_ex, bath, sp, normal, ze_flux, qx_flux, qy_flux);

    if (!Utilities::almost_equal(ze_flux, true_ze_flux)) {
        std::cerr << "Error in configuration " << configuration << " in surface elevation flux\n";
        std::cerr << "Got: " << ze_flux << " Should be:  " << true_ze_flux << "\n";
        error_found = true;
    }

    if (!Utilities::almost_equal(qx_flux, true_qx_flux)) {
        std::cerr << "Error in configuration " << configuration << " in x-momentum flux\n";
        std::cerr << "Got: " << qx_flux << " Should be:  " << true_qx_flux << "\n";
        error_found = true;
    }

    if (!Utilities::almost_equal(qy_flux, true_qy_flux)) {
        std::cerr << "Error in configuration 1 in y-momentum flux\n";
        std::cerr << "Got: " << qy_flux << " Should be:  " << true_qy_flux << "\n";
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

    // we test a few configurations
    {  // Configuration 1
        double ze_in = 1.0;
        double ze_ex = 1.5;
        double qx_in = 2.2;
        double qx_ex = 3.2;
        double qy_in = 0.1;
        double qy_ex = -0.2;

        std::vector<double> normal{1. / std::sqrt(2.), 1. / std::sqrt(2.)};

        double bath = 0;
        double sp   = 1;

        if (test_configuration(1,
                               ze_in,
                               ze_ex,
                               qx_in,
                               qx_ex,
                               qy_in,
                               qy_ex,
                               bath,
                               sp,
                               normal,
                               0.561276190610246,
                               7.062691284925730,
                               6.363512979114632)) {
            error_found = true;
        }
    }

    {  // Configuration 2
        double ze_in = 1.0;
        double ze_ex = 1.5;
        double qx_in = 0;
        double qx_ex = 0;
        double qy_in = 0;
        double qy_ex = 0;

        std::vector<double> normal{1. / std::sqrt(2.), -1. / std::sqrt(2.)};

        double bath = 0;
        double sp   = 1;

        if (test_configuration(2,
                               ze_in,
                               ze_ex,
                               qx_in,
                               qx_ex,
                               qy_in,
                               qy_ex,
                               bath,
                               sp,
                               normal,
                               -0.959003388940832,
                               5.636082987795024,
                               -5.636082987795024)) {
            error_found = true;
        }
    }

    {  // Configuration 3
        double ze_in = 1.0;
        double ze_ex = 1.5;
        double qx_in = 5.2;
        double qx_ex = 6.2;
        double qy_in = 7.1;
        double qy_ex = 6.2;

        std::vector<double> normal{1. / std::sqrt(2.), 1. / std::sqrt(2.)};

        double bath = 0;
        double sp   = 1;

        if (test_configuration(2,
                               ze_in,
                               ze_ex,
                               qx_in,
                               qx_ex,
                               qy_in,
                               qy_ex,
                               bath,
                               sp,
                               normal,
                               5.775392407336936,
                               40.455394948714215,
                               59.955967780083341)) {
            error_found = true;
        }
    }

    if (error_found) {
        return 1;
    }
    return 0;
}
