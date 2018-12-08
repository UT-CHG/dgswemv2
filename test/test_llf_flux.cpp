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

    { // pass using doubles
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
            std::cerr << "Error in configuration 1 in y-momentum flux\n";
            std::cerr << "Got: " << F_hat[2] << " Should be:  " << true_qy_flux << "\n";
            error_found = true;
        }
    }

    { //pass using input array
        DynVector<double> ze_in_(1), qx_in_(1), qy_in_(1), ze_ex_(1), qx_ex_(1), qy_ex_(1);
        ze_in_[0] = ze_in;
        qx_in_[0] = qx_in;
        qy_in_[0] = qy_in;
        ze_ex_[0] = ze_ex;
        qx_ex_[0] = qx_ex;
        qy_ex_[0] = qy_ex;

        std::array<DynVector<double>,3> aux_in;
        aux_in.fill(DynVector<double>(1));
        aux_in[0][0] = bath;
        aux_in[1][0] = ze_in + bath;
        aux_in[2][0] = sp;

        std::array<DynVector<double>,2> norm;
        norm.fill(DynVector<double>(1));
        norm[0][0] = normal[0];
        norm[1][0] = normal[1];

        std::array<DynVector<double>,3> F_hat;
        F_hat.fill(DynVector<double>(1));

        SWE::RKDG::LLF_flux( SWE::Global::g,
                             ze_in_,
                             qx_in_,
                             qy_in_,
                             ze_ex_,
                             qx_ex_,
                             qy_ex_,
                             aux_in,
                             norm,
                             F_hat);

        if (!Utilities::almost_equal(F_hat[0][0], true_ze_flux)) {
            std::cerr << "Error in configuration " << configuration << " in surface elevation flux\n";
            std::cerr << "Got: " << F_hat[0][0] << " Should be:  " << true_ze_flux << "\n";
            error_found = true;
        }

        if (!Utilities::almost_equal(F_hat[1][0], true_qx_flux)) {
            std::cerr << "Error in configuration " << configuration << " in x-momentum flux\n";
            std::cerr << "Got: " << F_hat[1][0] << " Should be:  " << true_qx_flux << "\n";
            error_found = true;
        }

        if (!Utilities::almost_equal(F_hat[2][0], true_qy_flux)) {
            std::cerr << "Error in configuration 1 in y-momentum flux\n";
            std::cerr << "Got: " << F_hat[2][0] << " Should be:  " << true_qy_flux << "\n";
            error_found = true;
        }

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

        StatVector<double, 2> normal{1. / std::sqrt(2.), 1. / std::sqrt(2.)};

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

        StatVector<double, 2> normal{1. / std::sqrt(2.), -1. / std::sqrt(2.)};

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

        StatVector<double, 2> normal{1. / std::sqrt(2.), 1. / std::sqrt(2.)};

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
