#ifndef SWE_MANUFACTURED_PROBLEM_HPP
#define SWE_MANUFACTURED_PROBLEM_HPP

#include "general_definitions.hpp"

#include "problem/SWE/swe_cuda_problem.hpp"
#include "problem/SWE/swe_kernels_preprocessor.hpp"
#include "problem/SWE/swe_cuda_kernels_processor.hpp"
#include "problem/SWE/swe_kernels_postprocessor.hpp"

#include "swe_true_src_functions.hpp"

namespace SWE {
struct ManufacturedProblem : CUDAProblem {
    template <typename ElementType>
    static void source_kernel(const Stepper& stepper, ElementType& elt);
};

template <typename ElementType>
void ManufacturedProblem::source_kernel(const Stepper& stepper, ElementType& elt) {
    const uint stage = stepper.get_stage();

    auto& state = elt.data.state[stage];
    auto& internal = elt.data.internal;

    double t = stepper.get_t_at_curr_stage();

    const auto source_ze = [t](Point<2>& pt) { return SWE::source_ze(t, pt); };

    const auto source_qx = [t](Point<2>& pt) { return SWE::source_qx(t, pt); };

    const auto source_qy = [t](Point<2>& pt) { return SWE::source_qy(t, pt); };

    elt.ComputeFgp(source_ze, internal.ze_source_term_at_gp);
    elt.ComputeFgp(source_qx, internal.qx_source_term_at_gp);
    elt.ComputeFgp(source_qy, internal.qy_source_term_at_gp);

    // note we assume that the values at gauss points have already been computed
    for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
        // compute contribution of hydrostatic pressure
        internal.qx_source_term_at_gp[gp] += Global::g * internal.bath_deriv_wrt_x_at_gp[gp] * internal.ze_at_gp[gp];
        internal.qy_source_term_at_gp[gp] += Global::g * internal.bath_deriv_wrt_y_at_gp[gp] * internal.ze_at_gp[gp];
    }

    for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
        state.rhs_ze[dof] += elt.IntegrationPhi(dof, internal.ze_source_term_at_gp);
        state.rhs_qx[dof] += elt.IntegrationPhi(dof, internal.qx_source_term_at_gp);
        state.rhs_qy[dof] += elt.IntegrationPhi(dof, internal.qy_source_term_at_gp);
    }
}
}

#endif
