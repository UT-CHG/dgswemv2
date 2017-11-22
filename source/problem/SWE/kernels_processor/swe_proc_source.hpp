#ifndef SWE_PROC_SOURCE_HPP
#define SWE_PROC_SOURCE_HPP

namespace SWE {
template <typename ElementType>
void Problem::source_kernel(const Stepper& stepper, ElementType& elt) {
    auto& wd_state = elt.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.get_stage();

        auto& state = elt.data.state[stage];
        auto& internal = elt.data.internal;

        double t = stepper.get_t_at_curr_stage();

        auto source_ze = [t](Point<2>& pt) { return SWE::source_ze(t, pt); };

        auto source_qx = [t](Point<2>& pt) { return SWE::source_qx(t, pt); };

        auto source_qy = [t](Point<2>& pt) { return SWE::source_qy(t, pt); };

        elt.ComputeFgp(source_ze, internal.ze_source_term_at_gp);
        elt.ComputeFgp(source_qx, internal.qx_source_term_at_gp);
        elt.ComputeFgp(source_qy, internal.qy_source_term_at_gp);

        // note we assume that the values at gauss points have already been computed
        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            // compute contribution of hydrostatic pressure
            internal.qx_source_term_at_gp[gp] +=
                Global::g * internal.bath_deriv_wrt_x_at_gp[gp] * internal.ze_at_gp[gp];
            internal.qy_source_term_at_gp[gp] +=
                Global::g * internal.bath_deriv_wrt_y_at_gp[gp] * internal.ze_at_gp[gp];

            double u_at_gp = internal.qx_at_gp[gp] / internal.h_at_gp[gp];
            double v_at_gp = internal.qy_at_gp[gp] / internal.h_at_gp[gp];

            // compute bottom friction contribution
            double bottom_friction_stress = Global::Cf * std::hypot(u_at_gp, v_at_gp) / internal.h_at_gp[gp];

            internal.qx_source_term_at_gp[gp] -= bottom_friction_stress * internal.qx_at_gp[gp];
            internal.qy_source_term_at_gp[gp] -= bottom_friction_stress * internal.qy_at_gp[gp];
        }

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.rhs_ze[dof] += elt.IntegrationPhi(dof, internal.ze_source_term_at_gp);
            state.rhs_qx[dof] += elt.IntegrationPhi(dof, internal.qx_source_term_at_gp);
            state.rhs_qy[dof] += elt.IntegrationPhi(dof, internal.qy_source_term_at_gp);
        }
    }
}
}

#endif
