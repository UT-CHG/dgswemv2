#ifndef IHDG_SWE_PROC_VOLUME_HPP
#define IHDG_SWE_PROC_VOLUME_HPP

namespace SWE {
namespace IHDG {
template <typename ElementType>
void Problem::prepare_volume_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state_back = elt.data.state[stage];
    auto& state      = elt.data.state[stage + 1];
    auto& internal   = elt.data.internal;
    auto& local      = elt.data.local;

    elt.ComputeUgp(state_back.ze, internal.ze_back_at_gp);
    elt.ComputeUgp(state_back.qx, internal.qx_back_at_gp);
    elt.ComputeUgp(state_back.qy, internal.qy_back_at_gp);

    elt.ComputeUgp(state.ze, internal.ze_at_gp);
    elt.ComputeUgp(state.qx, internal.qx_at_gp);
    elt.ComputeUgp(state.qy, internal.qy_at_gp);

    double u_at_gp, v_at_gp;
    double uuh_at_gp, vvh_at_gp, uvh_at_gp, pe_at_gp;

    for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
        internal.aux_at_gp[gp][SWE::Auxiliaries::h] =
            internal.ze_at_gp[gp] + internal.aux_at_gp[gp][SWE::Auxiliaries::bath];

        u_at_gp = internal.qx_at_gp[gp] / internal.aux_at_gp[gp][SWE::Auxiliaries::h];
        v_at_gp = internal.qy_at_gp[gp] / internal.aux_at_gp[gp][SWE::Auxiliaries::h];

        uuh_at_gp = u_at_gp * internal.qx_at_gp[gp];
        vvh_at_gp = v_at_gp * internal.qy_at_gp[gp];
        uvh_at_gp = u_at_gp * internal.qy_at_gp[gp];
        pe_at_gp  = Global::g * (0.5 * std::pow(internal.ze_at_gp[gp], 2) +
                                internal.ze_at_gp[gp] * internal.aux_at_gp[gp][SWE::Auxiliaries::bath]);

        local.delta_ze_kernel_at_gp[gp] = 1.0 / stepper.GetDT();  // there is a good reason this should be set only once
        local.delta_qx_kernel_at_gp[gp] = 1.0 / stepper.GetDT();  // there is a good reason this should be set only once
        local.delta_qy_kernel_at_gp[gp] = 1.0 / stepper.GetDT();  // there is a good reason this should be set only once

        /* assuming default initialization is 0.0 */

        // local.delta_ze_grad_x_kernel_at_gp[Variables::ze][gp] = 0.0;
        // local.delta_ze_grad_y_kernel_at_gp[Variables::ze][gp] = 0.0;
        local.delta_ze_grad_x_kernel_at_gp[Variables::qx][gp] =
            -std::pow(u_at_gp, 2) + Global::g * internal.aux_at_gp[gp][SWE::Auxiliaries::h];
        local.delta_ze_grad_y_kernel_at_gp[Variables::qx][gp] = -u_at_gp * v_at_gp;
        local.delta_ze_grad_x_kernel_at_gp[Variables::qy][gp] = -u_at_gp * v_at_gp;
        local.delta_ze_grad_y_kernel_at_gp[Variables::qy][gp] =
            -std::pow(v_at_gp, 2) + Global::g * internal.aux_at_gp[gp][SWE::Auxiliaries::h];

        /* there is a good reason this should be set only once */
        local.delta_qx_grad_x_kernel_at_gp[Variables::ze][gp] =
            1.0;  // there is a good reason this should be set only once
        // local.delta_qx_grad_y_kernel_at_gp[Variables::ze][gp] = 0.0;
        local.delta_qx_grad_x_kernel_at_gp[Variables::qx][gp] = 2 * u_at_gp;
        local.delta_qx_grad_y_kernel_at_gp[Variables::qx][gp] = v_at_gp;
        local.delta_qx_grad_x_kernel_at_gp[Variables::qy][gp] = v_at_gp;
        // local.delta_qx_grad_y_kernel_at_gp[Variables::qy][gp] = 0.0;

        // local.delta_qy_grad_x_kernel_at_gp[Variables::ze][gp] = 0.0;
        local.delta_qy_grad_y_kernel_at_gp[Variables::ze][gp] =
            1.0;  // there is a good reason this should be set only once
        // local.delta_qy_grad_x_kernel_at_gp[Variables::qx][gp] = 0.0;
        local.delta_qy_grad_y_kernel_at_gp[Variables::qx][gp] = u_at_gp;
        local.delta_qy_grad_x_kernel_at_gp[Variables::qy][gp] = u_at_gp;
        local.delta_qy_grad_y_kernel_at_gp[Variables::qy][gp] = 2 * v_at_gp;

        local.ze_rhs_kernel_at_gp[gp] = (internal.ze_back_at_gp[gp] - internal.ze_at_gp[gp]) / stepper.GetDT();
        local.qx_rhs_kernel_at_gp[gp] = (internal.qx_back_at_gp[gp] - internal.qx_at_gp[gp]) / stepper.GetDT();
        local.qy_rhs_kernel_at_gp[gp] = (internal.qy_back_at_gp[gp] - internal.qy_at_gp[gp]) / stepper.GetDT();

        local.ze_rhs_grad_x_kernel_at_gp[gp] = internal.qx_at_gp[gp];
        local.ze_rhs_grad_y_kernel_at_gp[gp] = internal.qy_at_gp[gp];

        local.qx_rhs_grad_x_kernel_at_gp[gp] = (uuh_at_gp + pe_at_gp);
        local.qx_rhs_grad_y_kernel_at_gp[gp] = uvh_at_gp;

        local.qy_rhs_grad_x_kernel_at_gp[gp] = uvh_at_gp;
        local.qy_rhs_grad_y_kernel_at_gp[gp] = vvh_at_gp + pe_at_gp;
    }
}
}
}

#endif
