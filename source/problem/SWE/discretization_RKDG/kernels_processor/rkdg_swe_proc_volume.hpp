#ifndef RKDG_SWE_PROC_VOLUME_HPP
#define RKDG_SWE_PROC_VOLUME_HPP

#include "problem/SWE/problem_flux/swe_flux.hpp"
namespace SWE {
namespace RKDG {
template <typename StepperType>
class VolumeKernel {
public:
    template <typename ElementType>
    constexpr static bool is_vectorized() { return true; }

    VolumeKernel(const StepperType& stepper_) : stepper(stepper_) {}

    template <typename SoA>
    void operator() (SoA& soa) const {
        const uint stage = stepper.GetStage();

//    auto& wd_state = elt.data.wet_dry_state;
        auto& state    = soa.data.state[stage];

//    for ( uint var = 0; var < SWE::n_variables; ++var ) {
//        set_constant(state.rhs[var], 0.0);
//    }

//    if (wd_state.wet) {
        auto& internal = soa.data.internal;

        internal.q_at_gp[SWE::Variables::ze] = soa.ComputeUgp(state.q[SWE::Variables::ze]);
        internal.q_at_gp[SWE::Variables::qx] = soa.ComputeUgp(state.q[SWE::Variables::qx]);
        internal.q_at_gp[SWE::Variables::qy] = soa.ComputeUgp(state.q[SWE::Variables::qy]);

        auto& ze   = internal.q_at_gp[SWE::Variables::ze];
        auto& qx   = internal.q_at_gp[SWE::Variables::qx];
        auto& qy   = internal.q_at_gp[SWE::Variables::qy];
        auto&  h   = internal.aux_at_gp[SWE::Auxiliaries::h];
        auto& bath = internal.aux_at_gp[SWE::Auxiliaries::bath];
        h = ze + bath;

        auto u = mat_cw_div(qx, h);
        auto v = mat_cw_div(qy, h);

        auto uuh = mat_cw_mult(u, qx);
        auto vvh = mat_cw_mult(v, qy);
        auto uvh = mat_cw_mult(u, qy);
        auto pe  = Global::g * mat_cw_mult(ze, 0.5 * ze + bath);

        //F_at_gp is a temporary variable. Since we're only integrating one coordinate at a time,
        //we simply reuse F_at_gp
        //compute ze terms
        internal.F_at_gp[GlobalCoord::x] =
            mat_cw_mult(internal.aux_at_gp[SWE::Auxiliaries::sp], qx);
        internal.F_at_gp[GlobalCoord::y] = qy;

        state.rhs[SWE::Variables::ze] = soa.IntegrationDPhi( internal.F_at_gp );

        //compute qx terms
        internal.F_at_gp[GlobalCoord::x] =
            mat_cw_mult(internal.aux_at_gp[SWE::Auxiliaries::sp], uuh + pe);

        internal.F_at_gp[GlobalCoord::y] = uvh;

        state.rhs[SWE::Variables::qx] = soa.IntegrationDPhi( internal.F_at_gp );

        //compute qy terms

        //Save one expression template evaluation by using the symmetric nature of the flux
        auto& uvh_ = internal.F_at_gp[GlobalCoord::y];
        internal.F_at_gp[GlobalCoord::x] =
            mat_cw_mult(internal.aux_at_gp[SWE::Auxiliaries::sp], uvh_ );

        internal.F_at_gp[GlobalCoord::y] = vvh + pe;

        state.rhs[SWE::Variables::qy] = soa.IntegrationDPhi( internal.F_at_gp );
//    }
    }

private:
    const StepperType& stepper;
};
}
}

#endif
