#ifndef RKDG_SWE_PROC_INTFACE_HPP
#define RKDG_SWE_PROC_INTFACE_HPP

namespace SWE {
namespace RKDG {
namespace detail {

template <typename IntfaceSoA>
struct ComputeUgpBdryHelper {
    constexpr bool is_vectorized() { return true; }

    ComputeUgpBdryHelper(uint stage_, const IntfaceSoA& interface_soa_, uint element_type_index_)
        : stage(stage_), interface_soa(interface_soa_) , element_type_index(element_type_index_)
    {}

    template <typename ElementSoA>
    void operator() (ElementSoA& soa) const {
        const auto& state = soa.data.state[stage].q;

        for ( uint side = 0; side < 3; ++side ) {
            const DynMatrix<double>& phi_gp_bdry = interface_soa.GetPhiGPBdry(element_type_index, side);

            for ( uint var =0; var < SWE::n_variables; ++var ) {
                soa.data.boundary[side].q_at_gp = soa.data.state[stage].q * phi_gp_bdry;
            }
        }
    }

private:
    const uint stage;

    const uint element_type_index;
    const IntfaceSoA& interface_soa;
};
}

class InterfaceKernel {
private:
    using StepperType = typename Problem::ProblemStepperType;

public:
    InterfaceKernel(StepperType& stepper_) : stepper(stepper_) {}

    template <typename InterfaceType>
    constexpr static bool is_vectorized() {
        return false;
//        return std::is_same<typename InterfaceType::specialization_t,
        //                          typename ISP::Internal>::value;
    }

    template <typename InterfaceType>
    typename std::enable_if<!is_vectorized<InterfaceType>()>::type operator()(InterfaceType& intface) const {
//        auto& wd_state_in = intface.data_in.wet_dry_state;
//        auto& wd_state_ex = intface.data_ex.wet_dry_state;

//    if (wd_state_in.wet || wd_state_ex.wet) {
        const uint stage = stepper.GetStage();

        auto& state_in    = intface.data_in.state[stage];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

        auto& state_ex    = intface.data_ex.state[stage];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        boundary_in.q_at_gp[SWE::Variables::ze] = intface.ComputeUgpIN(state_in.q[SWE::Variables::ze]);
        boundary_in.q_at_gp[SWE::Variables::qx] = intface.ComputeUgpIN(state_in.q[SWE::Variables::qx]);
        boundary_in.q_at_gp[SWE::Variables::qy] = intface.ComputeUgpIN(state_in.q[SWE::Variables::qy]);

        boundary_ex.q_at_gp[SWE::Variables::ze] = intface.ComputeUgpEX(state_ex.q[SWE::Variables::ze]);
        boundary_ex.q_at_gp[SWE::Variables::qx] = intface.ComputeUgpEX(state_ex.q[SWE::Variables::qx]);
        boundary_ex.q_at_gp[SWE::Variables::qy] = intface.ComputeUgpEX(state_ex.q[SWE::Variables::qy]);

        intface.specialization.ComputeFlux(intface);

        // now compute contributions to the righthand side
        state_in.rhs[SWE::Variables::ze] -= intface.IntegrationPhiIN(boundary_in.F_hat_at_gp[SWE::Variables::ze]);
        state_in.rhs[SWE::Variables::qx] -= intface.IntegrationPhiIN(boundary_in.F_hat_at_gp[SWE::Variables::qx]);
        state_in.rhs[SWE::Variables::qy] -= intface.IntegrationPhiIN(boundary_in.F_hat_at_gp[SWE::Variables::qy]);

        state_ex.rhs[SWE::Variables::ze] -= intface.IntegrationPhiEX(boundary_ex.F_hat_at_gp[SWE::Variables::ze]);
        state_ex.rhs[SWE::Variables::qx] -= intface.IntegrationPhiEX(boundary_ex.F_hat_at_gp[SWE::Variables::qx]);
        state_ex.rhs[SWE::Variables::qy] -= intface.IntegrationPhiEX(boundary_ex.F_hat_at_gp[SWE::Variables::qy]);

    }

    template <typename InterfaceType>
    typename std::enable_if<is_vectorized<InterfaceType>()>::type operator()(InterfaceType& soa) const {
        const uint stage = stepper.GetStage();

        using ElementContainers = typename InterfaceType::ElementContainers;
        ElementContainers& elt_data = soa.GetElementData();

        set_constant(soa.data.Ugp_in, 0.0);
        set_constant(soa.data.Ugp_ex, 0.0);

        uint element_type_index = 0u;
        Utilities::for_each_in_tuple(elt_data,[stage, &soa, &element_type_index](auto& elt_container) {
                //fixme
                //only compute boundary values once
                if ( std::is_same<typename InterfaceType::specialization_t, typename ISP::Internal>::value ) {
                    detail::ComputeUgpBdryHelper<InterfaceType> compute_ugp_bdry_fctor(stage, soa, element_type_index);

                    using is_vectorized_t = std::true_type;
                    elt_container.CallForEachElement(compute_ugp_bdry_fctor, is_vectorized_t{} );
                }

                //Assemble values at GP
                soa.data.Ugp_in += soa.gather_in(element_type_index);
                //flip gp order
                soa.data.Ugp_ex += soa.gather_ex(element_type_index);

                ++element_type_index;
            });

        //soa.specialization.ComputeFlux(soa);

        element_type_index = 0;
        Utilities::for_each_in_tuple(elt_data,[stage, &soa, &element_type_index](auto& elt_container) {
                auto& state = elt_container.data.state[stage];

                auto F_hat_ze_in = soa.scatter_in(element_type_index);
                auto F_hat_qx_in = soa.scatter_in(element_type_index);
                auto F_hat_qy_in = soa.scatter_in(element_type_index);

                state.rhs[SWE::Variables::ze] -= soa.IntegrationPhiIN(F_hat_ze_in);
                state.rhs[SWE::Variables::qx] -= soa.IntegrationPhiIN(F_hat_qx_in);
                state.rhs[SWE::Variables::qy] -= soa.IntegrationPhiIN(F_hat_qy_in);

                auto F_hat_ze_ex = soa.scatter_ex(element_type_index);
                auto F_hat_qx_ex = soa.scatter_ex(element_type_index);
                auto F_hat_qy_ex = soa.scatter_ex(element_type_index);

                state.rhs[SWE::Variables::ze] -= soa.IntegrationPhiEX(F_hat_ze_ex);
                state.rhs[SWE::Variables::qx] -= soa.IntegrationPhiEX(F_hat_qx_ex);
                state.rhs[SWE::Variables::qy] -= soa.IntegrationPhiEX(F_hat_qy_ex);
            });
    }

private:

    StepperType& stepper;
};
}
}

#endif
