#ifndef RKDG_SWE_PROC_INTFACE_HPP
#define RKDG_SWE_PROC_INTFACE_HPP

namespace SWE {
namespace RKDG {
namespace detail {
template <typename IntfaceSoA>
struct ComputeUgpBdryHelper {
    constexpr bool is_vectorized() { return true; }

    ComputeUgpBdryHelper(uint stage_, uint element_type_index_, IntfaceSoA& interface_soa_)
        : stage(stage_), element_type_index(element_type_index_), interface_soa(interface_soa_)
    {}

    template <typename ElementSoA>
    void operator() (ElementSoA& soa) const {
        const auto& q = soa.data.state[stage].q;

        for ( uint side = 0; side < 3; ++side ) {
            for ( uint var =0; var < SWE::n_variables; ++var ) {
                auto& q_at_gp = soa.data.boundary[side].q_at_gp[var];

                q_at_gp = interface_soa.ComputeUgpBdry(q[var], side, element_type_index);
                interface_soa.data.q_in_at_gp[var] += interface_soa.ScatterIn(element_type_index, side, q_at_gp);
                interface_soa.data.q_ex_at_gp[var] += interface_soa.ScatterEx(element_type_index, side, q_at_gp);
            }
        }
    }

private:
    const uint stage;
    const uint element_type_index;

    IntfaceSoA& interface_soa;
};

template <typename IntfaceSoA>
struct IntegratePhiBdryHelper {
    constexpr bool is_vectorized() { return true; }

    IntegratePhiBdryHelper(uint stage_, uint element_type_index_, IntfaceSoA& interface_soa_)
        : stage(stage_), element_type_index(element_type_index_), interface_soa(interface_soa_)
    {}

    template <typename ElementSoA>
    void operator() (ElementSoA& soa) const {
        auto& state = soa.data.state[stage];

        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            for ( uint side = 0; side < 3; ++side ) {
                soa.data.boundary[side].F_hat_at_gp[var] =
                    interface_soa.GatherIn(element_type_index, side, interface_soa.data.F_hat_at_gp[var]) -
                    interface_soa.GatherEx(element_type_index, side, interface_soa.data.F_hat_at_gp[var]);

                //Not every interface will contribute to this integral, however, the jacobians for these
                //interfaces are set to zero. Therefore, we will effectively mask these contributions to
                //the rhs.
                state.rhs[var] -= interface_soa.IntegratePhiBdry(soa.data.boundary[side].F_hat_at_gp[var],
                                                   side, element_type_index);
            }
        }
    }

private:
    const uint stage;
    const uint element_type_index;

    IntfaceSoA& interface_soa;
};
}

class InterfaceKernel {
private:
    using StepperType = typename Problem::ProblemStepperType;

public:
    InterfaceKernel(StepperType& stepper_) : stepper(stepper_) {}

    template <typename InterfaceType>
    constexpr static bool is_vectorized() {
        return std::is_same<typename InterfaceType::specialization_t,
                            typename ISP::Internal>::value;
    }

    template <typename InterfaceType>
    typename std::enable_if<!is_vectorized<InterfaceType>()>::type operator()(InterfaceType& intface) const {
//        auto& wd_state_in = intface.data_in.wet_dry_state;
//        auto& wd_state_ex = intface.data_ex.wet_dry_state;

//    if (wd_state_in.wet || wd_state_ex.wet) {

        auto& state_in    = intface.data_in.state[stepper.GetStage()];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

        auto& state_ex    = intface.data_ex.state[stepper.GetStage()];
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

        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            set_constant(soa.data.q_in_at_gp[var], 0.0);
            set_constant(soa.data.q_ex_at_gp[var], 0.0);
        }

        Utilities::for_each_in_tuple(elt_data,[stage, &soa](auto& elt_container) {
                constexpr uint element_type_index = Utilities::index<typename std::decay<decltype(elt_container)>::type,
                                                                     ElementContainers>::value;

                //fixme
                //only compute boundary values once
                //if ( std::is_same<typename InterfaceType::specialization_t, typename ISP::Internal>::value ) {
                detail::ComputeUgpBdryHelper<InterfaceType> compute_ugp_bdry_fctor(stage, element_type_index, soa);

                using is_vectorized_t = std::true_type;
                elt_container.CallForEachElement(compute_ugp_bdry_fctor, is_vectorized_t{} );
                //}
            });

        typename InterfaceType::specialization_t spec{};
        spec.ComputeFlux(soa);

        Utilities::for_each_in_tuple(elt_data,[stage, &soa, this](auto& elt_container) {
                constexpr uint element_type_index = Utilities::index<typename std::decay<decltype(elt_container)>::type,
                                                                     ElementContainers>::value;

                detail::IntegratePhiBdryHelper<InterfaceType> integrate_phi_bdry_fctor(stage, element_type_index, soa);
                //fixme
                //Same thing, once you have F_hat_at_gp only compute integrals once at the very end of the vectorized
                //interfaces
                using is_vectorized_t = std::true_type;
                elt_container.CallForEachElement(integrate_phi_bdry_fctor, is_vectorized_t{} );
            });
    }

private:

    StepperType& stepper;
};
}
}

#endif
