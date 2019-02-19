#ifndef IHDG_SWE_PROC_SOURCE_HPP
#define IHDG_SWE_PROC_SOURCE_HPP

#include "problem/SWE/problem_source/swe_source.hpp"

namespace SWE {
namespace IHDG {
template <typename StepperType, typename ElementType>
void Problem::init_source_kernel(const StepperType& stepper, ElementType& elt) {
    if (stepper.GetOrder() == 2) {
        auto& internal = elt.data.internal;

        SWE::get_source(stepper.GetTimeAtCurrentStage(), elt);

        for (uint dof_i = 0; dof_i < elt.data.get_ndof(); ++dof_i) {
            subvector(internal.rhs_prev, SWE::n_variables * dof_i, SWE::n_variables) +=
                elt.IntegrationPhi(dof_i, internal.source_at_gp);
        }
    }
}

template <typename StepperType, typename ElementType>
void Problem::local_source_kernel(const StepperType& stepper, ElementType& elt) {
    auto& internal = elt.data.internal;
    auto& source   = elt.data.source;

    set_constant(internal.dsource_dq_at_gp, 0.0);

    if (SWE::SourceTerms::bottom_friction) {
        double Cf      = SWE::SourceTerms::Cf;
        double dCf_dze = 0.0;

        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            // compute bottom friction contribution
            double u = internal.q_at_gp(SWE::Variables::qx, gp) / internal.aux_at_gp(SWE::Auxiliaries::h, gp);
            double v = internal.q_at_gp(SWE::Variables::qy, gp) / internal.aux_at_gp(SWE::Auxiliaries::h, gp);

            // compute manning friction factor
            if (source.manning) {
                Cf = source.g_manning_n_sq / std::pow(internal.aux_at_gp(SWE::Auxiliaries::h, gp), 1.0 / 3.0);
                dCf_dze =
                    -source.g_manning_n_sq / std::pow(internal.aux_at_gp(SWE::Auxiliaries::h, gp), 4.0 / 3.0) / 3.0;

                if (Cf < SWE::SourceTerms::Cf) {
                    Cf      = SWE::SourceTerms::Cf;
                    dCf_dze = 0.0;
                }
            }

            double bottom_friction_stress = Cf * std::hypot(u, v) / internal.aux_at_gp(SWE::Auxiliaries::h, gp);
            double C                      = Cf / std::hypot(u, v) / internal.aux_at_gp(SWE::Auxiliaries::h, gp);

            internal.dsource_dq_at_gp(JacobianVariables::qx_ze, gp) -=
                dCf_dze * bottom_friction_stress / Cf * internal.q_at_gp(SWE::Variables::qx, gp) -
                2 * bottom_friction_stress * u;
            internal.dsource_dq_at_gp(JacobianVariables::qx_qx, gp) -= bottom_friction_stress + C * u * u;
            internal.dsource_dq_at_gp(JacobianVariables::qx_qy, gp) -= C * u * v;

            internal.dsource_dq_at_gp(JacobianVariables::qy_ze, gp) -=
                dCf_dze * bottom_friction_stress / Cf * internal.q_at_gp(SWE::Variables::qy, gp) -
                2 * bottom_friction_stress * v;
            internal.dsource_dq_at_gp(JacobianVariables::qy_qx, gp) -= C * u * v;
            internal.dsource_dq_at_gp(JacobianVariables::qy_qy, gp) -= bottom_friction_stress + C * v * v;
        }
    }

    SWE::get_source(stepper.GetTimeAtCurrentStage() + stepper.GetDT(), elt);

    double theta = stepper.GetTheta();

    for (uint dof_i = 0; dof_i < elt.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < elt.data.get_ndof(); ++dof_j) {
            submatrix(internal.delta_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) +=
                reshape<double, SWE::n_variables>((1.0 - theta) *
                                                  elt.IntegrationPhiPhi(dof_j, dof_i, internal.dsource_dq_at_gp));
        }

        subvector(internal.rhs_local, SWE::n_variables * dof_i, SWE::n_variables) +=
            (1.0 - theta) * elt.IntegrationPhi(dof_i, internal.source_at_gp);
    }
}
}
}

#endif
