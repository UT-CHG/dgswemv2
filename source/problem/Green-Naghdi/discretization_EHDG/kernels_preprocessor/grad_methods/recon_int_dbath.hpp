#ifndef RECON_INT_DBATH_HPP
#define RECON_INT_DBATH_HPP

namespace GN {
namespace EHDG {
template <typename ProblemDiscretizationType, typename ProblemGlobalDataType>
void reconstruct_dbath(ProblemDiscretizationType& discretization, ProblemGlobalDataType& global_data) {
    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& derivative_in                                       = intface.data_in.derivative;
        auto& derivative_ex                                       = intface.data_ex.derivative;
        derivative_in.dbath_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.dbath_at_baryctr;
        derivative_ex.dbath_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.dbath_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& derivative                                  = bound.data.derivative;
        derivative.dbath_at_baryctr_neigh[bound.bound_id] = derivative.dbath_at_baryctr;
    });

    discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& derivative = dbound.data.derivative;
        std::vector<double> message(GN::n_dimensions);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dbath, message);
        for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
            derivative.dbath_at_baryctr_neigh[dbound.bound_id][dbath] = message[dbath];
        }
    });

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& state      = elt.data.state[0];
        auto& derivative = elt.data.derivative;
        auto& internal   = elt.data.internal;

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            const uint element_1 = derivative.a_elem[2 * bound];
            const uint element_2 = derivative.a_elem[2 * bound + 1];
            column(derivative.dbath_at_midpts, bound) =
                derivative.dbath_at_baryctr +
                (derivative.dbath_at_baryctr_neigh[element_1] - derivative.dbath_at_baryctr) * derivative.a[bound][0] +
                (derivative.dbath_at_baryctr_neigh[element_2] - derivative.dbath_at_baryctr) * derivative.a[bound][1];
        }
        derivative.dbath_lin = derivative.dbath_at_midpts * derivative.T;
        state.dbath          = elt.ProjectLinearToBasis(derivative.dbath_lin);
        internal.dbath_at_gp = elt.ComputeUgp(state.dbath);
    });
}

template <typename ProblemDiscretizationType, typename ProblemGlobalDataType>
void reconstruct_ddbath(ProblemDiscretizationType& discretization, ProblemGlobalDataType& global_data) {
    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& derivative_in                                        = intface.data_in.derivative;
        auto& derivative_ex                                        = intface.data_ex.derivative;
        derivative_in.ddbath_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.ddbath_at_baryctr;
        derivative_ex.ddbath_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.ddbath_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& derivative                                   = bound.data.derivative;
        derivative.ddbath_at_baryctr_neigh[bound.bound_id] = derivative.ddbath_at_baryctr;
    });

    discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& derivative = dbound.data.derivative;
        std::vector<double> message(GN::n_ddbath_terms);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dbath, message);
        for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
            derivative.ddbath_at_baryctr_neigh[dbound.bound_id][ddbath] = message[ddbath];
        }
    });

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& state      = elt.data.state[0];
        auto& derivative = elt.data.derivative;
        auto& internal   = elt.data.internal;

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            const uint element_1 = derivative.a_elem[2 * bound];
            const uint element_2 = derivative.a_elem[2 * bound + 1];
            column(derivative.ddbath_at_midpts, bound) =
                derivative.ddbath_at_baryctr +
                (derivative.ddbath_at_baryctr_neigh[element_1] - derivative.ddbath_at_baryctr) *
                    derivative.a[bound][0] +
                (derivative.ddbath_at_baryctr_neigh[element_2] - derivative.ddbath_at_baryctr) * derivative.a[bound][1];
        }
        derivative.ddbath_lin = derivative.ddbath_at_midpts * derivative.T;
        state.ddbath          = elt.ProjectLinearToBasis(derivative.ddbath_lin);
        internal.ddbath_at_gp = elt.ComputeUgp(state.ddbath);
    });
}

template <typename ProblemDiscretizationType, typename ProblemGlobalDataType>
void reconstruct_dddbath(ProblemDiscretizationType& discretization, ProblemGlobalDataType& global_data) {
    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& derivative_in                                         = intface.data_in.derivative;
        auto& derivative_ex                                         = intface.data_ex.derivative;
        derivative_in.dddbath_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.dddbath_at_baryctr;
        derivative_ex.dddbath_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.dddbath_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& derivative                                    = bound.data.derivative;
        derivative.dddbath_at_baryctr_neigh[bound.bound_id] = derivative.dddbath_at_baryctr;
    });

    discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& derivative = dbound.data.derivative;
        std::vector<double> message(GN::n_dddbath_terms);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dbath, message);
        for (uint dddbath = 0; dddbath < GN::n_dddbath_terms; ++dddbath) {
            derivative.dddbath_at_baryctr_neigh[dbound.bound_id][dddbath] = message[dddbath];
        }
    });

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& state      = elt.data.state[0];
        auto& derivative = elt.data.derivative;
        auto& internal   = elt.data.internal;

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            const uint element_1 = derivative.a_elem[2 * bound];
            const uint element_2 = derivative.a_elem[2 * bound + 1];
            column(derivative.dddbath_at_midpts, bound) =
                derivative.dddbath_at_baryctr +
                (derivative.dddbath_at_baryctr_neigh[element_1] - derivative.dddbath_at_baryctr) *
                    derivative.a[bound][0] +
                (derivative.dddbath_at_baryctr_neigh[element_2] - derivative.dddbath_at_baryctr) *
                    derivative.a[bound][1];
        }
        derivative.dddbath_lin = derivative.dddbath_at_midpts * derivative.T;
        state.dddbath          = elt.ProjectLinearToBasis(derivative.dddbath_lin);
        internal.dddbath_at_gp = elt.ComputeUgp(state.dddbath);
    });
}
}
}

#endif
