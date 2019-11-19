#ifndef LEASTSQUARES_DBATH_SERIAL_HPP
#define LEASTSQUARES_DBATH_SERIAL_HPP

namespace GN {
namespace EHDG {
template <typename ProblemDiscretizationType>
void compute_dbath_ls(ProblemDiscretizationType& discretization);
template <typename ProblemDiscretizationType>
void compute_ddbath_ls(ProblemDiscretizationType& discretization);
template <typename ProblemDiscretizationType>
void compute_dddbath_ls(ProblemDiscretizationType& discretization);

void Problem::compute_bathymetry_derivatives_serial(ProblemDiscretizationType& discretization,
                                                    ProblemGlobalDataType& global_data) {
    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& state                = elt.data.state[0];
        auto& derivative           = elt.data.derivative;
        derivative.bath_lin        = elt.ProjectBasisToLinear(row(state.aux, SWE::Auxiliaries::bath));
        derivative.bath_at_baryctr = elt.ComputeLinearUbaryctr(derivative.bath_lin);
        derivative.bath_at_midpts  = elt.ComputeLinearUmidpts(derivative.bath_lin);
    });
    compute_dbath_ls(discretization);
    reconstruct_dbath(discretization, global_data);

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& derivative            = elt.data.derivative;
        derivative.dbath_at_baryctr = elt.ComputeLinearUbaryctr(derivative.dbath_lin);
        derivative.dbath_at_midpts  = elt.ComputeLinearUmidpts(derivative.dbath_lin);
    });
    compute_ddbath_ls(discretization);
    reconstruct_ddbath(discretization, global_data);

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& derivative             = elt.data.derivative;
        derivative.ddbath_at_baryctr = elt.ComputeLinearUbaryctr(derivative.ddbath_lin);
        derivative.ddbath_at_midpts  = elt.ComputeLinearUmidpts(derivative.ddbath_lin);
    });
    compute_dddbath_ls(discretization);
    reconstruct_dddbath(discretization, global_data);
}

template <typename ProblemDiscretizationType>
void compute_dbath_ls(ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& derivative_in                                      = intface.data_in.derivative;
        auto& derivative_ex                                      = intface.data_ex.derivative;
        derivative_in.bath_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.bath_at_baryctr;
        derivative_ex.bath_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.bath_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& derivative                                 = bound.data.derivative;
        derivative.bath_at_baryctr_neigh[bound.bound_id] = derivative.bath_at_midpts[bound.bound_id];
    });

    discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& derivative = dbound.data.derivative;
        std::vector<double> message(1);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dbath, message);
        derivative.bath_at_baryctr_neigh[dbound.bound_id] = message[0];
    });

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& derivative = elt.data.derivative;
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            derivative.bath_at_midpts[bound] = derivative.bath_at_baryctr_neigh[bound] - derivative.bath_at_baryctr;
        }
        derivative.dbath_at_baryctr = derivative.P * transpose(derivative.bath_at_midpts);
    });
}

template <typename ProblemDiscretizationType>
void compute_ddbath_ls(ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& state_in      = intface.data_in.state[0];
        auto& state_ex      = intface.data_ex.state[0];
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;
        auto& boundary_in   = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex   = intface.data_ex.boundary[intface.bound_id_ex];

        boundary_in.dbath_hat_at_gp = intface.ComputeUgpIN(state_in.dbath);
        boundary_ex.dbath_hat_at_gp = intface.ComputeUgpEX(state_ex.dbath);

        const uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        for (uint gp = 0; gp < ngp; ++gp) {
            const uint gp_ex = ngp - gp - 1;
            for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                boundary_in.dbath_hat_at_gp(dbath, gp) =
                    (boundary_in.dbath_hat_at_gp(dbath, gp) + boundary_ex.dbath_hat_at_gp(dbath, gp_ex)) / 2.0;
                boundary_ex.dbath_hat_at_gp(dbath, gp_ex) = boundary_in.dbath_hat_at_gp(dbath, gp);
            }
        }

        derivative_in.dbath_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.dbath_at_baryctr;
        derivative_ex.dbath_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.dbath_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& state      = bound.data.state[0];
        auto& derivative = bound.data.derivative;
        auto& boundary   = bound.data.boundary[bound.bound_id];

        boundary.dbath_hat_at_gp = bound.ComputeUgp(state.dbath);

        derivative.dbath_at_baryctr_neigh[bound.bound_id] = column(derivative.dbath_at_midpts, bound.bound_id);
    });

    discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& state      = dbound.data.state[0];
        auto& derivative = dbound.data.derivative;
        auto& boundary   = dbound.data.boundary[dbound.bound_id];

        boundary.dbath_hat_at_gp = dbound.ComputeUgp(state.dbath);

        std::vector<double> message(GN::n_dimensions);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dbath, message);
        for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
            derivative.dbath_at_baryctr_neigh[dbound.bound_id][dbath] = message[dbath];
        }
    });

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& derivative = elt.data.derivative;
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            column(derivative.dbath_at_midpts, bound) =
                derivative.dbath_at_baryctr_neigh[bound] - derivative.dbath_at_baryctr;
        }
        derivative.ddbath_at_baryctr = flatten<double, GN::n_dimensions, GN::n_dimensions, SO::ColumnMajor>(
            derivative.P * transpose(derivative.dbath_at_midpts));
    });
}

template <typename ProblemDiscretizationType>
void compute_dddbath_ls(ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& derivative_in                                        = intface.data_in.derivative;
        auto& derivative_ex                                        = intface.data_ex.derivative;
        derivative_in.ddbath_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.ddbath_at_baryctr;
        derivative_ex.ddbath_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.ddbath_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& derivative                                   = bound.data.derivative;
        derivative.ddbath_at_baryctr_neigh[bound.bound_id] = column(derivative.ddbath_at_midpts, bound.bound_id);
    });

    discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& derivative = dbound.data.derivative;
        auto& boundary   = dbound.data.boundary[dbound.bound_id];

        const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
        std::vector<double> message(GN::n_ddbath_terms + GN::n_dimensions * ngp);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dbath, message);
        for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
            derivative.ddbath_at_baryctr_neigh[dbound.bound_id][ddbath] = message[ddbath];
        }
        for (uint gp = 0; gp < ngp; ++gp) {
            const uint gp_ex = ngp - gp - 1;
            for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                boundary.dbath_hat_at_gp(dbath, gp) = (boundary.dbath_hat_at_gp(dbath, gp) +
                                                       message[GN::n_ddbath_terms + GN::n_dimensions * gp_ex + dbath]) /
                                                      2.0;
            }
        }
    });

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& derivative = elt.data.derivative;
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            column(derivative.ddbath_at_midpts, bound) =
                derivative.ddbath_at_baryctr_neigh[bound] - derivative.ddbath_at_baryctr;
        }
        derivative.dddbath_at_baryctr = flatten<double, GN::n_dimensions, GN::n_ddbath_terms, SO::ColumnMajor>(
            derivative.P * transpose(derivative.ddbath_at_midpts));
    });
}
}
}

#endif
