#ifndef GREENGAUSS_DBATH_SERIAL_HPP
#define GREENGAUSS_DBATH_SERIAL_HPP

namespace GN {
namespace EHDG {
template <typename ProblemDiscretizationType>
void compute_dbath_gg(ProblemDiscretizationType& discretization);
template <typename ProblemDiscretizationType>
void compute_ddbath_gg(ProblemDiscretizationType& discretization);
template <typename ProblemDiscretizationType>
void compute_dddbath_gg(ProblemDiscretizationType& discretization);

void Problem::compute_bathymetry_derivatives_serial(ProblemDiscretizationType& discretization, ProblemGlobalDataType& global_data) {
    compute_dbath_gg(discretization);
    reconstruct_dbath(discretization, global_data);

    compute_ddbath_gg(discretization);
    reconstruct_ddbath(discretization, global_data);

    compute_dddbath_gg(discretization);
    reconstruct_dddbath(discretization, global_data);
}

template <typename ProblemDiscretizationType>
void compute_dbath_gg(ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& derivative = elt.data.derivative;
        set_constant(derivative.dbath_at_baryctr, 0);
    });

    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;
        auto& boundary_in   = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex   = intface.data_ex.boundary[intface.bound_id_ex];

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(derivative_in.dbath_at_baryctr, dir) +=
                1.0 / derivative_in.area *
                intface.IntegrationIN(vec_cw_mult(row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath),
                                                  row(intface.surface_normal_in, dir)));
            row(derivative_ex.dbath_at_baryctr, dir) +=
                1.0 / derivative_ex.area *
                intface.IntegrationEX(vec_cw_mult(row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath),
                                                  row(intface.surface_normal_ex, dir)));
        }
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& derivative = bound.data.derivative;
        auto& boundary   = bound.data.boundary[bound.bound_id];

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(derivative.dbath_at_baryctr, dir) +=
                1.0 / derivative.area *
                bound.Integration(
                    vec_cw_mult(row(boundary.aux_at_gp, SWE::Auxiliaries::bath), row(bound.surface_normal, dir)));
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& derivative = dbound.data.derivative;
        auto& boundary   = dbound.data.boundary[dbound.bound_id];

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(derivative.dbath_at_baryctr, dir) +=
                1.0 / derivative.area *
                dbound.Integration(
                    vec_cw_mult(row(boundary.aux_at_gp, SWE::Auxiliaries::bath), row(dbound.surface_normal, dir)));
        }
    });
}

template <typename ProblemDiscretizationType>
void compute_ddbath_gg(ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& derivative = elt.data.derivative;
        set_constant(derivative.ddbath_at_baryctr, 0);
    });

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

        for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative_in.ddbath_at_baryctr, GN::n_dimensions * dbath + dir) +=
                    1.0 / derivative_in.area *
                    intface.IntegrationIN(
                        vec_cw_mult(row(boundary_in.dbath_hat_at_gp, dbath), row(intface.surface_normal_in, dir)));
                row(derivative_ex.ddbath_at_baryctr, GN::n_dimensions * dbath + dir) +=
                    1.0 / derivative_ex.area *
                    intface.IntegrationEX(
                        vec_cw_mult(row(boundary_ex.dbath_hat_at_gp, dbath), row(intface.surface_normal_ex, dir)));
            }
        }
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& state      = bound.data.state[0];
        auto& derivative = bound.data.derivative;
        auto& boundary   = bound.data.boundary[bound.bound_id];

        boundary.dbath_hat_at_gp = bound.ComputeUgp(state.dbath);

        for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative.ddbath_at_baryctr, GN::n_dimensions * dbath + dir) +=
                    1.0 / derivative.area *
                    bound.Integration(
                        vec_cw_mult(row(boundary.dbath_hat_at_gp, dbath), row(bound.surface_normal, dir)));
            }
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& state      = dbound.data.state[0];
        auto& derivative = dbound.data.derivative;
        auto& boundary   = dbound.data.boundary[dbound.bound_id];

        boundary.dbath_hat_at_gp = dbound.ComputeUgp(state.dbath);

        for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative.ddbath_at_baryctr, GN::n_dimensions * dbath + dir) +=
                    1.0 / derivative.area *
                    dbound.Integration(
                        vec_cw_mult(row(boundary.dbath_hat_at_gp, dbath), row(dbound.surface_normal, dir)));
            }
        }
    });
}

template <typename ProblemDiscretizationType>
void compute_dddbath_gg(ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& derivative = elt.data.derivative;
        set_constant(derivative.dddbath_at_baryctr, 0);
    });

    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& state_in      = intface.data_in.state[0];
        auto& state_ex      = intface.data_ex.state[0];
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;

        const auto ddbath_at_gp_in = intface.ComputeUgpIN(state_in.ddbath);
        const auto ddbath_at_gp_ex = intface.ComputeUgpEX(state_ex.ddbath);

        for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative_in.dddbath_at_baryctr, GN::n_dimensions * ddbath + dir) +=
                    1.0 / derivative_in.area *
                    intface.IntegrationIN(
                        vec_cw_mult(row(ddbath_at_gp_in, ddbath), row(intface.surface_normal_in, dir)));
                row(derivative_ex.dddbath_at_baryctr, GN::n_dimensions * ddbath + dir) +=
                    1.0 / derivative_ex.area *
                    intface.IntegrationEX(
                        vec_cw_mult(row(ddbath_at_gp_ex, ddbath), row(intface.surface_normal_ex, dir)));
            }
        }
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& state      = bound.data.state[0];
        auto& derivative = bound.data.derivative;

        const auto ddbath_at_gp = bound.ComputeUgp(state.ddbath);

        for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative.dddbath_at_baryctr, GN::n_dimensions * ddbath + dir) +=
                    1.0 / derivative.area *
                    bound.Integration(vec_cw_mult(row(ddbath_at_gp, ddbath), row(bound.surface_normal, dir)));
            }
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& state      = dbound.data.state[0];
        auto& derivative = dbound.data.derivative;
        auto& boundary   = dbound.data.boundary[dbound.bound_id];

        const auto ddbath_at_gp = dbound.ComputeUgp(state.ddbath);

        for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative.dddbath_at_baryctr, GN::n_dimensions * ddbath + dir) +=
                    1.0 / derivative.area *
                    dbound.Integration(vec_cw_mult(row(ddbath_at_gp, ddbath), row(dbound.surface_normal, dir)));
            }
        }

        const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
        std::vector<double> message(GN::n_ddbath_terms + GN::n_dimensions * ngp);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dbath, message);
        for (uint gp = 0; gp < ngp; ++gp) {
            const uint gp_ex = ngp - gp - 1;
            for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                boundary.dbath_hat_at_gp(dbath, gp) = (boundary.dbath_hat_at_gp(dbath, gp) +
                                                       message[GN::n_ddbath_terms + GN::n_dimensions * gp_ex + dbath]) /
                                                      2.0;
            }
        }
    });
}
}
}

#endif
