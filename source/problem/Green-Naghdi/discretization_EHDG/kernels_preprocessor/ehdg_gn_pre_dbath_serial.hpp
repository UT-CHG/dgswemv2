#ifndef EHDG_GN_PRE_DBATH_SERIAL_HPP
#define EHDG_GN_PRE_DBATH_SERIAL_HPP

namespace GN {
namespace EHDG {
template <typename ProblemDiscretizationType>
void compute_dbath_avg(ProblemDiscretizationType& discretization);
template <typename ProblemDiscretizationType>
void compute_ddbath_avg(ProblemDiscretizationType& discretization);
template <typename ProblemDiscretizationType>
void compute_dddbath_avg(ProblemDiscretizationType& discretization);

template <typename ProblemDiscretizationType>
void compute_dbath(ProblemDiscretizationType& discretization);
template <typename ProblemDiscretizationType>
void compute_ddbath(ProblemDiscretizationType& discretization);
template <typename ProblemDiscretizationType>
void compute_dddbath(ProblemDiscretizationType& discretization);

void Problem::compute_bathymetry_derivatives_serial(ProblemDiscretizationType& discretization) {
    compute_dbath_avg(discretization);
    compute_dbath(discretization);

    compute_ddbath_avg(discretization);
    compute_ddbath(discretization);

    compute_dddbath_avg(discretization);
    compute_dddbath(discretization);
}

template <typename ProblemDiscretizationType>
void compute_dbath_avg(ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& derivative = elt.data.derivative;
        set_constant(derivative.dbath_at_baryctr, 0);
    });

    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& derivative_in    = intface.data_in.derivative;
        auto& derivative_ex    = intface.data_ex.derivative;
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(derivative_in.dbath_at_baryctr, dir) +=
                    1.0 / derivative_in.area * intface.IntegrationIN(
                            vec_cw_mult(
                                    row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath),
                                    row(intface.surface_normal_in, dir)));
            row(derivative_ex.dbath_at_baryctr, dir) +=
                    1.0 / derivative_ex.area * intface.IntegrationEX(
                            vec_cw_mult(
                                    row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath),
                                    row(intface.surface_normal_ex, dir)));
        }
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& derivative    = bound.data.derivative;
        auto& boundary = bound.data.boundary[bound.bound_id];

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(derivative.dbath_at_baryctr, dir) +=
                    1.0 / derivative.area * bound.Integration(
                            vec_cw_mult(
                                    row(boundary.aux_at_gp, SWE::Auxiliaries::bath),
                                    row(bound.surface_normal, dir)));
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& derivative    = dbound.data.derivative;
        auto& boundary = dbound.data.boundary[dbound.bound_id];

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(derivative.dbath_at_baryctr, dir) +=
                    1.0 / derivative.area * dbound.Integration(
                            vec_cw_mult(
                                    row(boundary.aux_at_gp, SWE::Auxiliaries::bath),
                                    row(dbound.surface_normal, dir)));
        }
    });
}

template <typename ProblemDiscretizationType>
void compute_dbath(ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;
        derivative_in.dbath_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.dbath_at_baryctr;
        derivative_ex.dbath_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.dbath_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& derivative = bound.data.derivative;
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
        auto& state = elt.data.state[0];
        auto& derivative = elt.data.derivative;
        auto& internal = elt.data.internal;

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            const uint element_1 = derivative.a_elem[2 * bound];
            const uint element_2 = derivative.a_elem[2 * bound + 1];
            column(derivative.dbath_at_midpts, bound) = derivative.dbath_at_baryctr +
                                                   (derivative.dbath_at_baryctr_neigh[element_1] - derivative.dbath_at_baryctr) * derivative.a[bound][0] +
                                                   (derivative.dbath_at_baryctr_neigh[element_2] - derivative.dbath_at_baryctr) * derivative.a[bound][1];
        }
        derivative.dbath_lin = derivative.dbath_at_midpts * derivative.T;
        state.dbath = elt.ProjectLinearToBasis(derivative.dbath_lin);
        internal.dbath_at_gp = elt.ComputeUgp(state.dbath);
    });
}

template <typename ProblemDiscretizationType>
void compute_ddbath_avg(ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& derivative = elt.data.derivative;
        set_constant(derivative.ddbath_at_baryctr, 0);
    });

    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& state_in    = intface.data_in.state[0];
        auto& state_ex    = intface.data_ex.state[0];
        auto& derivative_in    = intface.data_in.derivative;
        auto& derivative_ex    = intface.data_ex.derivative;
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        boundary_in.dbath_hat_at_gp = intface.ComputeUgpIN(state_in.dbath);
        boundary_ex.dbath_hat_at_gp = intface.ComputeUgpEX(state_ex.dbath);

        for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative_in.ddbath_at_baryctr, GN::n_dimensions * dbath + dir) +=
                        1.0 / derivative_in.area * intface.IntegrationIN(
                                vec_cw_mult(
                                        row(boundary_in.dbath_hat_at_gp, dbath),
                                        row(intface.surface_normal_in, dir)));
                row(derivative_ex.ddbath_at_baryctr, GN::n_dimensions * dbath + dir) +=
                        1.0 / derivative_ex.area * intface.IntegrationEX(
                                vec_cw_mult(
                                        row(boundary_ex.dbath_hat_at_gp, dbath),
                                        row(intface.surface_normal_ex, dir)));
            }
        }

        const uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        for (uint gp = 0; gp < ngp; ++gp) {
            const uint gp_ex = ngp - gp - 1;
            for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                boundary_in.dbath_hat_at_gp(dbath, gp) =
                        (boundary_in.dbath_hat_at_gp(dbath, gp) + boundary_ex.dbath_hat_at_gp(dbath, gp_ex)) / 2.0;
                boundary_ex.dbath_hat_at_gp(dbath, gp_ex) = boundary_in.dbath_hat_at_gp(dbath, gp);
            }
        }
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& state    = bound.data.state[0];
        auto& derivative    = bound.data.derivative;
        auto& boundary = bound.data.boundary[bound.bound_id];

        boundary.dbath_hat_at_gp = bound.ComputeUgp(state.dbath);

        for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative.ddbath_at_baryctr, GN::n_dimensions * dbath + dir) +=
                        1.0 / derivative.area * bound.Integration(
                                vec_cw_mult(
                                        row(boundary.dbath_hat_at_gp, dbath),
                                        row(bound.surface_normal, dir)));
            }
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& state    = dbound.data.state[0];
        auto& derivative    = dbound.data.derivative;
        auto& boundary = dbound.data.boundary[dbound.bound_id];

        boundary.dbath_hat_at_gp = dbound.ComputeUgp(state.dbath);

        for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative.ddbath_at_baryctr, GN::n_dimensions * dbath + dir) +=
                        1.0 / derivative.area * dbound.Integration(
                                vec_cw_mult(
                                        row(boundary.dbath_hat_at_gp, dbath),
                                        row(dbound.surface_normal, dir)));
            }
        }
    });
}

template <typename ProblemDiscretizationType>
void compute_ddbath(ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;
        derivative_in.ddbath_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.ddbath_at_baryctr;
        derivative_ex.ddbath_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.ddbath_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& derivative = bound.data.derivative;
        derivative.ddbath_at_baryctr_neigh[bound.bound_id] = derivative.ddbath_at_baryctr;
    });

    discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& derivative = dbound.data.derivative;
        auto& boundary = dbound.data.boundary[dbound.bound_id];

        const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
        std::vector<double> message(GN::n_ddbath_terms + GN::n_dimensions * ngp);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dbath, message);
        for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
            derivative.ddbath_at_baryctr_neigh[dbound.bound_id][ddbath] = message[ddbath];
        }
        for (uint gp = 0; gp < ngp; ++gp) {
            const uint gp_ex = ngp - gp - 1;
            for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                boundary.dbath_hat_at_gp(dbath, gp) =
                        (boundary.dbath_hat_at_gp(dbath, gp) + message[GN::n_ddbath_terms + GN::n_dimensions * gp_ex + dbath]) / 2.0;
            }
        }
    });

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& state = elt.data.state[0];
        auto& derivative = elt.data.derivative;
        auto& internal = elt.data.internal;

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            const uint element_1 = derivative.a_elem[2 * bound];
            const uint element_2 = derivative.a_elem[2 * bound + 1];
            column(derivative.ddbath_at_midpts, bound) = derivative.ddbath_at_baryctr +
                                                    (derivative.ddbath_at_baryctr_neigh[element_1] - derivative.ddbath_at_baryctr) * derivative.a[bound][0] +
                                                    (derivative.ddbath_at_baryctr_neigh[element_2] - derivative.ddbath_at_baryctr) * derivative.a[bound][1];
        }
        derivative.ddbath_lin = derivative.ddbath_at_midpts * derivative.T;
        state.ddbath = elt.ProjectLinearToBasis(derivative.ddbath_lin);
        internal.ddbath_at_gp = elt.ComputeUgp(state.ddbath);
    });
}

template <typename ProblemDiscretizationType>
void compute_dddbath_avg(ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& derivative = elt.data.derivative;
        set_constant(derivative.dddbath_at_baryctr, 0);
    });

    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& state_in    = intface.data_in.state[0];
        auto& state_ex    = intface.data_ex.state[0];
        auto& derivative_in    = intface.data_in.derivative;
        auto& derivative_ex    = intface.data_ex.derivative;

        auto ddbath_at_gp_in = intface.ComputeUgpIN(state_in.ddbath);
        auto ddbath_at_gp_ex = intface.ComputeUgpEX(state_ex.ddbath);

        for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative_in.dddbath_at_baryctr, GN::n_dimensions * ddbath + dir) +=
                        1.0 / derivative_in.area * intface.IntegrationIN(
                                vec_cw_mult(
                                        row(ddbath_at_gp_in, ddbath),
                                        row(intface.surface_normal_in, dir)));
                row(derivative_ex.dddbath_at_baryctr, GN::n_dimensions * ddbath + dir) +=
                        1.0 / derivative_ex.area * intface.IntegrationEX(
                                vec_cw_mult(
                                        row(ddbath_at_gp_ex, ddbath),
                                        row(intface.surface_normal_ex, dir)));
            }
        }
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& state    = bound.data.state[0];
        auto& derivative    = bound.data.derivative;

        auto ddbath_at_gp = bound.ComputeUgp(state.ddbath);

        for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative.dddbath_at_baryctr, GN::n_dimensions * ddbath + dir) +=
                        1.0 / derivative.area * bound.Integration(
                                vec_cw_mult(
                                        row(ddbath_at_gp, ddbath),
                                        row(bound.surface_normal, dir)));
            }
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& state    = dbound.data.state[0];
        auto& derivative    = dbound.data.derivative;

        auto ddbath_at_gp = dbound.ComputeUgp(state.ddbath);

        for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative.dddbath_at_baryctr, GN::n_dimensions * ddbath + dir) +=
                        1.0 / derivative.area * dbound.Integration(
                                vec_cw_mult(
                                        row(ddbath_at_gp, ddbath),
                                        row(dbound.surface_normal, dir)));
            }
        }
    });
}

template <typename ProblemDiscretizationType>
void compute_dddbath(ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;
        derivative_in.dddbath_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.dddbath_at_baryctr;
        derivative_ex.dddbath_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.dddbath_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& derivative = bound.data.derivative;
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
        auto& state = elt.data.state[0];
        auto& derivative = elt.data.derivative;
        auto& internal = elt.data.internal;

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            const uint element_1 = derivative.a_elem[2 * bound];
            const uint element_2 = derivative.a_elem[2 * bound + 1];
            column(derivative.dddbath_at_midpts, bound) = derivative.dddbath_at_baryctr +
                                                     (derivative.dddbath_at_baryctr_neigh[element_1] - derivative.dddbath_at_baryctr) * derivative.a[bound][0] +
                                                     (derivative.dddbath_at_baryctr_neigh[element_2] - derivative.dddbath_at_baryctr) * derivative.a[bound][1];
        }
        derivative.dddbath_lin = derivative.dddbath_at_midpts * derivative.T;
        state.dddbath = elt.ProjectLinearToBasis(derivative.dddbath_lin);
        internal.dddbath_at_gp = elt.ComputeUgp(state.dddbath);
     });
}
}
}

#endif
