#ifndef PROJECTION_DBATH_SERIAL_HPP
#define PROJECTION_DBATH_SERIAL_HPP

namespace GN {
namespace EHDG {
template <typename ProblemDiscretizationType>
void compute_dbath_rhs(ProblemDiscretizationType& discretization, const uint stage);
template <typename ProblemDiscretizationType>
void compute_ddbath_rhs(ProblemDiscretizationType& discretization, const uint stage);
template <typename ProblemDiscretizationType>
void compute_dddbath_rhs(ProblemDiscretizationType& discretization, const uint stage);

void Problem::compute_bathymetry_derivatives_serial(ProblemDiscretizationType& discretization,
                                                    ProblemGlobalDataType& global_data,
                                                    const uint stage) {
    compute_dbath_rhs(discretization, stage);
    discretization.mesh.CallForEachElement([stage](auto& elt) {
        auto& state = elt.data.state[stage];
        state.dbath = elt.ApplyMinv(elt.data.state[stage].dbath);
    });

    compute_ddbath_rhs(discretization, stage);
    discretization.mesh.CallForEachElement([stage](auto& elt) {
        auto& state  = elt.data.state[stage];
        state.ddbath = elt.ApplyMinv(elt.data.state[stage].ddbath);
    });

    compute_dddbath_rhs(discretization, stage);
    discretization.mesh.CallForEachElement([stage](auto& elt) {
        auto& state                     = elt.data.state[stage];
        state.dddbath                   = elt.ApplyMinv(state.dddbath);
        elt.data.internal.dddbath_at_gp = elt.ComputeUgp(state.dddbath);
    });
}

template <typename ProblemDiscretizationType>
void compute_dbath_rhs(ProblemDiscretizationType& discretization, const uint stage) {
    discretization.mesh.CallForEachElement([stage](auto& elt) {
        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(state.dbath, dir) = -elt.IntegrationDPhi(dir, row(internal.aux_at_gp, SWE::Auxiliaries::bath));
        }
    });

    discretization.mesh.CallForEachInterface([stage](auto& intface) {
        auto& state_in      = intface.data_in.state[stage];
        auto& state_ex      = intface.data_ex.state[stage];
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;
        auto& boundary_in   = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex   = intface.data_ex.boundary[intface.bound_id_ex];

        derivative_in.bath_hat_at_gp[intface.bound_id_in] = row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath);
        derivative_ex.bath_hat_at_gp[intface.bound_id_ex] = row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath);

        const uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        for (uint gp = 0; gp < ngp; ++gp) {
            const uint gp_ex = ngp - gp - 1;
            derivative_in.bath_hat_at_gp[intface.bound_id_in][gp] =
                (derivative_in.bath_hat_at_gp[intface.bound_id_in][gp] +
                 derivative_ex.bath_hat_at_gp[intface.bound_id_ex][gp_ex]) /
                2.0;
            derivative_ex.bath_hat_at_gp[intface.bound_id_ex][gp_ex] =
                derivative_in.bath_hat_at_gp[intface.bound_id_in][gp];
        }

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(state_in.dbath, dir) += intface.IntegrationPhiIN(
                vec_cw_mult(derivative_in.bath_hat_at_gp[intface.bound_id_in], row(intface.surface_normal_in, dir)));
            row(state_ex.dbath, dir) += intface.IntegrationPhiEX(
                vec_cw_mult(derivative_ex.bath_hat_at_gp[intface.bound_id_ex], row(intface.surface_normal_ex, dir)));
        }
    });

    discretization.mesh.CallForEachBoundary([stage](auto& bound) {
        auto& state      = bound.data.state[stage];
        auto& derivative = bound.data.derivative;
        auto& boundary   = bound.data.boundary[bound.bound_id];

        derivative.bath_hat_at_gp[bound.bound_id] = row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(state.dbath, dir) += bound.IntegrationPhi(
                vec_cw_mult(derivative.bath_hat_at_gp[bound.bound_id], row(bound.surface_normal, dir)));
        }
    });
}

template <typename ProblemDiscretizationType>
void compute_ddbath_rhs(ProblemDiscretizationType& discretization, const uint stage) {
    discretization.mesh.CallForEachElement([stage](auto& elt) {
        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;

        internal.dbath_at_gp = elt.ComputeUgp(state.dbath);

        for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.ddbath, GN::n_dimensions * dbath + dir) =
                    -elt.IntegrationDPhi(dir, row(internal.dbath_at_gp, dbath));
            }
        }
    });

    discretization.mesh.CallForEachInterface([stage](auto& intface) {
        auto& state_in    = intface.data_in.state[stage];
        auto& state_ex    = intface.data_ex.state[stage];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

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
                row(state_in.ddbath, GN::n_dimensions * dbath + dir) += intface.IntegrationPhiIN(
                    vec_cw_mult(row(boundary_in.dbath_hat_at_gp, dbath), row(intface.surface_normal_in, dir)));
                row(state_ex.ddbath, GN::n_dimensions * dbath + dir) += intface.IntegrationPhiEX(
                    vec_cw_mult(row(boundary_ex.dbath_hat_at_gp, dbath), row(intface.surface_normal_ex, dir)));
            }
        }
    });

    discretization.mesh.CallForEachBoundary([stage](auto& bound) {
        auto& state    = bound.data.state[stage];
        auto& boundary = bound.data.boundary[bound.bound_id];

        boundary.dbath_hat_at_gp = bound.ComputeUgp(state.dbath);

        for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.ddbath, GN::n_dimensions * dbath + dir) += bound.IntegrationPhi(
                    vec_cw_mult(row(boundary.dbath_hat_at_gp, dbath), row(bound.surface_normal, dir)));
            }
        }
    });
}

template <typename ProblemDiscretizationType>
void compute_dddbath_rhs(ProblemDiscretizationType& discretization, const uint stage) {
    discretization.mesh.CallForEachElement([stage](auto& elt) {
        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;

        internal.ddbath_at_gp = elt.ComputeUgp(state.ddbath);

        for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.dddbath, GN::n_dimensions * ddbath + dir) =
                    -elt.IntegrationDPhi(dir, row(internal.ddbath_at_gp, ddbath));
            }
        }
    });

    discretization.mesh.CallForEachInterface([stage](auto& intface) {
        auto& state_in      = intface.data_in.state[stage];
        auto& state_ex      = intface.data_ex.state[stage];
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;

        derivative_in.ddbath_hat_at_gp[intface.bound_id_in] = intface.ComputeUgpIN(state_in.ddbath);
        derivative_ex.ddbath_hat_at_gp[intface.bound_id_ex] = intface.ComputeUgpEX(state_ex.ddbath);

        const uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        for (uint gp = 0; gp < ngp; ++gp) {
            const uint gp_ex = ngp - gp - 1;
            for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
                derivative_in.ddbath_hat_at_gp[intface.bound_id_in](ddbath, gp) =
                    (derivative_in.ddbath_hat_at_gp[intface.bound_id_in](ddbath, gp) +
                     derivative_ex.ddbath_hat_at_gp[intface.bound_id_ex](ddbath, gp_ex)) /
                    2.0;
                derivative_ex.ddbath_hat_at_gp[intface.bound_id_ex](ddbath, gp_ex) =
                    derivative_in.ddbath_hat_at_gp[intface.bound_id_in](ddbath, gp);
            }
        }

        for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state_in.dddbath, GN::n_dimensions * ddbath + dir) += intface.IntegrationPhiIN(
                    vec_cw_mult(row(derivative_in.ddbath_hat_at_gp[intface.bound_id_in], ddbath),
                                row(intface.surface_normal_in, dir)));
                row(state_ex.dddbath, GN::n_dimensions * ddbath + dir) += intface.IntegrationPhiEX(
                    vec_cw_mult(row(derivative_ex.ddbath_hat_at_gp[intface.bound_id_ex], ddbath),
                                row(intface.surface_normal_ex, dir)));
            }
        }
    });

    discretization.mesh.CallForEachBoundary([stage](auto& bound) {
        auto& state      = bound.data.state[stage];
        auto& derivative = bound.data.derivative;

        derivative.ddbath_hat_at_gp[bound.bound_id] = bound.ComputeUgp(state.ddbath);

        for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.dddbath, GN::n_dimensions * ddbath + dir) += bound.IntegrationPhi(vec_cw_mult(
                    row(derivative.ddbath_hat_at_gp[bound.bound_id], ddbath), row(bound.surface_normal, dir)));
            }
        }
    });
}
}
}

#endif