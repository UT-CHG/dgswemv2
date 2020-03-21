#ifndef SWE_TROUBLE_CHECK_HPP
#define SWE_TROUBLE_CHECK_HpP

namespace SWE {
template <typename DiscretizationType, typename StepperType>
void check_trouble(DiscretizationType& discretization, const StepperType& stepper, uint comm_type) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        elt.data.slope_limit_state.I         = 0.0;
        elt.data.slope_limit_state.perimeter = 0.0;
        elt.data.slope_limit_state.troubled  = false;
    });

    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        auto& state_in    = intface.data_in.state[stepper.GetStage()];
        auto& state_ex    = intface.data_ex.state[stepper.GetStage()];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];
        auto& sl_state_in = intface.data_in.slope_limit_state;
        auto& sl_state_ex = intface.data_ex.slope_limit_state;

        if (intface.data_in.wet_dry_state.wet || intface.data_ex.wet_dry_state.wet) {
            boundary_in.q_at_gp = intface.ComputeUgpIN(state_in.q);
            boundary_ex.q_at_gp = intface.ComputeUgpEX(state_ex.q);
            row(boundary_in.aux_at_gp, SWE::Auxiliaries::h) =
                row(boundary_in.q_at_gp, SWE::Variables::ze) + row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath);
            row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h) =
                row(boundary_ex.q_at_gp, SWE::Variables::ze) + row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath);

            const auto qn_in =
                vec_cw_mult(row(intface.surface_normal_in, GlobalCoord::x), row(boundary_in.q_at_gp, GlobalCoord::x)) +
                vec_cw_mult(row(intface.surface_normal_in, GlobalCoord::y), row(boundary_in.q_at_gp, GlobalCoord::y));

            const auto qn_ex =
                vec_cw_mult(row(intface.surface_normal_ex, GlobalCoord::x), row(boundary_ex.q_at_gp, GlobalCoord::x)) +
                vec_cw_mult(row(intface.surface_normal_ex, GlobalCoord::y), row(boundary_ex.q_at_gp, GlobalCoord::y));

            if (intface.IntegrationIN(qn_in) < 0.0) {
                const uint ngp                              = intface.data_in.get_ngp_boundary(intface.bound_id_in);
                sl_state_in.hdif_at_gp[intface.bound_id_in] = row(boundary_in.aux_at_gp, SWE::Auxiliaries::h);
                for (uint gp = 0; gp < ngp; ++gp) {
                    const uint gp_ex = ngp - gp - 1;
                    sl_state_in.hdif_at_gp[intface.bound_id_in][gp] -=
                        boundary_ex.aux_at_gp(SWE::Auxiliaries::h, gp_ex);
                }
                sl_state_in.perimeter += sl_state_in.lengths[intface.bound_id_in];
                sl_state_in.I += std::abs(intface.IntegrationIN(sl_state_in.hdif_at_gp[intface.bound_id_in]));
            }

            if (intface.IntegrationEX(qn_ex) < 0.0) {
                const uint ngp                              = intface.data_ex.get_ngp_boundary(intface.bound_id_ex);
                sl_state_ex.hdif_at_gp[intface.bound_id_ex] = row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h);
                for (uint gp = 0; gp < ngp; ++gp) {
                    const uint gp_ex = ngp - gp - 1;
                    sl_state_ex.hdif_at_gp[intface.bound_id_ex][gp] -=
                        boundary_in.aux_at_gp(SWE::Auxiliaries::h, gp_ex);
                }
                sl_state_ex.perimeter += sl_state_ex.lengths[intface.bound_id_ex];
                sl_state_ex.I += std::abs(intface.IntegrationEX(sl_state_ex.hdif_at_gp[intface.bound_id_ex]));
            }
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([&stepper, comm_type](auto& dbound) {
        auto& boundary = dbound.data.boundary[dbound.bound_id];
        auto& sl_state = dbound.data.slope_limit_state;

        if (dbound.data.wet_dry_state.wet) {
            const auto qn =
                vec_cw_mult(row(dbound.surface_normal, GlobalCoord::x), row(boundary.q_at_gp, GlobalCoord::x)) +
                vec_cw_mult(row(dbound.surface_normal, GlobalCoord::y), row(boundary.q_at_gp, GlobalCoord::y));

            if (dbound.Integration(qn) < 0.0) {
                const uint ngp                       = dbound.data.get_ngp_boundary(dbound.bound_id);
                sl_state.hdif_at_gp[dbound.bound_id] = row(boundary.aux_at_gp, SWE::Auxiliaries::h);
                std::vector<double> message(1 + SWE::n_variables + ngp);
                dbound.boundary_condition.exchanger.GetFromReceiveBuffer(comm_type, message);
                for (uint gp = 0; gp < ngp; ++gp) {
                    const uint gp_ex = ngp - gp - 1;
                    sl_state.hdif_at_gp[dbound.bound_id][gp] -= message[1 + SWE::n_variables + gp_ex];
                }
                sl_state.perimeter += sl_state.lengths[dbound.bound_id];
                sl_state.I += std::abs(dbound.Integration(sl_state.hdif_at_gp[dbound.bound_id]));
            }
        }
    });

    std::set<uint> troubled_area;
    const uint p = discretization.mesh.GetP();
    discretization.mesh.CallForEachElement([&stepper, &troubled_area, p](auto& elt) {
        auto& sl_state = elt.data.slope_limit_state;
        if (elt.data.wet_dry_state.wet) {
            if (sl_state.perimeter != 0.0) {
                double h_norm = 0.0;
                for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
                    h_norm =
                        std::max(h_norm, sl_state.q_at_vrtx(SWE::Variables::ze, vrtx) + sl_state.bath_at_vrtx[vrtx]);
                }
                sl_state.I /= (std::pow(sl_state.radius, (p + 1.0) / 2.0) * sl_state.perimeter * h_norm);
                if (sl_state.I > 1.0) {
                    troubled_area.insert(elt.GetID());
                }
            }
        }
    });

    for (uint pass = 0; pass < 6; ++pass) {
        std::set<uint> troubled_area_new;
        discretization.mesh.CallForEachElement([&troubled_area, &troubled_area_new](auto& elt) {
            if (troubled_area.find(elt.GetID()) != troubled_area.end()) {
                elt.data.slope_limit_state.troubled = true;
                troubled_area_new.insert(elt.GetNeighborID().begin(), elt.GetNeighborID().end());
            }
        });
        troubled_area = troubled_area_new;
    }
}
}

#endif