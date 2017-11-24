#ifndef SWE_PRE_INIT_WD_DATA_HPP
#define SWE_PRE_INIT_WD_DATA_HPP

namespace SWE {
void Problem::initialize_wd_data_kernel(ProblemMeshType& mesh) {
    mesh.CallForEachElement([](auto& elt) {
        auto& state = elt.data.state[0];
        auto& internal = elt.data.internal;
        auto& wd_state = elt.data.wet_dry_state;

        elt.ComputeUvrtx(state.bath, wd_state.bath_at_vrtx);
        wd_state.bath_min = *std::min_element(wd_state.bath_at_vrtx.begin(), wd_state.bath_at_vrtx.end());

        elt.ComputeUvrtx(state.ze, wd_state.ze_at_vrtx);

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            wd_state.h_at_vrtx[vrtx] = wd_state.ze_at_vrtx[vrtx] + wd_state.bath_at_vrtx[vrtx];
        }

        bool set_wet_element = true;

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); vrtx++) {
            if (wd_state.h_at_vrtx[vrtx] <= Global::h_o) {
                wd_state.ze_at_vrtx[vrtx] = Global::h_o - wd_state.bath_at_vrtx[vrtx];

                set_wet_element = false;
            }
        }

        if (set_wet_element) {
            wd_state.wet = true;
        } else {
            wd_state.wet = false;

            state.ze = elt.L2Projection(wd_state.ze_at_vrtx);
            std::fill(state.qx.begin(), state.qx.end(), 0.0);
            std::fill(state.qy.begin(), state.qy.end(), 0.0);

            std::fill(state.rhs_ze.begin(), state.rhs_ze.end(), 0.0);
            std::fill(state.rhs_qx.begin(), state.rhs_qx.end(), 0.0);
            std::fill(state.rhs_qy.begin(), state.rhs_qy.end(), 0.0);
        }

        elt.ComputeUgp(state.ze, internal.ze_at_gp);

        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            internal.h_at_gp[gp] = internal.ze_at_gp[gp] + internal.bath_at_gp[gp];
        }

        wd_state.water_volume = elt.Integration(internal.h_at_gp);
    });
}
}

#endif