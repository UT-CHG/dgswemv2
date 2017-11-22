#ifndef SWE_PRE_INIT_SL_DATA_HPP
#define SWE_PRE_INIT_SL_DATA_HPP

namespace SWE {
void Problem::initialize_sl_data_kernel(ProblemMeshType& mesh) {
    mesh.CallForEachElement([](auto& elt) {
        auto& state = elt.data.state[0];
        auto& sl_state = elt.data.slope_limit_state;

        elt.ComputeUbaryctr(state.bath, sl_state.bath_at_baryctr);
        elt.ComputeUmidpts(state.bath, sl_state.bath_at_midpts);

        sl_state.baryctr_coord = elt.GetShape().GetBarycentricCoordinates();
        sl_state.midpts_coord = elt.GetShape().GetMidpointCoordinates();

        for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
            sl_state.surface_normal[bound] = elt.GetShape().GetSurfaceNormal(bound, std::vector<Point<2>>(0))[0];
        }
    });

    mesh.CallForEachInterface([](auto& intface) {
        auto& sl_state_in = intface.data_in.slope_limit_state;
        auto& sl_state_ex = intface.data_ex.slope_limit_state;

        sl_state_in.bath_at_baryctr_neigh[intface.bound_id_in] = sl_state_ex.bath_at_baryctr;
        sl_state_ex.bath_at_baryctr_neigh[intface.bound_id_ex] = sl_state_in.bath_at_baryctr;

        sl_state_in.baryctr_coord_neigh[intface.bound_id_in] = sl_state_ex.baryctr_coord;
        sl_state_ex.baryctr_coord_neigh[intface.bound_id_ex] = sl_state_in.baryctr_coord;
    });

    mesh.CallForEachBoundary([](auto& bound) {
        auto& sl_state = bound.data.slope_limit_state;

        sl_state.bath_at_baryctr_neigh[bound.bound_id] = sl_state.bath_at_baryctr;

        sl_state.baryctr_coord_neigh[bound.bound_id][GlobalCoord::x] =
            2.0 * sl_state.midpts_coord[bound.bound_id][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
        sl_state.baryctr_coord_neigh[bound.bound_id][GlobalCoord::y] =
            2.0 * sl_state.midpts_coord[bound.bound_id][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];
    });

    mesh.CallForEachElement([](auto& elt) {
        auto& sl_state = elt.data.slope_limit_state;

        Array2D<double> A = Array2D<double>(2, std::vector<double>(2));
        std::vector<double> b(2);

        for (uint bound = 0; bound < elt.data.get_nbound(); bound++) {
            uint element_1 = bound;
            uint element_2 = (bound + 1) % 3;

            A[0][0] = sl_state.baryctr_coord_neigh[element_1][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
            A[1][0] = sl_state.baryctr_coord_neigh[element_1][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];
            A[0][1] = sl_state.baryctr_coord_neigh[element_2][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
            A[1][1] = sl_state.baryctr_coord_neigh[element_2][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

            b[0] = sl_state.midpts_coord[bound][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
            b[1] = sl_state.midpts_coord[bound][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

            double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];

            sl_state.alpha_1[bound] = (A[1][1] * b[0] - A[0][1] * b[1]) / det;
            sl_state.alpha_2[bound] = (-A[1][0] * b[0] + A[0][0] * b[1]) / det;

            sl_state.r_sq[bound] = std::pow(b[0], 2.0) + std::pow(b[1], 2.0);
        }
    });
}
}

#endif