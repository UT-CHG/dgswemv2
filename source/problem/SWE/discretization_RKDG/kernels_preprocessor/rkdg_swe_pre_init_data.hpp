#ifndef RKDG_SWE_PRE_INIT_DATA_HPP
#define RKDG_SWE_PRE_INIT_DATA_HPP

#include "utilities/file_exists.hpp"

namespace SWE {
namespace RKDG {
void Problem::initialize_data_serial(ProblemMeshType& mesh, const ProblemInputType& problem_specific_input) {
    SWE::initialize_data(mesh, problem_specific_input);

    // WETTING-DRYING INITIALIZE
    mesh.CallForEachElement([](auto& elt) {
        auto& shape = elt.GetShape();

        auto& state    = elt.data.state[0];
        auto& wd_state = elt.data.wet_dry_state;

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
            wd_state.bath_at_vrtx[vrtx] = shape.nodal_coordinates[vrtx][GlobalCoord::z];
        }

        wd_state.bath_min = *std::min_element(wd_state.bath_at_vrtx.begin(), wd_state.bath_at_vrtx.end());

        wd_state.q_lin = elt.ProjectBasisToLinear(state.q);

        wd_state.q_at_vrtx = elt.ComputeLinearUvrtx(wd_state.q_lin);

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
            wd_state.h_at_vrtx[vrtx] = wd_state.q_at_vrtx(SWE::Variables::ze, vrtx) + wd_state.bath_at_vrtx[vrtx];
        }

        bool set_wet_element = true;

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
            if (wd_state.h_at_vrtx[vrtx] <= PostProcessing::h_o + PostProcessing::h_o_threshold) {
                wd_state.q_at_vrtx(SWE::Variables::ze, vrtx) = PostProcessing::h_o - wd_state.bath_at_vrtx[vrtx];

                set_wet_element = false;
            }
        }

        if (set_wet_element) {
            wd_state.wet = true;
        } else {
            wd_state.wet = false;

            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
                wd_state.q_at_vrtx(SWE::Variables::qx, vrtx) = 0.0;
                wd_state.q_at_vrtx(SWE::Variables::qy, vrtx) = 0.0;
            }

            state.q = elt.ProjectLinearToBasis(elt.data.get_ndof(), wd_state.q_at_vrtx);

            set_constant(state.rhs, 0.0);
        }
    });

    // SLOPE LIMIT INITIALIZE
    mesh.CallForEachElement([](auto& elt) {
        auto& shape = elt.GetShape();

        auto& sl_state = elt.data.slope_limit_state;

        DynRowVector<double> bath_lin(elt.data.get_nvrtx());

        for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
            bath_lin[vrtx] = shape.nodal_coordinates[vrtx][GlobalCoord::z];
        }

        sl_state.bath_at_baryctr = elt.ComputeLinearUbaryctr(bath_lin);
        sl_state.bath_at_midpts  = elt.ComputeLinearUmidpts(bath_lin);

        sl_state.baryctr_coord = elt.GetShape().GetBarycentricCoordinates();
        sl_state.midpts_coord  = elt.GetShape().GetMidpointCoordinates();

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            sl_state.surface_normal[bound] = elt.GetShape().GetSurfaceNormal(bound, std::vector<Point<2>>(0))[0];
        }
    });

    mesh.CallForEachInterface([](auto& intface) {
        auto& sl_state_in = intface.data_in.slope_limit_state;
        auto& sl_state_ex = intface.data_ex.slope_limit_state;

        sl_state_in.baryctr_coord_neigh[intface.bound_id_in] = sl_state_ex.baryctr_coord;
        sl_state_ex.baryctr_coord_neigh[intface.bound_id_ex] = sl_state_in.baryctr_coord;
    });

    mesh.CallForEachBoundary([](auto& bound) {
        auto& sl_state = bound.data.slope_limit_state;

        sl_state.baryctr_coord_neigh[bound.bound_id][GlobalCoord::x] =
            2.0 * sl_state.midpts_coord[bound.bound_id][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
        sl_state.baryctr_coord_neigh[bound.bound_id][GlobalCoord::y] =
            2.0 * sl_state.midpts_coord[bound.bound_id][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];
    });

    mesh.CallForEachElement([](auto& elt) {
        auto& sl_state = elt.data.slope_limit_state;

        Array2D<double> A = Array2D<double>(2, std::vector<double>(2));
        std::vector<double> b(2);

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            uint element_1 = bound;
            uint element_2 = (bound + 1) % elt.data.get_nbound();

            if (!is_distributed(elt.GetBoundaryType()[element_1]) &&
                !is_distributed(elt.GetBoundaryType()[element_2])) {
                A[0][0] =
                    sl_state.baryctr_coord_neigh[element_1][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                A[1][0] =
                    sl_state.baryctr_coord_neigh[element_1][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];
                A[0][1] =
                    sl_state.baryctr_coord_neigh[element_2][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                A[1][1] =
                    sl_state.baryctr_coord_neigh[element_2][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

                b[0] = sl_state.midpts_coord[bound][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                b[1] = sl_state.midpts_coord[bound][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

                double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];

                sl_state.alpha_1[bound] = (A[1][1] * b[0] - A[0][1] * b[1]) / det;
                sl_state.alpha_2[bound] = (-A[1][0] * b[0] + A[0][0] * b[1]) / det;

                sl_state.r_sq[bound] = std::pow(b[0], 2.0) + std::pow(b[1], 2.0);
            }
        }
    });
}

void Problem::initialize_data_parallel_pre_send(ProblemMeshType& mesh, const ProblemInputType& problem_specific_input) {
    Problem::initialize_data_serial(mesh, problem_specific_input);

    mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& sl_state = dbound.data.slope_limit_state;

        // Construct message to exterior state
        std::vector<double> message;

        message.reserve(SWE::n_dimensions);

        for (uint dim = 0; dim < SWE::n_dimensions; ++dim) {
            message.push_back(sl_state.baryctr_coord[dim]);
        }

        // Set message to send buffer
        dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::baryctr_coord, message);
    });
}

void Problem::initialize_data_parallel_post_receive(ProblemMeshType& mesh) {
    mesh.CallForEachDistributedBoundary([](auto& dbound) {
        auto& sl_state = dbound.data.slope_limit_state;

        std::vector<double> message;

        message.resize(SWE::n_dimensions);

        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::baryctr_coord, message);

        for (uint dim = 0; dim < SWE::n_dimensions; ++dim) {
            sl_state.baryctr_coord_neigh[dbound.bound_id][dim] = message[dim];
        }
    });

    mesh.CallForEachElement([](auto& elt) {
        auto& sl_state = elt.data.slope_limit_state;

        Array2D<double> A = Array2D<double>(2, std::vector<double>(2));
        std::vector<double> b(2);

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            uint element_1 = bound;
            uint element_2 = (bound + 1) % elt.data.get_nbound();

            if (is_distributed(elt.GetBoundaryType()[element_1]) || is_distributed(elt.GetBoundaryType()[element_2])) {
                A[0][0] =
                    sl_state.baryctr_coord_neigh[element_1][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                A[1][0] =
                    sl_state.baryctr_coord_neigh[element_1][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];
                A[0][1] =
                    sl_state.baryctr_coord_neigh[element_2][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                A[1][1] =
                    sl_state.baryctr_coord_neigh[element_2][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

                b[0] = sl_state.midpts_coord[bound][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                b[1] = sl_state.midpts_coord[bound][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

                double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];

                sl_state.alpha_1[bound] = (A[1][1] * b[0] - A[0][1] * b[1]) / det;
                sl_state.alpha_2[bound] = (-A[1][0] * b[0] + A[0][0] * b[1]) / det;

                sl_state.r_sq[bound] = std::pow(b[0], 2.0) + std::pow(b[1], 2.0);
            }
        }
    });
}
}
}

#endif