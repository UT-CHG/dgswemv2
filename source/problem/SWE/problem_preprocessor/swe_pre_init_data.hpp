#ifndef SWE_PRE_INIT_DATA_HPP
#define SWE_PRE_INIT_DATA_HPP

#include "utilities/file_exists.hpp"
#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"

namespace SWE {
namespace detail {

template <typename IntfaceSoA>
struct AuxAssembler {
    template <typename AccessorType>
    constexpr static bool is_vectorized() { return true; }

    AuxAssembler(uint element_type_index_, IntfaceSoA& interface_soa_)
        : element_type_index(element_type_index_), interface_soa(interface_soa_)
    {}

    template <typename ElementSoA>
    void operator() (ElementSoA& soa) const {
        for ( uint var = 0; var < SWE::n_auxiliaries; ++var ) {

            set_constant(interface_soa.data.aux_at_gp[var], 0.0);

            for (uint side = 0; side < 3; ++side ) {
                const auto& aux = soa.data.boundary[side].aux_at_gp[var];
                interface_soa.data.aux_at_gp[var] += interface_soa.ScatterIn(element_type_index, side, aux );
            }
        }
    }

private:
    uint element_type_index;
    IntfaceSoA& interface_soa;
};

struct InitializeInterfaceSoA {
    template <typename AccessorType>
    constexpr static bool is_vectorized() { return true; }

    template <typename SoA>
    void operator() (SoA& soa) const {
        static_assert(Utilities::is_SoA<SoA>::value);

        if ( soa.Getninterfaces() == 0 ) return;

        uint element_type_index = 0u;
        Utilities::for_each_in_tuple(soa.GetElementData(), [&soa, &element_type_index](auto& elt_container) {
                AuxAssembler<SoA> aux_assembler(element_type_index++, soa);

                using is_vectorized_t = std::true_type;
                elt_container.CallForEachElement(aux_assembler, is_vectorized_t{} );
            });

    }

};
}

template <typename MeshType, typename ProblemSpecificInputType>
void initialize_data_serial(MeshType& mesh, const ProblemSpecificInputType& problem_specific_input) {
    mesh.CallForEachElement([&problem_specific_input](auto& elt) {
        elt.data.initialize();

        auto& shape = elt.GetShape();

        auto& state    = elt.data.state[0];
        auto& internal = elt.data.internal;

        uint nnode = elt.data.get_nnode();
        uint ngp   = elt.data.get_ngp_internal();

        DynRowVector<double> bathymetry(nnode);

        for (uint node_id = 0; node_id < nnode; ++node_id) {
            bathymetry[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::z];
        }

        state.aux = elt.L2ProjectionNode(bathymetry);

        internal.aux_at_gp[SWE::Auxiliaries::bath] = elt.ComputeNodalUgp(bathymetry);
        row(internal.dbath_at_gp, GlobalCoord::x)       = elt.ComputeNodalDUgp(GlobalCoord::x, bathymetry);
        row(internal.dbath_at_gp, GlobalCoord::y)       = elt.ComputeNodalDUgp(GlobalCoord::y, bathymetry);

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);

            for (uint node_id = 0; node_id < nnode; ++node_id) {
                y_node[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::y];
            }

            DynRowVector<double> y_at_gp(ngp);

            y_at_gp = elt.ComputeNodalUgp(y_node);

            double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            double R         = problem_specific_input.spherical_projection.R;

            for (uint gp = 0; gp < ngp; ++gp) {
                internal.aux_at_gp[SWE::Auxiliaries::sp][gp] = cos_phi_o / std::cos(y_at_gp[gp] / R);
            }
        } else {
            for (uint gp = 0; gp < ngp; ++gp) {
                internal.aux_at_gp[SWE::Auxiliaries::sp][gp] = 1.0;
            }
        }

        if (problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Constant ||
            problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Default) {
            uint nnode = elt.GetShape().nodal_coordinates.size();

            DynRowVector<double> u_node(nnode);

            set_constant(u_node, problem_specific_input.initial_conditions.ze_initial);
            state.q[SWE::Variables::ze] = elt.L2ProjectionNode(u_node);

            set_constant(u_node, problem_specific_input.initial_conditions.qx_initial);
            state.q[SWE::Variables::qx] = elt.L2ProjectionNode(u_node);

            set_constant(u_node, problem_specific_input.initial_conditions.qy_initial);
            state.q[SWE::Variables::qy] = elt.L2ProjectionNode(u_node);

        } else if (problem_specific_input.initial_conditions.type == SWE::InitialConditionsType::Function) {
            DynMatrix<double> tmp = elt.L2ProjectionF([](Point<2>& pt) { return SWE::ic_q(0, pt); });
            for ( uint var = 0; var < SWE::n_variables; ++var ) {
                state.q[var] = row(tmp,var);
            }
        }
    });

    mesh.CallForEachInterface([&problem_specific_input](auto& intface) {
        auto& shape_in = intface.GetShapeIN();
        auto& shape_ex = intface.GetShapeEX();

        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        uint nnode = intface.data_in.get_nnode();
        uint ngp   = intface.data_in.get_ngp_boundary(intface.bound_id_in);

        DynRowVector<double> bathymetry(nnode);

        for (uint node_id = 0; node_id < nnode; ++node_id) {
            bathymetry[node_id] = shape_in.nodal_coordinates[node_id][GlobalCoord::z];
        }

        boundary_in.aux_at_gp[SWE::Auxiliaries::bath] = intface.ComputeNodalUgpIN(bathymetry);

        // bathymetry is not necessarily continuous across weir nodes
        for (uint node_id = 0; node_id < nnode; ++node_id) {
            bathymetry[node_id] = shape_ex.nodal_coordinates[node_id][GlobalCoord::z];
        }

        boundary_ex.aux_at_gp[SWE::Auxiliaries::bath] = intface.ComputeNodalUgpEX(bathymetry);

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);

            for (uint node_id = 0; node_id < intface.data_in.get_nnode(); ++node_id) {
                y_node[node_id] = shape_in.nodal_coordinates[node_id][GlobalCoord::y];
            }

            DynRowVector<double> y_at_gp(ngp);

            y_at_gp = intface.ComputeNodalUgpIN(y_node);

            double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            double R         = problem_specific_input.spherical_projection.R;

            for ( uint gp = 0; gp < ngp; ++gp ) {
                boundary_in.aux_at_gp[SWE::Auxiliaries::sp][gp] = cos_phi_o / cos_vec(y_at_gp[gp] / R);
                boundary_ex.aux_at_gp[SWE::Auxiliaries::sp][gp] = cos_phi_o / cos_vec(y_at_gp[gp] / R);
            }

        } else {

            set_constant(boundary_in.aux_at_gp[SWE::Auxiliaries::sp], 1.0);
            set_constant(boundary_ex.aux_at_gp[SWE::Auxiliaries::sp], 1.0);

        }
    });

    mesh.CallForEachInterface( detail::InitializeInterfaceSoA{} );

    mesh.CallForEachBoundary([&problem_specific_input](auto& bound) {
        auto& shape = bound.GetShape();

        auto& boundary = bound.data.boundary[bound.bound_id];

        uint nnode = bound.data.get_nnode();
        uint ngp   = bound.data.get_ngp_boundary(bound.bound_id);

        DynRowVector<double> bathymetry(nnode);

        for (uint node_id = 0; node_id < nnode; ++node_id) {
            bathymetry[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::z];
        }

        boundary.aux_at_gp[SWE::Auxiliaries::bath] = bound.ComputeNodalUgp(bathymetry);

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);

            for (uint node_id = 0; node_id < bound.data.get_nnode(); ++node_id) {
                y_node[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::y];
            }

            DynRowVector<double> y_at_gp(ngp);

            y_at_gp = bound.ComputeNodalUgp(y_node);

            double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            double R         = problem_specific_input.spherical_projection.R;

            for (uint gp = 0; gp < ngp; ++gp) {
                boundary.aux_at_gp[SWE::Auxiliaries::sp][gp] = cos_phi_o / std::cos(y_at_gp[gp] / R);
            }
        } else {
            for (uint gp = 0; gp < ngp; ++gp) {
                boundary.aux_at_gp[SWE::Auxiliaries::sp][gp] = 1.0;
            }
        }
    });

    mesh.CallForEachDistributedBoundary([&problem_specific_input](auto& dbound) {
        auto& shape = dbound.GetShape();

        auto& boundary = dbound.data.boundary[dbound.bound_id];

        uint nnode = dbound.data.get_nnode();
        uint ngp   = dbound.data.get_ngp_boundary(dbound.bound_id);

        DynRowVector<double> bathymetry(nnode);

        for (uint node_id = 0; node_id < nnode; ++node_id) {
            bathymetry[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::z];
        }

        boundary.aux_at_gp[SWE::Auxiliaries::bath] = dbound.ComputeNodalUgp(bathymetry);

        if (problem_specific_input.spherical_projection.type == SWE::SphericalProjectionType::Enable) {
            DynRowVector<double> y_node(nnode);

            for (uint node_id = 0; node_id < dbound.data.get_nnode(); ++node_id) {
                y_node[node_id] = shape.nodal_coordinates[node_id][GlobalCoord::y];
            }

            DynRowVector<double> y_at_gp(ngp);

            y_at_gp = dbound.ComputeNodalUgp(y_node);

            double cos_phi_o = std::cos(problem_specific_input.spherical_projection.latitude_o * PI / 180.0);
            double R         = problem_specific_input.spherical_projection.R;

            for (uint gp = 0; gp < ngp; ++gp) {
                boundary.aux_at_gp[SWE::Auxiliaries::sp][gp] = cos_phi_o / std::cos(y_at_gp[gp] / R);
            }
        } else {
            for (uint gp = 0; gp < ngp; ++gp) {
                boundary.aux_at_gp[SWE::Auxiliaries::sp][gp] = 1.0;
            }
        }
    });

    // SOURCE TERMS INITIALIZE
    // Manning N bottom friction
    if (problem_specific_input.bottom_friction.type == SWE::BottomFrictionType::Manning) {
        if (!Utilities::file_exists(problem_specific_input.bottom_friction.manning_data_file)) {
            throw std::logic_error("Fatal Error: manning factor data file " +
                                   problem_specific_input.bottom_friction.manning_data_file + " was not found!\n");
        }

        std::ifstream manning_file(problem_specific_input.bottom_friction.manning_data_file);

        uint node_ID;
        double manning_n;
        std::map<uint, double> node_manning_n;

        std::string line;
        while (std::getline(manning_file, line)) {
            std::istringstream input_string(line);

            if (!(input_string >> node_ID >> manning_n))
                break;

            node_manning_n[node_ID] = manning_n;
        }

        mesh.CallForEachElement([&node_manning_n](auto& elt) {
            const std::vector<uint>& node_ID = elt.GetNodeID();

            for (uint node = 0; node < elt.data.get_nnode(); ++node) {
                elt.data.source.manning_n[node] = node_manning_n[node_ID[node]];
            }

            elt.data.source.manning = true;
            elt.data.source.g_manning_n_sq =
                Global::g *
                std::pow(std::accumulate(elt.data.source.manning_n.begin(), elt.data.source.manning_n.end(), 0.0) /
                             elt.data.source.manning_n.size(),
                         2);
        });
    }

    // Coriolis effect
    if (problem_specific_input.coriolis.type == SWE::CoriolisType::Enable) {
        mesh.CallForEachElement([&problem_specific_input](auto& elt) {
            double y_avg = 0;

            for (uint node = 0; node < elt.data.get_nnode(); ++node) {
                y_avg += elt.GetShape().nodal_coordinates[node][GlobalCoord::y];
            }

            y_avg /= elt.data.get_nnode();

            // If spherical projection is not specified we use default value for R which is R earth
            double latitude_avg = y_avg / problem_specific_input.spherical_projection.R;

            elt.data.source.coriolis_f = 2 * 7.2921E-5 * std::sin(latitude_avg);
        });
    }

    // WETTING-DRYING INITIALIZE
    if (SWE::PostProcessing::wetting_drying) {
        mesh.CallForEachElement([](auto& elt) {
            auto& shape = elt.GetShape();

            auto& state    = elt.data.state[0];
            auto& wd_state = elt.data.wet_dry_state;

            for (uint vrtx = 0; vrtx < elt.data.get_nvrtx(); ++vrtx) {
                wd_state.bath_at_vrtx[vrtx] = shape.nodal_coordinates[vrtx][GlobalCoord::z];
            }

            wd_state.bath_min = *std::min_element(wd_state.bath_at_vrtx.begin(), wd_state.bath_at_vrtx.end());

            for ( uint var = 0; var < n_variables; ++var ) {
                row(wd_state.q_lin, var) = elt.ProjectBasisToLinear(state.q[var]);
            }

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

                for ( uint var = 0; var < n_variables; ++var) {
                    state.q[var] = elt.ProjectLinearToBasis(row(wd_state.q_at_vrtx,var));
                    set_constant(state.rhs[var], 0.0);
                }
            }
        });
    }

    // SLOPE LIMIT INITIALIZE
    if (SWE::PostProcessing::slope_limiting) {
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

            StatVector<double, SWE::n_dimensions> median;

            for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                median[GlobalCoord::x] =
                    sl_state.midpts_coord[bound][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                median[GlobalCoord::y] =
                    sl_state.midpts_coord[bound][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

                sl_state.median[bound] = median / norm(median);
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

            StatMatrix<double, 2, 2> A;
            StatVector<double, 2> b;

            for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                uint element_1 = bound;
                uint element_2 = (bound + 1) % elt.data.get_nbound();

                if (!is_distributed(elt.GetBoundaryType()[element_1]) &&
                    !is_distributed(elt.GetBoundaryType()[element_2])) {
                    A(0, 0) = sl_state.baryctr_coord_neigh[element_1][GlobalCoord::x] -
                              sl_state.baryctr_coord[GlobalCoord::x];
                    A(1, 0) = sl_state.baryctr_coord_neigh[element_1][GlobalCoord::y] -
                              sl_state.baryctr_coord[GlobalCoord::y];
                    A(0, 1) = sl_state.baryctr_coord_neigh[element_2][GlobalCoord::x] -
                              sl_state.baryctr_coord[GlobalCoord::x];
                    A(1, 1) = sl_state.baryctr_coord_neigh[element_2][GlobalCoord::y] -
                              sl_state.baryctr_coord[GlobalCoord::y];

                    b[0] = sl_state.midpts_coord[bound][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                    b[1] = sl_state.midpts_coord[bound][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

                    sl_state.alpha[bound] = inverse(A) * b;
                    sl_state.r_sq[bound]  = sq_norm(b);
                }
            }
        });
    }
}

template <typename MeshType, typename ProblemSpecificInputType>
void initialize_data_parallel_pre_send(MeshType& mesh,
                                       const ProblemSpecificInputType& problem_specific_input,
                                       uint comm_type) {
    initialize_data_serial(mesh, problem_specific_input);

    if (SWE::PostProcessing::slope_limiting) {
        mesh.CallForEachDistributedBoundary([comm_type](auto& dbound) {
            auto& sl_state = dbound.data.slope_limit_state;

            // Construct message to exterior state
            std::vector<double> message;

            message.reserve(SWE::n_dimensions);

            for (uint dim = 0; dim < SWE::n_dimensions; ++dim) {
                message.push_back(sl_state.baryctr_coord[dim]);
            }

            // Set message to send buffer
            dbound.boundary_condition.exchanger.SetToSendBuffer(comm_type, message);
        });
    }
}

template <typename MeshType>
void initialize_data_parallel_post_receive(MeshType& mesh, uint comm_type) {
    if (SWE::PostProcessing::slope_limiting) {
        mesh.CallForEachDistributedBoundary([comm_type](auto& dbound) {
            auto& sl_state = dbound.data.slope_limit_state;

            std::vector<double> message;

            message.resize(SWE::n_dimensions);

            dbound.boundary_condition.exchanger.GetFromReceiveBuffer(comm_type, message);

            for (uint dim = 0; dim < SWE::n_dimensions; ++dim) {
                sl_state.baryctr_coord_neigh[dbound.bound_id][dim] = message[dim];
            }
        });

        mesh.CallForEachElement([](auto& elt) {
            auto& sl_state = elt.data.slope_limit_state;

            StatMatrix<double, 2, 2> A;
            StatVector<double, 2> b;

            for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
                uint element_1 = bound;
                uint element_2 = (bound + 1) % elt.data.get_nbound();

                if (is_distributed(elt.GetBoundaryType()[element_1]) ||
                    is_distributed(elt.GetBoundaryType()[element_2])) {
                    A(0, 0) = sl_state.baryctr_coord_neigh[element_1][GlobalCoord::x] -
                              sl_state.baryctr_coord[GlobalCoord::x];
                    A(1, 0) = sl_state.baryctr_coord_neigh[element_1][GlobalCoord::y] -
                              sl_state.baryctr_coord[GlobalCoord::y];
                    A(0, 1) = sl_state.baryctr_coord_neigh[element_2][GlobalCoord::x] -
                              sl_state.baryctr_coord[GlobalCoord::x];
                    A(1, 1) = sl_state.baryctr_coord_neigh[element_2][GlobalCoord::y] -
                              sl_state.baryctr_coord[GlobalCoord::y];

                    b[0] = sl_state.midpts_coord[bound][GlobalCoord::x] - sl_state.baryctr_coord[GlobalCoord::x];
                    b[1] = sl_state.midpts_coord[bound][GlobalCoord::y] - sl_state.baryctr_coord[GlobalCoord::y];

                    sl_state.alpha[bound] = inverse(A) * b;
                    sl_state.r_sq[bound]  = sq_norm(b);
                }
            }
        });
    }
}
}

#endif