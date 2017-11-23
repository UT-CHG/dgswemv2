#ifndef SWE_PRE_INIT_DATA_HPP
#define SWE_PRE_INIT_DATA_HPP

#include "swe_pre_init_wd_data.hpp"
#include "swe_pre_init_sl_data.hpp"

namespace SWE {
void Problem::initialize_data_kernel(ProblemMeshType& mesh,
                                     const MeshMetaData& mesh_data,
                                     const ProblemInputType& problem_specific_input) {
    mesh.CallForEachElement([](auto& elt) { elt.data.initialize(); });

    std::unordered_map<uint, std::vector<double>> bathymetry;

    std::vector<double> bathymetry_temp;
    for (const auto& elt : mesh_data.elements) {
        auto nodal_coordinates = mesh_data.GetNodalCoordinates(elt.first);

        for (auto& node_coordinate : nodal_coordinates) {
            bathymetry_temp.push_back(node_coordinate[GlobalCoord::z]);
        }

        bathymetry.insert({elt.first, bathymetry_temp});

        bathymetry_temp.clear();
    }

    mesh.CallForEachElement([&bathymetry](auto& elt) {
        uint id = elt.GetID();

        auto& state = elt.data.state[0];
        auto& internal = elt.data.internal;

        if (!bathymetry.count(id)) {
            throw std::logic_error("Error: could not find bathymetry for element with id: " + id);
        }
        
        state.bath = elt.L2Projection(bathymetry[id]);

        elt.ComputeUgp(state.bath, internal.bath_at_gp);

        elt.ComputeDUgp(GlobalCoord::x, state.bath, internal.bath_deriv_wrt_x_at_gp);
        elt.ComputeDUgp(GlobalCoord::y, state.bath, internal.bath_deriv_wrt_y_at_gp);

        auto ze_init = [](Point<2>& pt) { return SWE::ic_ze(0, pt); };

        state.ze = elt.L2Projection(ze_init);

        auto qx_init = [](Point<2>& pt) { return SWE::ic_qx(0, pt); };

        state.qx = elt.L2Projection(qx_init);

        auto qy_init = [](Point<2>& pt) { return SWE::ic_qy(0, pt); };

        state.qy = elt.L2Projection(qy_init);
    });

    mesh.CallForEachInterface([](auto& intface) {
        auto& state_in = intface.data_in.state[0];
        auto& state_ex = intface.data_ex.state[0];

        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        intface.ComputeUgpIN(state_in.bath, boundary_in.bath_at_gp);
        intface.ComputeUgpEX(state_ex.bath, boundary_ex.bath_at_gp);
    });

    mesh.CallForEachBoundary([](auto& bound) {
        auto& state = bound.data.state[0];
        auto& boundary = bound.data.boundary[bound.bound_id];

        bound.ComputeUgp(state.bath, boundary.bath_at_gp);
    });

    mesh.CallForEachDistributedBoundary([](auto& dbound) {
        dbound.ComputeUgp(dbound.data.state[0].bath, dbound.data.boundary[dbound.bound_id].bath_at_gp);
    });

    Problem::initialize_wd_data_kernel(mesh);

    Problem::initialize_sl_data_kernel(mesh);
}
}

#endif