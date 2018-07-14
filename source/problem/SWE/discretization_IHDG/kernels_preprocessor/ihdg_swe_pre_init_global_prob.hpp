#ifndef IHDG_SWE_PRE_INIT_GLOBAL_PROB_HPP
#define IHDG_SWE_PRE_INIT_GLOBAL_PROB_HPP

namespace SWE {
namespace IHDG {
template <typename SimulationType>
void Problem::initialize_global_problem(SimulationType* simulation) {
    uint dof_local  = simulation->mesh.GetNumberElements() * 9;  // hardcode for p=1
    uint dof_global = simulation->mesh_skeleton.GetNumberEdgeInterfaces() * 6 +
                      simulation->mesh_skeleton.GetNumberEdgeBoundaries() * 6;  // hardcode for p=1

    simulation->delta_local_inv.resize(dof_local, dof_local);
    simulation->delta_hat_local.resize(dof_local, dof_global);
    simulation->rhs_local.resize(dof_local);
    simulation->rhs_local_prev.resize(dof_local);

    simulation->delta_global.resize(dof_global, dof_local);
    simulation->delta_hat_global.resize(dof_global, dof_global);
    simulation->rhs_global.resize(dof_global);
    simulation->rhs_global_prev.resize(dof_global);

    simulation->global.resize(dof_global, dof_global);

    simulation->mesh.CallForEachElement([](auto& elt) {
        auto& internal = elt.data.internal;

        // Set IDs for global matrix construction
        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); bound_id++) {
            elt.data.boundary[bound_id].elt_ID = elt.GetID();
        }

        // Initialize delta_local and rhs_local containers
        internal.delta_local.resize(SWE::n_variables * elt.data.get_ndof(), SWE::n_variables * elt.data.get_ndof());
        internal.rhs_local.resize(SWE::n_variables * elt.data.get_ndof());
    });

    simulation->mesh_skeleton.CallForEachEdgeInterface([](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        // Set IDs for global matrix construction
        boundary_in.edg_ID = edge_int.GetID();
        boundary_ex.edg_ID = edge_int.GetID();

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof(),
                                              SWE::n_variables * edge_int.edge_data.get_ndof());
        edge_internal.rhs_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof());

        // Initialize delta_hat_local containers
        boundary_in.delta_hat_local.resize(SWE::n_variables * edge_int.interface.data_in.get_ndof(),
                                           SWE::n_variables * edge_int.edge_data.get_ndof());
        boundary_ex.delta_hat_local.resize(SWE::n_variables * edge_int.interface.data_ex.get_ndof(),
                                           SWE::n_variables * edge_int.edge_data.get_ndof());

        // Initialize delta_global containers
        boundary_in.delta_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof(),
                                        SWE::n_variables * edge_int.interface.data_in.get_ndof());
        boundary_ex.delta_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof(),
                                        SWE::n_variables * edge_int.interface.data_ex.get_ndof());
    });

    simulation->mesh_skeleton.CallForEachEdgeBoundary([](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        // Set IDs for global matrix construction
        boundary.edg_ID = edge_bound.GetID();

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof(),
                                              SWE::n_variables * edge_bound.edge_data.get_ndof());
        edge_internal.rhs_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof());

        // Initialize delta_hat_local container
        boundary.delta_hat_local.resize(SWE::n_variables * edge_bound.boundary.data.get_ndof(),
                                        SWE::n_variables * edge_bound.edge_data.get_ndof());

        // Initialize delta_global container
        boundary.delta_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof(),
                                     SWE::n_variables * edge_bound.boundary.data.get_ndof());
    });
}
}
}

#endif