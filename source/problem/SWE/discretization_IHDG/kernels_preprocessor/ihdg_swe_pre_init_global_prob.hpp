#ifndef IHDG_SWE_PRE_INIT_GLOBAL_PROB_HPP
#define IHDG_SWE_PRE_INIT_GLOBAL_PROB_HPP

namespace SWE {
namespace IHDG {
template <typename SimulationType>
void Problem::initialize_global_problem(SimulationType* simulation) {
    uint local_dof_offset  = 0;
    uint global_dof_offset = 0;

    simulation->mesh.CallForEachElement([&local_dof_offset](auto& elt) {
        auto& internal = elt.data.internal;

        // Set offsets for global matrix construction
        internal.local_dof_offset = local_dof_offset;
        local_dof_offset += elt.data.get_ndof() * SWE::n_variables;

        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); bound_id++) {
            elt.data.boundary[bound_id].local_dof_offset = internal.local_dof_offset;
        }

        // Initialize delta_local and rhs_local containers
        internal.delta_local.resize(SWE::n_variables * elt.data.get_ndof(), SWE::n_variables * elt.data.get_ndof());
        internal.rhs_local.resize(SWE::n_variables * elt.data.get_ndof());
    });

    simulation->mesh_skeleton.CallForEachEdgeInterface([&global_dof_offset](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        // Set offsets for global matrix construction
        edge_internal.global_dof_offset = global_dof_offset;
        global_dof_offset += edge_int.edge_data.get_ndof() * SWE::n_variables;

        boundary_in.global_dof_offset = global_dof_offset;
        boundary_ex.global_dof_offset = global_dof_offset;

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

    simulation->mesh_skeleton.CallForEachEdgeBoundary([&global_dof_offset](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        // Set offsets for global matrix construction
        edge_internal.global_dof_offset = global_dof_offset;
        global_dof_offset += edge_bound.edge_data.get_ndof() * SWE::n_variables;

        boundary.global_dof_offset = global_dof_offset;

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

    simulation->delta_local_inv.resize(local_dof_offset, local_dof_offset);
    simulation->delta_hat_local.resize(local_dof_offset, global_dof_offset);
    simulation->rhs_local.resize(local_dof_offset);

    simulation->delta_global.resize(global_dof_offset, local_dof_offset);
    simulation->delta_hat_global.resize(global_dof_offset, global_dof_offset);
    simulation->rhs_global.resize(global_dof_offset);
}
}
}

#endif