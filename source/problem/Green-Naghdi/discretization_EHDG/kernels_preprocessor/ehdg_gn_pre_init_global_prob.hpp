#ifndef EHDG_GN_PRE_INIT_GLOBAL_PROB_HPP
#define EHDG_GN_PRE_INIT_GLOBAL_PROB_HPP

namespace GN {
namespace EHDG {
void Problem::initialize_global_problem(HDGDiscretization<Problem>* discretization) {
    discretization->mesh.CallForEachElement([](auto& elt) {
        auto& internal = elt.data.internal;

        // Initialize w1 containers
        internal.w1_w1.resize(GN::n_dimensions * elt.data.get_ndof(), GN::n_dimensions * elt.data.get_ndof());
        internal.w1_w2.resize(GN::n_dimensions * elt.data.get_ndof(), elt.data.get_ndof());
        internal.w1_rhs.resize(GN::n_dimensions * elt.data.get_ndof());

        // Initialize w2 containers
        internal.w2_w1.resize(elt.data.get_ndof(), GN::n_dimensions * elt.data.get_ndof());
        internal.w2_w2.resize(elt.data.get_ndof(), elt.data.get_ndof());
    });

    discretization->mesh_skeleton.CallForEachEdgeInterface([](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(GN::n_variables * edge_int.edge_data.get_ndof(),
                                              GN::n_variables * edge_int.edge_data.get_ndof());
        edge_internal.rhs_global.resize(GN::n_variables * edge_int.edge_data.get_ndof());
    });

    discretization->mesh_skeleton.CallForEachEdgeBoundary([](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(GN::n_variables * edge_bound.edge_data.get_ndof(),
                                              GN::n_variables * edge_bound.edge_data.get_ndof());
        edge_internal.rhs_global.resize(GN::n_variables * edge_bound.edge_data.get_ndof());
    });

    discretization->mesh_skeleton.CallForEachEdgeDistributed([](auto& edge_dbound) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(GN::n_variables * edge_dbound.edge_data.get_ndof(),
                                              GN::n_variables * edge_dbound.edge_data.get_ndof());
        edge_internal.rhs_global.resize(GN::n_variables * edge_dbound.edge_data.get_ndof());
    });
}
}
}

#endif