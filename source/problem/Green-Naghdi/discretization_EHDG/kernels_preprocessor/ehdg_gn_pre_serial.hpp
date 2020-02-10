#ifndef EHDG_GN_PRE_SERIAL_HPP
#define EHDG_GN_PRE_SERIAL_HPP

#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_init_data.hpp"
#include "ehdg_gn_pre_dbath_serial.hpp"

namespace GN {
namespace EHDG {
void Problem::preprocessor_serial(ProblemDiscretizationType& discretization,
                                  ProblemGlobalDataType& global_data,
                                  const ProblemStepperType& stepper,
                                  const ProblemInputType& problem_specific_input) {
    SWE_SIM::Problem::preprocessor_serial(
        discretization, global_data, stepper.GetFirstStepper(), problem_specific_input);

    GN::initialize_data_serial(discretization.mesh);

    uint dc_global_dof_offset = 0;
    Problem::initialize_global_dc_problem_serial(discretization, dc_global_dof_offset);

    const uint n_global_dofs = (discretization.mesh.GetP() + 1) * GN::n_dimensions;
    global_data.w1_hat_w1_hat.resize(n_global_dofs * dc_global_dof_offset, n_global_dofs * dc_global_dof_offset);
    global_data.w1_hat_rhs.resize(n_global_dofs * dc_global_dof_offset);

#if defined(B_RECON_AVG) || defined(D_RECON_AVG)
    std::set<uint> nodeIDs;
    discretization.mesh.CallForEachElement(
        [&nodeIDs](auto& elt) { nodeIDs.insert(elt.GetNodeID().begin(), elt.GetNodeID().end()); });
    uint max_nodeID = *std::max_element(nodeIDs.begin(), nodeIDs.end());

    global_data.derivatives_at_node = DynVector<double>((max_nodeID + 1) * GN::n_dddbath_terms);
    std::vector<uint> node_mult((max_nodeID + 1), 0);

    discretization.mesh.CallForEachElement([&node_mult](auto& elt) {
        auto& derivative = elt.data.derivative;
        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            derivative.local_nodeID[node] = elt.GetNodeID()[node];
            ++node_mult[derivative.local_nodeID[node]];
        }
    });

    discretization.mesh.CallForEachElement([&node_mult](auto& elt) {
        auto& derivative = elt.data.derivative;
        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            derivative.node_mult[node] = node_mult[derivative.local_nodeID[node]];
        }
    });
#endif

    Problem::compute_bathymetry_derivatives_serial(discretization, global_data);

    uint n_stages = stepper.GetFirstStepper().GetNumStages() > stepper.GetSecondStepper().GetNumStages()
                        ? stepper.GetFirstStepper().GetNumStages()
                        : stepper.GetSecondStepper().GetNumStages();

    discretization.mesh.CallForEachElement([n_stages](auto& elt) { elt.data.resize(n_stages + 1); });
}
}
}

#endif