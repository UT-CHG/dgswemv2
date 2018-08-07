#ifndef EHDG_GN_PROC_SERIAL_STAGE_HPP
#define EHDG_GN_PROC_SERIAL_STAGE_HPP

#include "general_definitions.hpp"

#include "ehdg_gn_kernels_processor.hpp"

namespace GN {
namespace EHDG {
void Problem::serial_stage_kernel(const RKStepper& stepper, ProblemDiscretizationType& discretization) {
    Problem::serial_swe_stage_kernel(stepper, discretization);

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& curr_state = elt.data.state[stage];
        auto& next_state = elt.data.state[stage + 1];

        curr_state.q = next_state.q;
    });

    Problem::serial_dispersive_correction_kernel(stepper, discretization);

    Problem::serial_swe_stage_kernel(stepper, discretization);

    /*discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::update_kernel(stepper, elt); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        bool nan_found = Problem::scrutinize_solution_kernel(stepper, elt);

        if (nan_found)
            abort();
    });*/
}

void Problem::serial_swe_stage_kernel(const RKStepper& stepper, ProblemDiscretizationType& discretization) {
    /* Global Step */
    discretization.mesh.CallForEachInterface(
        [&stepper](auto& intface) { Problem::global_swe_interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary(
        [&stepper](auto& bound) { Problem::global_swe_boundary_kernel(stepper, bound); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&stepper](auto& edge_int) { Problem::global_swe_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&stepper](auto& edge_bound) { Problem::global_swe_edge_boundary_kernel(stepper, edge_bound); });
    /* Global Step */

    /* Local Step */
    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::local_swe_volume_kernel(stepper, elt); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::local_swe_source_kernel(stepper, elt); });

    discretization.mesh.CallForEachInterface(
        [&stepper](auto& intface) { Problem::local_swe_interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary(
        [&stepper](auto& bound) { Problem::local_swe_boundary_kernel(stepper, bound); });
    /* Local Step */

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& curr_state = elt.data.state[stage];
        auto& next_state = elt.data.state[stage + 1];

        curr_state.solution = elt.ApplyMinv(curr_state.rhs);

        next_state.q = curr_state.q + stepper.GetDT() / 2.0 * curr_state.solution;
    });
}

void Problem::serial_dispersive_correction_kernel(const RKStepper& stepper, ProblemDiscretizationType& discretization) {
    // Compute du, ddu
    Problem::serial_derivatives_kernel(stepper, discretization);

    // Compute dbath, ddbath, dddbath (if bath is not a variable then this goes into data initialization)
    Problem::serial_bathymetry_derivatives_kernel(stepper, discretization);

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::local_dc_volume_kernel(stepper, elt); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::local_dc_source_kernel(stepper, elt); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&stepper](auto& edge_int) { Problem::local_dc_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&stepper](auto& edge_bound) { Problem::local_dc_edge_boundary_kernel(stepper, edge_bound); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&stepper](auto& edge_int) { Problem::global_dc_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&stepper](auto& edge_bound) { Problem::global_dc_edge_boundary_kernel(stepper, edge_bound); });

    Problem::solve_global_dc_problem(stepper, discretization);

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;

        auto h = row(internal.aux_at_gp, GN::Auxiliaries::h);

        auto dze_dx = elt.ComputeDUgp(GlobalCoord::x, row(state.q, GN::Variables::ze));
        auto dze_dy = elt.ComputeDUgp(GlobalCoord::y, row(state.q, GN::Variables::ze));

        row(internal.source_at_gp, GN::Variables::qx) =
            Global::g / NDParameters::alpha * cwise_multiplication(dze_dx, h);
        row(internal.source_at_gp, GN::Variables::qy) =
            Global::g / NDParameters::alpha * cwise_multiplication(dze_dy, h);

        row(internal.source_at_gp, GN::Variables::qx) -= elt.ComputeUgp(row(state.w1, GlobalCoord::x));
        row(internal.source_at_gp, GN::Variables::qy) -= elt.ComputeUgp(row(state.w1, GlobalCoord::y));

        row(state.rhs, GN::Variables::qx) = elt.IntegrationPhi(row(internal.source_at_gp, GN::Variables::qx));
        row(state.rhs, GN::Variables::qy) = elt.IntegrationPhi(row(internal.source_at_gp, GN::Variables::qy));

        row(state.solution, GN::Variables::qx) = elt.ApplyMinv(row(state.rhs, GN::Variables::qx));
        row(state.solution, GN::Variables::qy) = elt.ApplyMinv(row(state.rhs, GN::Variables::qy));

        row(state.q, GN::Variables::qx) += stepper.GetDT() * row(state.solution, GN::Variables::qx);
        row(state.q, GN::Variables::qy) += stepper.GetDT() * row(state.solution, GN::Variables::qy);
    });
}
}
}

#endif