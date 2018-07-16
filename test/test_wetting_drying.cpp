#include "general_definitions.hpp"
#include "utilities/almost_equal.hpp"
#include "geometry/mesh_definitions.hpp"
#include "preprocessor/input_parameters.hpp"

#include "simulation/stepper/rk_stepper.hpp"

#include "problem/SWE/problem_function_files/swe_source_functions.hpp"
#include "problem/SWE/discretization_RKDG/rkdg_swe_problem.hpp"
#include "problem/SWE/discretization_RKDG/kernels_postprocessor/rkdg_swe_post_wet_dry.hpp"

int main() {
    using Utilities::almost_equal;
    bool error_found = false;

    using MasterType  = Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>;
    using ShapeType   = Shape::StraightTriangle;
    using ElementType = Geometry::Element<2, MasterType, ShapeType, SWE::RKDG::Data>;

    // the whole test is designed for h_0 = 0.01
    SWE::PostProcessing::h_o = 0.01;

    // make an equilateral triangle
    DynVector<Point<3>> vrtxs(3);
    vrtxs[0] = {-0.5, 0., 0.};
    vrtxs[1] = {0.5, 0., 0.};
    vrtxs[2] = {0, std::sqrt(3.) / 2., 0.};

    MasterType master(10);

    ElementType triangle(0,
                         master,
                         std::move(vrtxs),
                         std::move(DynVector<uint>(0)),
                         std::move(DynVector<uint>(0)),
                         std::move(DynVector<unsigned char>(0)));

    triangle.data.initialize();
    triangle.data.resize(2);

    StepperInput stepper_input;

    stepper_input.nstages = 1;
    stepper_input.order   = 1;
    stepper_input.dt      = 1;

    RKStepper stepper(stepper_input);

    auto& wd_state = triangle.data.wet_dry_state;
    auto& state    = triangle.data.state[1];

    wd_state.bath_at_vrtx[0] = 1.;
    wd_state.bath_at_vrtx[1] = 2.;
    wd_state.bath_at_vrtx[2] = 3.;

    wd_state.bath_min = 1.;

    // Completely dry element
    wd_state.q_at_vrtx[0][SWE::Variables::ze] = SWE::PostProcessing::h_o / 2.0 - wd_state.bath_at_vrtx[0];
    wd_state.q_at_vrtx[1][SWE::Variables::ze] = SWE::PostProcessing::h_o / 4.0 - wd_state.bath_at_vrtx[1];
    wd_state.q_at_vrtx[2][SWE::Variables::ze] = SWE::PostProcessing::h_o / 6.0 - wd_state.bath_at_vrtx[2];

    wd_state.q_at_vrtx[0][SWE::Variables::qx] = 1.;
    wd_state.q_at_vrtx[1][SWE::Variables::qx] = 2.;
    wd_state.q_at_vrtx[2][SWE::Variables::qx] = 3.;

    wd_state.q_at_vrtx[0][SWE::Variables::qy] = -1.;
    wd_state.q_at_vrtx[1][SWE::Variables::qy] = -2.;
    wd_state.q_at_vrtx[2][SWE::Variables::qy] = -3.;

    triangle.L2Projection(wd_state.q_at_vrtx, state.q);

    SWE::RKDG::Problem::wetting_drying_kernel(stepper, triangle);

    triangle.ProjectBasisToLinear(state.q, wd_state.q_lin);
    triangle.ComputeLinearUvrtx(wd_state.q_lin, wd_state.q_at_vrtx);

    if (wd_state.wet) {
        error_found = true;
        printf("Completely dry element is not dry!\n");
    }

    double h_avg =
        std::accumulate(wd_state.h_at_vrtx.begin(), wd_state.h_at_vrtx.end(), 0.0) / triangle.data.get_nvrtx();

    for (uint vrtx = 0; vrtx < triangle.data.get_nvrtx(); vrtx++) {
        if (!almost_equal(wd_state.q_at_vrtx[vrtx][SWE::Variables::ze], h_avg - wd_state.bath_at_vrtx[vrtx], 1.e+4)) {
            error_found = true;
            printf("Error in setting ze at vrtx %d. Set value: %f. Correct value: %f.\n",
                   vrtx,
                   wd_state.q_at_vrtx[vrtx][SWE::Variables::ze],
                   h_avg - wd_state.bath_at_vrtx[vrtx]);
        }
        if (!almost_equal(wd_state.q_at_vrtx[vrtx][SWE::Variables::qx], 0.0)) {
            error_found = true;
            printf("Error in setting qx at vrtx %d. Set value: %f. Correct value: %f.\n",
                   vrtx,
                   wd_state.q_at_vrtx[vrtx][SWE::Variables::qx],
                   0.0);
        }
        if (!almost_equal(wd_state.q_at_vrtx[vrtx][SWE::Variables::qy], 0.0)) {
            error_found = true;
            printf("Error in setting qy at vrtx %d. Set value: %f. Correct value: %f.\n",
                   vrtx,
                   wd_state.q_at_vrtx[vrtx][SWE::Variables::qy],
                   0.0);
        }
    }

    // Completely wet element
    wd_state.q_at_vrtx[0][SWE::Variables::ze] = SWE::PostProcessing::h_o;
    wd_state.q_at_vrtx[1][SWE::Variables::ze] = SWE::PostProcessing::h_o;
    wd_state.q_at_vrtx[2][SWE::Variables::ze] = SWE::PostProcessing::h_o;

    wd_state.q_at_vrtx[0][SWE::Variables::qx] = 1.;
    wd_state.q_at_vrtx[1][SWE::Variables::qx] = 2.;
    wd_state.q_at_vrtx[2][SWE::Variables::qx] = 3.;

    wd_state.q_at_vrtx[0][SWE::Variables::qy] = -1.;
    wd_state.q_at_vrtx[1][SWE::Variables::qy] = -2.;
    wd_state.q_at_vrtx[2][SWE::Variables::qy] = -3.;

    triangle.L2Projection(wd_state.q_at_vrtx, state.q);

    SWE::RKDG::Problem::wetting_drying_kernel(stepper, triangle);

    triangle.ProjectBasisToLinear(state.q, wd_state.q_lin);
    triangle.ComputeLinearUvrtx(wd_state.q_lin, wd_state.q_at_vrtx);

    if (!wd_state.wet) {
        error_found = true;
        printf("Completely wet element is not wet!\n");
    }

    for (uint vrtx = 0; vrtx < triangle.data.get_nvrtx(); vrtx++) {
        if (!almost_equal(wd_state.q_at_vrtx[vrtx][SWE::Variables::ze], SWE::PostProcessing::h_o, 1.e+4)) {
            error_found = true;
            printf("Error in setting ze at vrtx %d. Set value: %f. Correct value: %f.\n",
                   vrtx,
                   wd_state.q_at_vrtx[vrtx][SWE::Variables::ze],
                   SWE::PostProcessing::h_o);
        }
        if (!almost_equal(wd_state.q_at_vrtx[vrtx][SWE::Variables::qx], (double)(vrtx + 1), 1.e+4)) {
            error_found = true;
            printf("Error in setting qx at vrtx %d. Set value: %f. Correct value: %f.\n",
                   vrtx,
                   wd_state.q_at_vrtx[vrtx][SWE::Variables::qx],
                   (double)(vrtx + 1));
        }
        if (!almost_equal(wd_state.q_at_vrtx[vrtx][SWE::Variables::qy], -(double)(vrtx + 1), 1.e+4)) {
            error_found = true;
            printf("Error in setting qy at vrtx %d. Set value: %f. Correct value: %f.\n",
                   vrtx,
                   wd_state.q_at_vrtx[vrtx][SWE::Variables::qy],
                   -(double)(vrtx + 1));
        }
    }

    // Dry element on average
    wd_state.q_at_vrtx[0][SWE::Variables::ze] = SWE::PostProcessing::h_o / 2.0 - wd_state.bath_at_vrtx[0];
    wd_state.q_at_vrtx[1][SWE::Variables::ze] = SWE::PostProcessing::h_o / 2.0 - wd_state.bath_at_vrtx[1];
    wd_state.q_at_vrtx[2][SWE::Variables::ze] = 2.0 * SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[2];

    wd_state.q_at_vrtx[0][SWE::Variables::qx] = 1.;
    wd_state.q_at_vrtx[1][SWE::Variables::qx] = 2.;
    wd_state.q_at_vrtx[2][SWE::Variables::qx] = 3.;

    wd_state.q_at_vrtx[0][SWE::Variables::qy] = -1.;
    wd_state.q_at_vrtx[1][SWE::Variables::qy] = -2.;
    wd_state.q_at_vrtx[2][SWE::Variables::qy] = -3.;

    triangle.L2Projection(wd_state.q_at_vrtx, state.q);

    SWE::RKDG::Problem::wetting_drying_kernel(stepper, triangle);

    triangle.ProjectBasisToLinear(state.q, wd_state.q_lin);
    triangle.ComputeLinearUvrtx(wd_state.q_lin, wd_state.q_at_vrtx);

    if (wd_state.wet) {
        error_found = true;
        printf("Dry on average element is not dry!\n");
    }

    h_avg = std::accumulate(wd_state.h_at_vrtx.begin(), wd_state.h_at_vrtx.end(), 0.0) / triangle.data.get_nvrtx();

    for (uint vrtx = 0; vrtx < triangle.data.get_nvrtx(); vrtx++) {
        if (!almost_equal(wd_state.q_at_vrtx[vrtx][SWE::Variables::ze], h_avg - wd_state.bath_at_vrtx[vrtx], 1.e+4)) {
            error_found = true;
            printf("Error in setting ze at vrtx %d. Set value: %f. Correct value: %f.\n",
                   vrtx,
                   wd_state.q_at_vrtx[vrtx][SWE::Variables::ze],
                   h_avg - wd_state.bath_at_vrtx[vrtx]);
        }
        if (!almost_equal(wd_state.q_at_vrtx[vrtx][SWE::Variables::qx], 0.0)) {
            error_found = true;
            printf("Error in setting qx at vrtx %d. Set value: %f. Correct value: %f.\n",
                   vrtx,
                   wd_state.q_at_vrtx[vrtx][SWE::Variables::qx],
                   0.0);
        }
        if (!almost_equal(wd_state.q_at_vrtx[vrtx][SWE::Variables::qy], 0.0)) {
            error_found = true;
            printf("Error in setting qy at vrtx %d. Set value: %f. Correct value: %f.\n",
                   vrtx,
                   wd_state.q_at_vrtx[vrtx][SWE::Variables::qy],
                   0.0);
        }
    }

    // Some nodes dry element of flood-type
    wd_state.wet = false;  // To do check_element

    wd_state.q_at_vrtx[0][SWE::Variables::ze] = SWE::PostProcessing::h_o / 2.0 - wd_state.bath_at_vrtx[0];
    wd_state.q_at_vrtx[1][SWE::Variables::ze] = SWE::PostProcessing::h_o / 2.0 - wd_state.bath_at_vrtx[1];
    wd_state.q_at_vrtx[2][SWE::Variables::ze] = 3.5 * SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[2];

    wd_state.q_at_vrtx[0][SWE::Variables::qx] = 1.;
    wd_state.q_at_vrtx[1][SWE::Variables::qx] = 2.;
    wd_state.q_at_vrtx[2][SWE::Variables::qx] = 3.;

    wd_state.q_at_vrtx[0][SWE::Variables::qy] = -1.;
    wd_state.q_at_vrtx[1][SWE::Variables::qy] = -2.;
    wd_state.q_at_vrtx[2][SWE::Variables::qy] = -3.;

    triangle.L2Projection(wd_state.q_at_vrtx, state.q);

    SWE::RKDG::Problem::wetting_drying_kernel(stepper, triangle);

    triangle.ProjectBasisToLinear(state.q, wd_state.q_lin);
    triangle.ComputeLinearUvrtx(wd_state.q_lin, wd_state.q_at_vrtx);

    if (wd_state.wet) {
        error_found = true;
        printf("flood-type element is not dry!\n");
    }

    if (!almost_equal(
            wd_state.q_at_vrtx[0][SWE::Variables::ze], SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[0], 1.e+5)) {
        error_found = true;
        printf("Error in setting ze at vrtx %d. Set value: %f. Correct value: %f.\n",
               0,
               wd_state.q_at_vrtx[0][SWE::Variables::ze],
               SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[0]);
    }
    if (!almost_equal(
            wd_state.q_at_vrtx[1][SWE::Variables::ze], SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[1], 1.e+5)) {
        error_found = true;
        printf("Error in setting ze at vrtx %d. Set value: %f. Correct value: %f.\n",
               1,
               wd_state.q_at_vrtx[1][SWE::Variables::ze],
               SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[1]);
    }
    if (!almost_equal(wd_state.q_at_vrtx[2][SWE::Variables::ze],
                      2.5 * SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[2],
                      1.e+5)) {
        error_found = true;
        printf("Error in setting ze at vrtx %d. Set value: %f. Correct value: %f.\n",
               2,
               wd_state.q_at_vrtx[2][SWE::Variables::ze],
               2.5 * SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[2]);
    }

    for (uint vrtx = 0; vrtx < triangle.data.get_nvrtx(); vrtx++) {
        if (!almost_equal(wd_state.q_at_vrtx[vrtx][SWE::Variables::qx], 0.0)) {
            error_found = true;
            printf("Error in setting qx at vrtx %d. Set value: %f. Correct value: %f.\n",
                   vrtx,
                   wd_state.q_at_vrtx[vrtx][SWE::Variables::qx],
                   0.0);
        }
        if (!almost_equal(wd_state.q_at_vrtx[vrtx][SWE::Variables::qy], 0.0)) {
            error_found = true;
            printf("Error in setting qy at vrtx %d. Set value: %f. Correct value: %f.\n",
                   vrtx,
                   wd_state.q_at_vrtx[vrtx][SWE::Variables::qy],
                   0.0);
        }
    }

    // Some nodes dry element of dam-break-type
    wd_state.wet = false;  // To do check_element

    wd_state.q_at_vrtx[0][SWE::Variables::ze] = 3.5 * SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[0];
    wd_state.q_at_vrtx[1][SWE::Variables::ze] = SWE::PostProcessing::h_o / 2.0 - wd_state.bath_at_vrtx[1];
    wd_state.q_at_vrtx[2][SWE::Variables::ze] = SWE::PostProcessing::h_o / 2.0 - wd_state.bath_at_vrtx[2];

    wd_state.q_at_vrtx[0][SWE::Variables::qx] = 1.;
    wd_state.q_at_vrtx[1][SWE::Variables::qx] = 2.;
    wd_state.q_at_vrtx[2][SWE::Variables::qx] = 3.;

    wd_state.q_at_vrtx[0][SWE::Variables::qy] = -1.;
    wd_state.q_at_vrtx[1][SWE::Variables::qy] = -2.;
    wd_state.q_at_vrtx[2][SWE::Variables::qy] = -3.;

    triangle.L2Projection(wd_state.q_at_vrtx, state.q);

    SWE::RKDG::Problem::wetting_drying_kernel(stepper, triangle);

    triangle.ProjectBasisToLinear(state.q, wd_state.q_lin);
    triangle.ComputeLinearUvrtx(wd_state.q_lin, wd_state.q_at_vrtx);

    if (!wd_state.wet) {
        error_found = true;
        printf("dam-break-type element is not wet!\n");
    }

    if (!almost_equal(wd_state.q_at_vrtx[0][SWE::Variables::ze],
                      2.5 * SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[0],
                      1.e+5)) {
        error_found = true;
        printf("Error in setting ze at vrtx %d. Set value: %f. Correct value: %f.\n",
               0,
               wd_state.q_at_vrtx[0][SWE::Variables::ze],
               2.5 * SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[0]);
    }
    if (!almost_equal(
            wd_state.q_at_vrtx[1][SWE::Variables::ze], SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[1], 1.e+5)) {
        error_found = true;
        printf("Error in setting ze at vrtx %d. Set value: %f. Correct value: %f.\n",
               1,
               wd_state.q_at_vrtx[1][SWE::Variables::ze],
               SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[1]);
    }
    if (!almost_equal(
            wd_state.q_at_vrtx[2][SWE::Variables::ze], SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[2], 1.e+5)) {
        error_found = true;
        printf("Error in setting ze at vrtx %d. Set value: %f. Correct value: %f.\n",
               2,
               wd_state.q_at_vrtx[2][SWE::Variables::ze],
               SWE::PostProcessing::h_o - wd_state.bath_at_vrtx[2]);
    }

    if (!almost_equal(wd_state.q_at_vrtx[0][SWE::Variables::qx], 6.0, 1.e+5)) {
        error_found = true;
        printf("Error in setting qx at vrtx %d. Set value: %f. Correct value: %f.\n",
               0,
               wd_state.q_at_vrtx[0][SWE::Variables::qx],
               6.0);
    }
    if (!almost_equal(wd_state.q_at_vrtx[1][SWE::Variables::qx], 0.0, 1.e+5)) {
        error_found = true;
        printf("Error in setting qx at vrtx %d. Set value: %f. Correct value: %f.\n",
               1,
               wd_state.q_at_vrtx[1][SWE::Variables::qx],
               0.0);
    }
    if (!almost_equal(wd_state.q_at_vrtx[2][SWE::Variables::qx], 0.0, 1.e+5)) {
        error_found = true;
        printf("Error in setting qx at vrtx %d. Set value: %f. Correct value: %f.\n",
               2,
               wd_state.q_at_vrtx[2][SWE::Variables::qx],
               0.0);
    }

    if (!almost_equal(wd_state.q_at_vrtx[0][SWE::Variables::qy], -6.0, 1.e+5)) {
        error_found = true;
        printf("Error in setting qy at vrtx %d. Set value: %f. Correct value: %f.\n",
               0,
               wd_state.q_at_vrtx[0][SWE::Variables::qy],
               -6.0);
    }
    if (!almost_equal(wd_state.q_at_vrtx[1][SWE::Variables::qy], 0.0, 1.e+5)) {
        error_found = true;
        printf("Error in setting qy at vrtx %d. Set value: %f. Correct value: %f.\n",
               1,
               wd_state.q_at_vrtx[1][SWE::Variables::qy],
               0.0);
    }
    if (!almost_equal(wd_state.q_at_vrtx[2][SWE::Variables::qy], 0.0, 1.e+5)) {
        error_found = true;
        printf("Error in setting qy at vrtx %d. Set value: %f. Correct value: %f.\n",
               2,
               wd_state.q_at_vrtx[2][SWE::Variables::qy],
               0.0);
    }

    if (error_found) {
        return 1;
    }

    return 0;
}