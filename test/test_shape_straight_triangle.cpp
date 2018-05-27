#include "shape/shapes_2D.hpp"
#include "utilities/almost_equal.hpp"

int main() {
    using Utilities::almost_equal;
    bool error_found = false;

    // make an equilateral triangle
    std::vector<Point<2>> vrtxs(3);
    vrtxs[0] = {-0.5, 0.};
    vrtxs[1] = {0.5, 0.};
    vrtxs[2] = {0, std::sqrt(3.) / 2.};

    Shape::StraightTriangle triangle(vrtxs);

    // check Jdet
    double Jdet_true = std::sqrt(3.) / 8.;
    double Jdet_comp = triangle.GetJdet(std::vector<Point<2>>(0))[0];

    if (!almost_equal(Jdet_true, Jdet_comp)) {
        std::cerr << "Error GetJdet\n";
        error_found = true;
    }

    // check Jinv
    Array2D<double> Jinv_true = {{2., -2. / std::sqrt(3.)}, {0, 4. / std::sqrt(3.)}};
    Array3D<double> Jinv_comp = triangle.GetJinv(std::vector<Point<2>>(0));

    if (!almost_equal(Jinv_true[0][0], Jinv_comp[0][0][0]) || !almost_equal(Jinv_true[0][1], Jinv_comp[0][1][0]) ||
        !almost_equal(Jinv_true[1][0], Jinv_comp[1][0][0]) || !almost_equal(Jinv_true[1][1], Jinv_comp[1][1][0])) {
        std::cerr << "Error GetJinv\n";
        error_found = true;
    }

    // check SurfaceJ
    for (uint i = 0; i < 3; ++i) {
        double SurfaceJ_true = 0.5;
        double SurfaceJ_comp = triangle.GetSurfaceJ(i, std::vector<Point<2>>(0))[0];

        if (!almost_equal(SurfaceJ_true, SurfaceJ_comp)) {
            std::cerr << "Error GetSurfaceJ\n";
            error_found = true;
        }
    }

    // check SurfaceNormal
    Array2D<double> SurfaceNormal_true = {{std::sqrt(3.) / 2., 0.5}, {-std::sqrt(3.) / 2., 0.5}, {0, -1.}};

    for (uint i = 0; i < 3; ++i) {
        Array2D<double> SurfaceNormal_comp = triangle.GetSurfaceNormal(i, std::vector<Point<2>>(0));

        if (!almost_equal(SurfaceNormal_true[i][GlobalCoord::x], SurfaceNormal_comp[0][GlobalCoord::x]) ||
            !almost_equal(SurfaceNormal_true[i][GlobalCoord::y], SurfaceNormal_comp[0][GlobalCoord::y])) {
            std::cerr << "Error GetSurfaceNormal\n";
            error_found = true;
        }
    }

    // check Interpolation
    std::vector<double> nodal_vals = {-2., 2., 3.};

    std::vector<Point<2>> interpolation_pts = {
        {-1, -1}, {1, -1}, {-1, 1}, {0, 0}, {-1, 0}, {0, -1}, {-1 / 3., -1 / 3.}};

    std::vector<double> interpolation_true = {-2., 2., 3., 2.5, 0.5, 0, 1.};

    std::vector<double> interpolation_comp = triangle.InterpolateNodalValues(nodal_vals, interpolation_pts);

    for (uint i = 0; i < 6; i++) {
        if (!almost_equal(interpolation_true[i], interpolation_comp[i])) {
            std::cerr << "Error InterpolateNodalValues\n";
            error_found = true;
        }
    }

    // check Interpolation of boundary nodal values to boundaries
    std::vector<double> bound_nodal_vals(2);

    std::vector<Point<1>> bound_interpolation_pts = {{-1}, {-0.5}, {0}, {0.5}, {1}};

    Array2D<double> bound_interpolation_true = {
        {2., 2.25, 2.5, 2.75, 3.}, {3., 1.75, 0.5, -0.75, -2}, {-2., -1., 0., 1., 2.}};

    Array2D<double> bound_interpolation_comp;
    bound_interpolation_comp.resize(3);

    for (uint bound_id = 0; bound_id < 3; bound_id++) {
        bound_nodal_vals[0] = nodal_vals[(bound_id + 1) % 3];
        bound_nodal_vals[1] = nodal_vals[(bound_id + 2) % 3];

        bound_interpolation_comp[bound_id] =
            triangle.InterpolateBoundaryNodalValues(bound_id, bound_nodal_vals, bound_interpolation_pts);
    }

    for (uint bound_id = 0; bound_id < 3; bound_id++) {
        for (uint i = 0; i < 5; i++) {
            if (!almost_equal(bound_interpolation_true[bound_id][i], bound_interpolation_comp[bound_id][i])) {
                std::cerr << "Error InterpolateNodalValues\n";
                error_found = true;
            }
        }
    }

    // check Interpolation Derivative
    // take linear function f = x + 1/sqrt(3) * y
    nodal_vals = {-0.5, 0.5, 0.5};

    Array2D<double> interpolation_derivative_true = {{1.0}, {1.0 / std::sqrt(3.)}};

    Array2D<double> interpolation_derivative_comp =
        triangle.InterpolateNodalValuesDerivatives(nodal_vals, std::vector<Point<2>>{{0.0, 0.0}});

    if (!almost_equal(interpolation_derivative_true[0][0], interpolation_derivative_comp[0][0]) ||
        !almost_equal(interpolation_derivative_true[1][0], interpolation_derivative_comp[1][0])) {
        std::cerr << "Error InterpolateNodalValuesDerivatives\n";
        error_found = true;
    }

    if (error_found) {
        return 1;
    }

    return 0;
}