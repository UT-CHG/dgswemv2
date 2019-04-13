#include "shape/shapes_2D.hpp"
#include "utilities/almost_equal.hpp"

int main() {
    using Utilities::almost_equal;
    bool error_found = false;

    // make an equilateral triangle
    AlignedVector<Point<3>> vrtxs(3);
    vrtxs[0] = {-0.5, 0., 0.};
    vrtxs[1] = {0.5, 0., 0.};
    vrtxs[2] = {0, std::sqrt(3.) / 2., 0.};

    Shape::StraightTriangle triangle(std::move(vrtxs));

    // check Jdet
    double Jdet_true = std::sqrt(3.) / 8.;
    double Jdet_comp = triangle.GetJdet(AlignedVector<Point<2>>(0))[0];

    if (!almost_equal(Jdet_true, Jdet_comp)) {
        std::cerr << "Error GetJdet\n";
        error_found = true;
    }

    // check Jinv
    StatMatrix<double, 2, 2> Jinv_true;
    Jinv_true(0, 0) = 2.;
    Jinv_true(0, 1) = -2. / std::sqrt(3.);
    Jinv_true(1, 0) = 0.0;
    Jinv_true(1, 1) = 4. / std::sqrt(3.);

    StatMatrix<double, 2, 2> Jinv_comp = triangle.GetJinv(AlignedVector<Point<2>>(0))[0];

    if (!almost_equal(Jinv_true(0, 0), Jinv_comp(0, 0)) || !almost_equal(Jinv_true(0, 1), Jinv_comp(0, 1)) ||
        !almost_equal(Jinv_true(1, 0), Jinv_comp(1, 0)) || !almost_equal(Jinv_true(1, 1), Jinv_comp(1, 1))) {
        std::cerr << "Error GetJinv\n";
        error_found = true;
    }

    // check SurfaceJ
    for (uint i = 0; i < 3; ++i) {
        double SurfaceJ_true = 0.5;
        double SurfaceJ_comp = triangle.GetSurfaceJ(i, AlignedVector<Point<2>>(0))[0];

        if (!almost_equal(SurfaceJ_true, SurfaceJ_comp)) {
            std::cerr << "Error GetSurfaceJ\n";
            error_found = true;
        }
    }

    // check SurfaceNormal
    AlignedVector<StatVector<double, 2>> SurfaceNormal_true(3);

    SurfaceNormal_true[0][GlobalCoord::x] = std::sqrt(3.) / 2.;
    SurfaceNormal_true[0][GlobalCoord::y] = 0.5;

    SurfaceNormal_true[1][GlobalCoord::x] = -std::sqrt(3.) / 2.;
    SurfaceNormal_true[1][GlobalCoord::y] = 0.5;

    SurfaceNormal_true[2][GlobalCoord::x] = 0.0;
    SurfaceNormal_true[2][GlobalCoord::y] = -1.0;

    for (uint i = 0; i < 3; ++i) {
        StatVector<double, 2> SurfaceNormal_comp = triangle.GetSurfaceNormal(i, AlignedVector<Point<2>>(0))[0];

        if (!almost_equal(SurfaceNormal_true[i][GlobalCoord::x], SurfaceNormal_comp[GlobalCoord::x]) ||
            !almost_equal(SurfaceNormal_true[i][GlobalCoord::y], SurfaceNormal_comp[GlobalCoord::y])) {
            std::cerr << "Error GetSurfaceNormal\n";
            error_found = true;
        }
    }

    // check GetPsi
    std::vector<double> nodal_vals = {-2., 2., 3.};

    AlignedVector<Point<2>> interpolation_pts(7);
    interpolation_pts[0] = {-1, -1};
    interpolation_pts[1] = {1, -1};
    interpolation_pts[2] = {-1, 1};
    interpolation_pts[3] = {0, 0};
    interpolation_pts[4] = {-1, 0};
    interpolation_pts[5] = {0, -1};
    interpolation_pts[6] = {-1 / 3., -1 / 3.};

    std::vector<double> interpolation_true = {-2., 2., 3., 2.5, 0.5, 0, 1.};

    DynMatrix<double> psi_interp = triangle.GetPsi(interpolation_pts);

    std::vector<double> interpolation_comp(7);

    for (uint i = 0; i < 6; ++i) {
        interpolation_comp[i] =
            psi_interp(0, i) * nodal_vals[0] + psi_interp(1, i) * nodal_vals[1] + psi_interp(2, i) * nodal_vals[2];

        if (!almost_equal(interpolation_true[i], interpolation_comp[i])) {
            std::cerr << "Error in GetPsi\n";
            error_found = true;
        }
    }

    // check GetBoudaryPsi
    std::vector<double> bound_nodal_vals(2);

    AlignedVector<Point<1>> bound_interpolation_pts(5);
    bound_interpolation_pts[0] = Point<1>{-1};
    bound_interpolation_pts[1] = Point<1>{-0.5};
    bound_interpolation_pts[2] = Point<1>{0};
    bound_interpolation_pts[3] = Point<1>{0.5};
    bound_interpolation_pts[4] = Point<1>{1};

    Array2D<double> bound_interpolation_true = {
        {2., 2.25, 2.5, 2.75, 3.}, {3., 1.75, 0.5, -0.75, -2}, {-2., -1., 0., 1., 2.}};

    Array2D<double> bound_interpolation_comp(3, std::vector<double>(5));

    DynMatrix<double> psi_bound_interp;

    for (uint bound_id = 0; bound_id < 3; ++bound_id) {
        bound_nodal_vals[0] = nodal_vals[(bound_id + 1) % 3];
        bound_nodal_vals[1] = nodal_vals[(bound_id + 2) % 3];

        psi_bound_interp = triangle.GetBoundaryPsi(bound_id, bound_interpolation_pts);

        for (uint i = 0; i < 5; ++i) {
            bound_interpolation_comp[bound_id][i] =
                psi_bound_interp(0, i) * bound_nodal_vals[0] + psi_bound_interp(1, i) * bound_nodal_vals[1];
        }
    }

    for (uint bound_id = 0; bound_id < 3; ++bound_id) {
        for (uint i = 0; i < 5; ++i) {
            if (!almost_equal(bound_interpolation_true[bound_id][i], bound_interpolation_comp[bound_id][i])) {
                std::cerr << "Error in GetBoundaryPsi\n";
                error_found = true;
            }
        }
    }

    // check GetDPsi
    // take linear function f = x + 1/sqrt(3) * y
    nodal_vals = {-0.5, 0.5, 0.5};

    Array2D<double> interpolation_derivative_true = {{1.0}, {1.0 / std::sqrt(3.)}};

    std::array<DynMatrix<double>, 2> dpsi_interp = triangle.GetDPsi(AlignedVector<Point<2>>(1, Point<2>{0, 0}));

    Array2D<double> interpolation_derivative_comp(2, std::vector<double>(1));

    interpolation_derivative_comp[0][0] = dpsi_interp[0](0, 0) * nodal_vals[0] + dpsi_interp[0](1, 0) * nodal_vals[1] +
                                          dpsi_interp[0](2, 0) * nodal_vals[2];

    interpolation_derivative_comp[1][0] = dpsi_interp[1](0, 0) * nodal_vals[0] + dpsi_interp[1](1, 0) * nodal_vals[1] +
                                          dpsi_interp[1](2, 0) * nodal_vals[2];

    if (!almost_equal(interpolation_derivative_true[0][0], interpolation_derivative_comp[0][0]) ||
        !almost_equal(interpolation_derivative_true[1][0], interpolation_derivative_comp[1][0])) {
        std::cerr << "Error in GetDPsi\n";
        error_found = true;
    }

    // Check Contains Point
    if (!triangle.ContainsPoint(Point<2>{0.0, 0.5})) {
        std::cerr << "Error in ContainsPoint\n";
        error_found = true;
    }

    if (!triangle.ContainsPoint(Point<2>{0.5, 0.0})) {
        std::cerr << "Error in ContainsPoint\n";
        error_found = true;
    }

    if (triangle.ContainsPoint(Point<2>{0.501, 0.0})) {
        std::cerr << "Error in ContainsPoint\n";
        error_found = true;
    }

    // Check global/local transform of coordinates
    AlignedVector<Point<2>> transformation_pts(4);
    transformation_pts[0] = {-1. / 3., -1. / 3.};  // barycenter
    transformation_pts[1] = {0., 0.};              // midpt 0
    transformation_pts[2] = {-1., 0.};             // midpt 1
    transformation_pts[3] = {0., -1.};             // midpt 2

    AlignedVector<Point<2>> transformation_true(4);
    transformation_true[0] = {0.0, std::sqrt(3.) / 6.};    // barycenter
    transformation_true[1] = {0.25, std::sqrt(3.) / 4.};   // midpt 0
    transformation_true[2] = {-0.25, std::sqrt(3.) / 4.};  // midpt 1
    transformation_true[3] = {0., 0.};                     // midpt 2

    AlignedVector<Point<2>> transformation = triangle.LocalToGlobalCoordinates(transformation_pts);

    for (uint pt = 0; pt < 4; ++pt) {
        for (uint dir = 0; dir < 2; ++dir) {
            if (!almost_equal(transformation[pt][dir], transformation_true[pt][dir])) {
                std::cerr << "Error in LocalToGlobalCoordinates\n";
                std::cerr << transformation[pt][dir] << ' ' << transformation_true[pt][dir] << '\n';
                error_found = true;
            }
        }
    }

    transformation = triangle.GlobalToLocalCoordinates(transformation);

    for (uint pt = 0; pt < 4; ++pt) {
        for (uint dir = 0; dir < 2; ++dir) {
            if (!almost_equal(transformation[pt][dir], transformation_pts[pt][dir])) {
                std::cerr << "Error in LocalToGlobalCoordinates\n";
                std::cerr << transformation[pt][dir] << ' ' << transformation_pts[pt][dir] << '\n';
                error_found = true;
            }
        }
    }

    if (error_found) {
        return 1;
    }

    return 0;
}