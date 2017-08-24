#include "general_definitions.hpp"
#include "utilities/almost_equal.hpp"
#include "geometry/mesh_definitions.hpp"
#include "problem/SWE/swe_problem.hpp"

int main() {
    using Utilities::almost_equal;
    bool error_found = false;

    using MasterType = Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>;
    using ShapeType = Shape::StraightTriangle;
    using ElementType = Geometry::Element<2, MasterType, ShapeType, SWE::Data>;

    using RawBoundaryType = Geometry::RawBoundary<1, SWE::Data>;
    using BoundaryType = Geometry::Boundary<1, Integration::GaussLegendre_1D, SWE::Data, SWE::Land>;
    using InterfaceType = Geometry::Interface<1, Integration::GaussLegendre_1D, SWE::Data>;

    // make an equilateral triangle
    std::vector<Point<2>> vrtxs(3);
    vrtxs[0] = {-0.5, 0.};
    vrtxs[1] = {0.5, 0.};
    vrtxs[2] = {0, std::sqrt(3.) / 2.};

    MasterType master(10);
    ShapeType shape(vrtxs);

    ElementType triangle(0,
                         master,
                         vrtxs,
                         std::vector<uint>{DEFAULT_ID, DEFAULT_ID, DEFAULT_ID},
                         std::vector<unsigned char>{SWE::BoundaryConditions::land, SWE::BoundaryConditions::land,
                                                    SWE::BoundaryConditions::land});

    std::map<unsigned char, std::vector<RawBoundaryType>> pre_boundaries;
    std::map<uint, std::map<uint, RawBoundaryType>> pre_interfaces;
    std::map<uint, std::map<uint, RawBoundaryType>> pre_distributed_interfaces;

    // generate boundaries
    triangle.CreateRawBoundaries(pre_interfaces, pre_boundaries, pre_distributed_interfaces);

    std::vector<BoundaryType> boundaries;

    for (uint i = 0; i < 3; i++) {
        boundaries.emplace_back(BoundaryType(pre_boundaries[SWE::BoundaryConditions::land].at(i)));
    }

    // Check Integrations
    Integration::Dunavant_2D integ_2D;
    std::vector<Point<2>> gp_2D = integ_2D.GetRule(20).second;

    std::vector<double> x = shape.InterpolateNodalValues({-0.5, 0.5, 0}, gp_2D);
    std::vector<double> y = shape.InterpolateNodalValues({0, 0, std::sqrt(3.) / 2.}, gp_2D);

    Array2D<double> F_vals_int(2);
    F_vals_int[GlobalCoord::x].resize(triangle.data.get_ngp_internal());
    F_vals_int[GlobalCoord::y].resize(triangle.data.get_ngp_internal());

    std::vector<double> divF_vals_int(triangle.data.get_ngp_internal());

    for (uint gp = 0; gp < triangle.data.get_ngp_internal(); gp++) {
        F_vals_int[GlobalCoord::x][gp] = std::pow(x[gp], 2);
        F_vals_int[GlobalCoord::y][gp] = std::pow(y[gp], 2);

        divF_vals_int[gp] = 2 * x[gp] + 2 * y[gp];
    }

    Integration::GaussLegendre_1D integ_1D;
    std::vector<Point<1>> gp_1D = integ_1D.GetRule(20).second;
    std::vector<Point<2>> gp_bound;

    Array2D<double> F_vals_bound(2);
    F_vals_bound[0].resize(triangle.data.get_ngp_boundary());
    F_vals_bound[1].resize(triangle.data.get_ngp_boundary());

    Array2D<double> Fn_vals_bound(3);
    Fn_vals_bound[0].resize(triangle.data.get_ngp_boundary());
    Fn_vals_bound[1].resize(triangle.data.get_ngp_boundary());
    Fn_vals_bound[2].resize(triangle.data.get_ngp_boundary());

    for (uint n_bound = 0; n_bound < 3; n_bound++) {
        gp_bound = master.BoundaryToMasterCoordinates(n_bound, gp_1D);

        x = shape.InterpolateNodalValues({-0.5, 0.5, 0}, gp_bound);
        y = shape.InterpolateNodalValues({0, 0, std::sqrt(3.) / 2.}, gp_bound);

        for (uint gp = 0; gp < triangle.data.get_ngp_boundary(); gp++) {
            F_vals_bound[GlobalCoord::x][gp] = std::pow(x[gp], 2);
            F_vals_bound[GlobalCoord::y][gp] = std::pow(y[gp], 2);

            Fn_vals_bound[n_bound][gp] =
                F_vals_bound[GlobalCoord::x][gp] * boundaries[n_bound].surface_normal[gp][GlobalCoord::x] +
                F_vals_bound[GlobalCoord::y][gp] * boundaries[n_bound].surface_normal[gp][GlobalCoord::y];
        }
    }

    double divergence;

    for (uint dof = 0; dof < 66; dof++) {
        divergence = triangle.IntegrationPhi(dof, divF_vals_int) +
                     triangle.IntegrationDPhi(GlobalCoord::x, dof, F_vals_int[GlobalCoord::x]) +
                     triangle.IntegrationDPhi(GlobalCoord::y, dof, F_vals_int[GlobalCoord::y]) -
                     boundaries[0].IntegrationPhi(dof, Fn_vals_bound[0]) -
                     boundaries[1].IntegrationPhi(dof, Fn_vals_bound[1]) -
                     boundaries[2].IntegrationPhi(dof, Fn_vals_bound[2]);

        if (!almost_equal(divergence, 0.0, 1.e+03)) {
            error_found = true;

            std::cerr << "Error found in boundary in IntegrationPhi" << std::endl;
        }
    }

    std::vector<double> mod_vals(triangle.data.get_ndof());
    std::vector<double> gp_vals(triangle.data.get_ngp_boundary());
    std::vector<double> unit(triangle.data.get_ngp_boundary(), 1.0);

    // Check ComputeUgp
    for (uint dof = 0; dof < 66; dof++) {
        std::fill(mod_vals.begin(), mod_vals.end(), 0.0);
        mod_vals[dof] = 1.0;

        for (uint n_bound = 0; n_bound < 3; n_bound++) {
            boundaries[n_bound].ComputeUgp(mod_vals, gp_vals);
            if (!almost_equal(boundaries[n_bound].IntegrationPhi(dof, unit),
                              boundaries[n_bound].IntegrationPhi(0, gp_vals),
                              1.e+03)) {
                error_found = true;

                std::cerr << "Error found in boundary in ComputeUgp" << std::endl;
            }
        }
    }

    // generate Interfaces
    std::vector<InterfaceType> interfaces;

    for (uint i = 0; i < 3; i++) {
        interfaces.emplace_back(InterfaceType(pre_boundaries[SWE::BoundaryConditions::land].at(i),
                                              pre_boundaries[SWE::BoundaryConditions::land].at(i)));
    }

    // Check IntegrationPhiIN/EX
    for (uint dof = 0; dof < 66; dof++) {
        std::fill(mod_vals.begin(), mod_vals.end(), 0.0);
        mod_vals[dof] = 1.0;

        for (uint n_intface = 0; n_intface < 3; n_intface++) {
            if (!almost_equal(boundaries[n_intface].IntegrationPhi(dof, Fn_vals_bound[n_intface]),
                              interfaces[n_intface].IntegrationPhiIN(dof, Fn_vals_bound[n_intface]))) {
                error_found = true;

                std::cerr << "Error found in interface in IntegrationPhiIN" << std::endl;
            }

            if (!almost_equal(boundaries[n_intface].IntegrationPhi(dof, Fn_vals_bound[n_intface]),
                              interfaces[n_intface].IntegrationPhiEX(dof, Fn_vals_bound[n_intface]))) {
                error_found = true;

                std::cerr << "Error found in interface in IntegrationPhiEX" << std::endl;
            }
        }
    }

    // Check ComputeUgpIN/EX
    for (uint dof = 0; dof < 66; dof++) {
        std::fill(mod_vals.begin(), mod_vals.end(), 0.0);
        mod_vals[dof] = 1.0;

        for (uint n_intface = 0; n_intface < 3; n_intface++) {
            interfaces[n_intface].ComputeUgpIN(mod_vals, gp_vals);
            if (!almost_equal(interfaces[n_intface].IntegrationPhiIN(dof, unit),
                              interfaces[n_intface].IntegrationPhiIN(0, gp_vals),
                              1.e+03)) {
                error_found = true;

                std::cerr << "Error found in interface in ComputeUgpIN" << std::endl;
            }

            interfaces[n_intface].ComputeUgpEX(mod_vals, gp_vals);
            if (!almost_equal(interfaces[n_intface].IntegrationPhiEX(dof, unit),
                              interfaces[n_intface].IntegrationPhiEX(0, gp_vals),
                              1.e+03)) {
                error_found = true;

                std::cerr << "Error found in interface in ComputeUgpEX" << std::endl;
            }
        }
    }

    if (error_found) {
        return 1;
    }

    return 0;
}