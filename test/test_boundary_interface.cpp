#include "general_definitions.hpp"
#include "utilities/almost_equal.hpp"
#include "geometry/mesh_definitions.hpp"
#include "problem/SWE/swe_problem.hpp"

int main() {
    using Utilities::almost_equal;
    bool error_found = false;

    using MasterType  = Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>;
    using ShapeType   = Shape::StraightTriangle;
    using ElementType = Geometry::Element<2, MasterType, ShapeType, SWE::Data>;

    using RawBoundaryType = Geometry::RawBoundary<1, SWE::Data>;
    using BoundaryType    = Geometry::Boundary<1, Integration::GaussLegendre_1D, SWE::Data, SWE::BC::Land>;
    using InterfaceType   = Geometry::Interface<1, Integration::GaussLegendre_1D, SWE::Data, SWE::IS::Interface>;

    // make an equilateral triangle
    std::vector<Point<2>> vrtxs(3);
    vrtxs[0] = {-0.5, 0.};
    vrtxs[1] = {0.5, 0.};
    vrtxs[2] = {0, std::sqrt(3.) / 2.};

    MasterType master(10);
    ShapeType shape(vrtxs);

    ElementType triangle(
        0,
        master,
        vrtxs,
        std::vector<uint>{0, 0, 0},
        std::vector<uint>{DEFAULT_ID, DEFAULT_ID, DEFAULT_ID},
        std::vector<unsigned char>{
            SWE::BoundaryConditions::land, SWE::BoundaryConditions::land, SWE::BoundaryConditions::land});

    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>> raw_boundary;

    // generate boundaries
    triangle.CreateRawBoundaries(raw_boundary);

    std::vector<BoundaryType> boundaries;

    for (auto& rb : raw_boundary[SWE::BoundaryConditions::land]) {
        boundaries.emplace_back(BoundaryType(rb.second));
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
    std::vector<Point<1>> gp_1D = integ_1D.GetRule(21).second;
    std::vector<Point<2>> gp_bound;

    Array2D<double> F_vals_bound(2);
    F_vals_bound[0].resize(triangle.data.get_ngp_boundary(0));
    F_vals_bound[1].resize(triangle.data.get_ngp_boundary(0));

    Array2D<double> Fn_vals_bound(3);
    Fn_vals_bound[0].resize(triangle.data.get_ngp_boundary(0));
    Fn_vals_bound[1].resize(triangle.data.get_ngp_boundary(0));
    Fn_vals_bound[2].resize(triangle.data.get_ngp_boundary(0));

    for (uint n_bound = 0; n_bound < 3; n_bound++) {
        gp_bound = master.BoundaryToMasterCoordinates(n_bound, gp_1D);

        x = shape.InterpolateNodalValues({-0.5, 0.5, 0}, gp_bound);
        y = shape.InterpolateNodalValues({0, 0, std::sqrt(3.) / 2.}, gp_bound);

        for (uint gp = 0; gp < triangle.data.get_ngp_boundary(0); gp++) {
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
    std::vector<double> gp_vals(triangle.data.get_ngp_boundary(0));
    std::vector<double> unit(triangle.data.get_ngp_boundary(0), 1.0);

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

    // Check ComputeNodalUgp, ComputeBoundaryNodalUgp and Integration
    std::vector<double> nodal_vals = {1.0, 2.0, 3.0};
    std::vector<double> bound_nodal_vals(2);

    std::vector<double> nodal_vals_gp(gp_1D.size());
    std::vector<double> nodal_vals_gp_comp(gp_1D.size());

    for (uint n_bound = 0; n_bound < 3; n_bound++) {
        gp_bound = master.BoundaryToMasterCoordinates(n_bound, gp_1D);

        nodal_vals_gp = shape.InterpolateNodalValues(nodal_vals, gp_bound);  // Assuming this is correctly evaluated
        boundaries[n_bound].ComputeNodalUgp(nodal_vals, nodal_vals_gp_comp);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            if (!almost_equal(nodal_vals_gp[gp], nodal_vals_gp_comp[gp])) {
                error_found = true;

                std::cerr << "Error found in boundary in ComputeNodalUgp" << std::endl;
            }
        }

        bound_nodal_vals[0] = nodal_vals[(n_bound + 1) % 3];
        bound_nodal_vals[1] = nodal_vals[(n_bound + 2) % 3];

        nodal_vals_gp = shape.InterpolateBoundaryNodalValues(
            n_bound, bound_nodal_vals, gp_1D);  // Assuming this is correctly evaluated
        boundaries[n_bound].ComputeBoundaryNodalUgp(bound_nodal_vals, nodal_vals_gp_comp);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            if (!almost_equal(nodal_vals_gp[gp], nodal_vals_gp_comp[gp])) {
                error_found = true;

                std::cerr << "Error found in boundary in ComputeBoundaryNodalUgp" << std::endl;
            }
        }

        double exact_integration = (nodal_vals[(n_bound + 1) % 3] + nodal_vals[(n_bound + 2) % 3]) / 2.0;

        if (!almost_equal(exact_integration, boundaries[n_bound].Integration(nodal_vals_gp_comp))) {
            error_found = true;

            std::cerr << "Error found in boundary in Integration" << std::endl;
        }
    }

    // generate Interfaces
    std::vector<InterfaceType> interfaces;
    for (auto& rb : raw_boundary[SWE::BoundaryConditions::land]) {
        interfaces.emplace_back(InterfaceType(rb.second, rb.second));
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

    // Check ComputeNodalUgp, ComputeBoundaryNodalUgp and Integration
    for (uint n_intface = 0; n_intface < 3; n_intface++) {
        gp_bound = master.BoundaryToMasterCoordinates(n_intface, gp_1D);

        nodal_vals_gp = shape.InterpolateNodalValues(nodal_vals, gp_bound);
        interfaces[n_intface].ComputeNodalUgpIN(nodal_vals, nodal_vals_gp_comp);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            if (!almost_equal(nodal_vals_gp[gp], nodal_vals_gp_comp[gp])) {
                error_found = true;

                std::cerr << "Error found in boundary in ComputeNodalUgpIN" << std::endl;
            }
        }

        bound_nodal_vals[0] = nodal_vals[(n_intface + 1) % 3];
        bound_nodal_vals[1] = nodal_vals[(n_intface + 2) % 3];

        nodal_vals_gp = shape.InterpolateBoundaryNodalValues(
            n_intface, bound_nodal_vals, gp_1D);  // Assuming this is correctly evaluated
        interfaces[n_intface].ComputeBoundaryNodalUgpIN(bound_nodal_vals, nodal_vals_gp_comp);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            if (!almost_equal(nodal_vals_gp[gp], nodal_vals_gp_comp[gp])) {
                error_found = true;

                std::cerr << "Error found in boundary in ComputeBoundaryNodalUgpIN" << std::endl;
            }
        }

        double exact_integration = (nodal_vals[(n_intface + 1) % 3] + nodal_vals[(n_intface + 2) % 3]) / 2.0;

        if (!almost_equal(exact_integration, interfaces[n_intface].IntegrationIN(nodal_vals_gp_comp))) {
            error_found = true;

            std::cerr << "Error found in boundary in IntegrationIN" << std::endl;
        }

        interfaces[n_intface].ComputeNodalUgpEX(nodal_vals, nodal_vals_gp_comp);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            if (!almost_equal(nodal_vals_gp[gp], nodal_vals_gp_comp[gp])) {
                error_found = true;

                std::cerr << "Error found in boundary in ComputeNodalUgpEX" << std::endl;
            }
        }

        nodal_vals_gp = shape.InterpolateBoundaryNodalValues(
            n_intface, bound_nodal_vals, gp_1D);  // Assuming this is correctly evaluated
        interfaces[n_intface].ComputeBoundaryNodalUgpEX(bound_nodal_vals, nodal_vals_gp_comp);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            if (!almost_equal(nodal_vals_gp[gp], nodal_vals_gp_comp[gp])) {
                error_found = true;

                std::cerr << "Error found in boundary in ComputeBoundaryNodalUgpEX" << std::endl;
            }
        }

        if (!almost_equal(exact_integration, interfaces[n_intface].IntegrationEX(nodal_vals_gp_comp))) {
            error_found = true;

            std::cerr << "Error found in boundary in IntegrationEX" << std::endl;
        }
    }

    if (error_found) {
        return 1;
    }

    return 0;
}