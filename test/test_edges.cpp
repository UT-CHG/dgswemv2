#include "general_definitions.hpp"
#include "utilities/almost_equal.hpp"
#include "geometry/mesh_definitions.hpp"
#include "problem/SWE/discretization_EHDG/ehdg_swe_problem.hpp"

int main() {
    // Warning!!!
    // This test will not work unless you use Legendre_1D basis for your Edge elements!!!

    using Utilities::almost_equal;
    bool error_found = false;

    using MasterType  = Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>;
    using ShapeType   = Shape::StraightTriangle;
    using ElementType = Geometry::Element<2, MasterType, ShapeType, SWE::EHDG::Data>;

    using RawBoundaryType = Geometry::RawBoundary<1, SWE::EHDG::Data>;
    using BoundaryType    = Geometry::Boundary<1, Integration::GaussLegendre_1D, SWE::EHDG::Data, SWE::EHDG::BC::Land>;
    using InterfaceType =
        Geometry::Interface<1, Integration::GaussLegendre_1D, SWE::EHDG::Data, SWE::EHDG::IS::Internal>;
    using EdgeInternalType =
        typename std::tuple_element<0, Geometry::EdgeInternalTypeTuple<SWE::EHDG::Data, SWE::EHDG::EdgeData>>::type;
    using EdgeBoundaryType =
        typename std::tuple_element<0, Geometry::EdgeBoundaryTypeTuple<SWE::EHDG::Data, SWE::EHDG::EdgeData>>::type;

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
        std::vector<unsigned char>{SWE::BoundaryTypes::land, SWE::BoundaryTypes::land, SWE::BoundaryTypes::land});

    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>> raw_boundary;

    triangle.CreateRawBoundaries(raw_boundary);

    Integration::GaussLegendre_1D integration;

    std::vector<double> u(11);                           // ndof for p = 10
    std::vector<double> u_gp(integration.GetNumGP(21));  // ngp for 2*p+1

    // generate edge boundaries
    std::vector<BoundaryType> boundaries;

    for (auto& rb : raw_boundary[SWE::BoundaryTypes::land]) {
        boundaries.emplace_back(BoundaryType(rb.second));
    }

    std::vector<EdgeBoundaryType> edge_boundaries;

    for (auto& bound : boundaries) {
        edge_boundaries.emplace_back(EdgeBoundaryType(bound));
    }

    for (uint dof = 0; dof < 11; dof++) {
        std::fill(u.begin(), u.end(), 0.0);
        u[dof] = 1.0;

        for (uint n_bound = 0; n_bound < 3; n_bound++) {
            // just checking orthogonality of Legendre 1D poly
            edge_boundaries[n_bound].ComputeUgp(u, u_gp);

            double inner_product = edge_boundaries[n_bound].IntegrationLambda(dof, u_gp);

            if (!almost_equal(inner_product, 1.0 / (2 * dof + 1))) {
                error_found = true;
            }
        }
    }

    // generate edge intenrals
    std::vector<InterfaceType> interfaces;

    for (auto& rb : raw_boundary[SWE::BoundaryTypes::land]) {
        interfaces.emplace_back(InterfaceType(rb.second, rb.second));
    }

    std::vector<EdgeInternalType> edge_internals;

    for (auto& intface : interfaces) {
        edge_internals.emplace_back(EdgeInternalType(intface));
    }

    for (uint dof = 0; dof < 11; dof++) {
        std::fill(u.begin(), u.end(), 0.0);
        u[dof] = 1.0;

        for (uint n_int = 0; n_int < 3; n_int++) {
            // just checking orthogonality of Legendre 1D poly
            edge_internals[n_int].ComputeUgp(u, u_gp);

            double inner_product = edge_internals[n_int].IntegrationLambda(dof, u_gp);

            if (!almost_equal(inner_product, 1.0 / (2 * dof + 1))) {
                error_found = true;
            }
        }
    }

    if (error_found) {
        return 1;
    }

    return 0;
}