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

    using EdgeBoundaryTypes =
        Geometry::EdgeBoundaryTypeTuple<SWE::EHDG::EdgeData,
                                        Geometry::BoundaryTypeTuple<SWE::EHDG::Data, SWE::EHDG::BC::Land>>::Type;
    using EdgeBoundaryType = typename std::tuple_element<0, EdgeBoundaryTypes>::type;

    using EdgeInterfaceTypes =
        Geometry::EdgeInterfaceTypeTuple<SWE::EHDG::EdgeData,
                                         Geometry::InterfaceTypeTuple<SWE::EHDG::Data, SWE::EHDG::IS::Internal>>::Type;
    using EdgeInterfaceType = typename std::tuple_element<0, EdgeInterfaceTypes>::type;

    // make an equilateral triangle
    std::vector<Point<3>> vrtxs(3);
    vrtxs[0] = {-0.5, 0., 0.};
    vrtxs[1] = {0.5, 0., 0.};
    vrtxs[2] = {0, std::sqrt(3.) / 2., 0.};

    MasterType master(10);

    ElementType triangle(0,
                         master,
                         std::move(vrtxs),
                         std::move(std::vector<uint>{0, 0, 0}),
                         std::move(std::vector<uint>{DEFAULT_ID, DEFAULT_ID, DEFAULT_ID}),
                         std::move(std::vector<unsigned char>{
                             SWE::BoundaryTypes::land, SWE::BoundaryTypes::land, SWE::BoundaryTypes::land}));

    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>> raw_boundary;

    triangle.CreateRawBoundaries(raw_boundary);

    Integration::GaussLegendre_1D integration;

    std::vector<double> u(11);                           // ndof for p = 10
    std::vector<double> u_gp(integration.GetNumGP(21));  // ngp for 2*p+1

    std::vector<double> u_proj{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0};
    std::vector<double> u_proj_res(11);
    std::vector<double> u_proj_gp(integration.GetNumGP(21));  // ngp for 2*p+1

    // generate edge boundaries
    std::vector<BoundaryType> boundaries;

    for (auto& rb : raw_boundary[SWE::BoundaryTypes::land]) {
        boundaries.emplace_back(BoundaryType(std::move(rb.second)));
    }

    std::vector<EdgeBoundaryType> edge_boundaries;

    for (auto& bound : boundaries) {
        edge_boundaries.emplace_back(EdgeBoundaryType(bound));
    }

    for (uint n_bound = 0; n_bound < 3; n_bound++) {
        edge_boundaries[n_bound].ComputeUgp(u_proj, u_proj_gp);
        edge_boundaries[n_bound].L2Projection(u_proj_gp, u_proj_res);

        for (uint dof = 0; dof < 11; dof++) {
            if (!almost_equal(u_proj[dof], u_proj_res[dof])) {
                error_found = true;
            }
        }
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

            std::fill(u_gp.begin(), u_gp.end(), 1.0);

            inner_product = edge_boundaries[n_bound].IntegrationLambdaLambda(dof, dof, u_gp);

            if (!almost_equal(inner_product, 1.0 / (2 * dof + 1))) {
                error_found = true;
            }
        }
    }

    // generate edge intenrals
    std::vector<InterfaceType> interfaces;

    for (auto& rb : raw_boundary[SWE::BoundaryTypes::land]) {
        interfaces.emplace_back(InterfaceType(std::move(rb.second), std::move(rb.second)));
    }

    std::vector<EdgeInterfaceType> edge_interfaces;

    for (auto& intface : interfaces) {
        edge_interfaces.emplace_back(EdgeInterfaceType(intface));
    }

    for (uint n_int = 0; n_int < 3; n_int++) {
        edge_interfaces[n_int].ComputeUgp(u_proj, u_proj_gp);
        edge_interfaces[n_int].L2Projection(u_proj_gp, u_proj_res);

        for (uint dof = 0; dof < 11; dof++) {
            if (!almost_equal(u_proj[dof], u_proj_res[dof])) {
                error_found = true;
            }
        }
    }

    for (uint dof = 0; dof < 11; dof++) {
        std::fill(u.begin(), u.end(), 0.0);
        u[dof] = 1.0;

        for (uint n_int = 0; n_int < 3; n_int++) {
            // just checking orthogonality of Legendre 1D poly
            edge_interfaces[n_int].ComputeUgp(u, u_gp);

            double inner_product = edge_interfaces[n_int].IntegrationLambda(dof, u_gp);

            if (!almost_equal(inner_product, 1.0 / (2 * dof + 1))) {
                error_found = true;
            }

            std::fill(u_gp.begin(), u_gp.end(), 1.0);

            inner_product = edge_interfaces[n_int].IntegrationLambdaLambda(dof, dof, u_gp);

            if (!almost_equal(inner_product, 1.0 / (2 * dof + 1))) {
                error_found = true;
            }
        }
    }

    std::vector<double> u_phi(66);                           // ndof for p = 10
    std::vector<double> u_phi_gp(integration.GetNumGP(21));  // ngp for 2*p+1
    std::vector<double> unit(integration.GetNumGP(21), 1.0);

    for (uint n_bound = 0; n_bound < 1; n_bound++) {
        for (uint dof = 0; dof < 66; dof++) {
            std::fill(u_phi.begin(), u_phi.end(), 0.0);
            u_phi[dof] = 1.0;

            edge_boundaries[n_bound].boundary.ComputeUgp(u_phi, u_phi_gp);

            for (uint doff = 0; doff < 11; doff++) {
                if (!almost_equal(edge_boundaries[n_bound].IntegrationPhiLambda(dof, doff, unit),
                                  edge_boundaries[n_bound].IntegrationLambda(doff, u_phi_gp),
                                  1.e+03)) {
                    error_found = true;

                    std::cout << edge_boundaries[n_bound].IntegrationPhiLambda(dof, doff, unit) << ' '
                              << edge_boundaries[n_bound].IntegrationLambda(doff, u_phi_gp) << std::endl;

                    std::cerr << "Error found in boundary edge in IntegrationPhiLambda" << std::endl;
                }
            }
        }
    }

    for (uint n_int = 0; n_int < 3; n_int++) {
        for (uint dof = 0; dof < 66; dof++) {
            std::fill(u_phi.begin(), u_phi.end(), 0.0);
            u_phi[dof] = 1.0;

            edge_interfaces[n_int].interface.ComputeUgpIN(u_phi, u_phi_gp);

            for (uint doff = 0; doff < 11; doff++) {
                if (!almost_equal(edge_interfaces[n_int].IntegrationPhiLambdaIN(dof, doff, unit),
                                  edge_interfaces[n_int].IntegrationLambda(doff, u_phi_gp),
                                  1.e+03)) {
                    error_found = true;

                    std::cerr << "Error found in interface edge in IntegrationPhiLambdaIN" << std::endl;
                }
            }

            for (uint doff = 0; doff < 11; doff++) {
                if (!almost_equal(std::abs(edge_interfaces[n_int].IntegrationPhiLambdaIN(dof, doff, unit)),
                                  std::abs(edge_interfaces[n_int].IntegrationPhiLambdaEX(dof, doff, unit)),
                                  1.e+03)) {
                    error_found = true;

                    std::cerr << "Error found in interface edge in IntegrationPhiLambdaEX" << std::endl;
                }
            }
        }
    }

    if (error_found) {
        return 1;
    }

    return 0;
}