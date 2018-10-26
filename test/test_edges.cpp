#include "general_definitions.hpp"
#include "utilities/almost_equal.hpp"
#include "geometry/mesh_definitions.hpp"
#include "preprocessor/input_parameters.hpp"
#include "problem/SWE/problem_function_files/swe_true_solution_functions.hpp"
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
        Geometry::Interface<1, Integration::GaussLegendre_1D, SWE::EHDG::Data, SWE::EHDG::ISP::Internal>;

    using EdgeBoundaryTypes =
        Geometry::EdgeBoundaryTypeTuple<SWE::EHDG::EdgeData,
                                        Geometry::BoundaryTypeTuple<SWE::EHDG::Data, SWE::EHDG::BC::Land>>::Type;
    using EdgeBoundaryType = typename std::tuple_element<0, EdgeBoundaryTypes>::type;

    using EdgeInterfaceTypes =
        Geometry::EdgeInterfaceTypeTuple<SWE::EHDG::EdgeData,
                                         Geometry::InterfaceTypeTuple<SWE::EHDG::Data, SWE::EHDG::ISP::Internal>>::Type;
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

    DynMatrix<double> u(1, 11);                           // ndof for p = 10
    DynMatrix<double> u_gp(1, integration.GetNumGP(21));  // ngp for 2*p+1

    DynMatrix<double> u_proj(1, 11);

    u_proj(0, 0)  = 1.0;
    u_proj(0, 1)  = 2.0;
    u_proj(0, 2)  = 3.0;
    u_proj(0, 3)  = 4.0;
    u_proj(0, 4)  = 5.0;
    u_proj(0, 5)  = 6.0;
    u_proj(0, 6)  = 7.0;
    u_proj(0, 7)  = 8.0;
    u_proj(0, 8)  = 9.0;
    u_proj(0, 9)  = 10.0;
    u_proj(0, 10) = 11.0;

    DynMatrix<double> u_proj_res(1, 11);
    DynMatrix<double> u_proj_gp(1, integration.GetNumGP(21));  // ngp for 2*p+1

    // generate edge boundaries
    AlignedVector<BoundaryType> boundaries;

    for (auto& rb : raw_boundary[SWE::BoundaryTypes::land]) {
        boundaries.emplace_back(BoundaryType(std::move(rb.second)));
    }

    AlignedVector<EdgeBoundaryType> edge_boundaries;

    for (auto& bound : boundaries) {
        edge_boundaries.emplace_back(EdgeBoundaryType(bound));
    }

    for (uint n_bound = 0; n_bound < 3; ++n_bound) {
        u_proj_gp  = edge_boundaries[n_bound].ComputeUgp(u_proj);
        u_proj_res = edge_boundaries[n_bound].L2Projection(u_proj_gp);

        for (uint dof = 0; dof < 11; ++dof) {
            if (!almost_equal(u_proj(0, dof), u_proj_res(0, dof))) {
                error_found = true;

                std::cout << "Error found edge bound L2Projection" << std::endl;
            }
        }
    }

    for (uint dof = 0; dof < 11; ++dof) {
        set_constant(u, 0.0);
        row(u, 0)[dof] = 1.0;

        for (uint n_bound = 0; n_bound < 3; ++n_bound) {
            // just checking orthogonality of Legendre 1D poly
            u_gp = edge_boundaries[n_bound].ComputeUgp(u);

            double inner_product = edge_boundaries[n_bound].IntegrationLambda(dof, u_gp)[0];

            if (!almost_equal(inner_product, 1.0 / (2 * dof + 1))) {
                error_found = true;

                std::cout << "Error found edge bound int lambda" << std::endl;
            }

            set_constant(u_gp, 1.0);
            ;

            inner_product = edge_boundaries[n_bound].IntegrationLambdaLambda(dof, dof, u_gp)[0];

            if (!almost_equal(inner_product, 1.0 / (2 * dof + 1))) {
                error_found = true;

                std::cout << "Error found edge bound int lambda lambda" << std::endl;
            }
        }
    }

    // generate edge intenrals
    AlignedVector<InterfaceType> interfaces;

    for (auto& rb : raw_boundary[SWE::BoundaryTypes::land]) {
        interfaces.emplace_back(InterfaceType(std::move(rb.second), std::move(rb.second)));
    }

    AlignedVector<EdgeInterfaceType> edge_interfaces;

    for (auto& intface : interfaces) {
        edge_interfaces.emplace_back(EdgeInterfaceType(intface));
    }

    for (uint n_int = 0; n_int < 3; n_int++) {
        u_proj_gp  = edge_interfaces[n_int].ComputeUgp(u_proj);
        u_proj_res = edge_interfaces[n_int].L2Projection(u_proj_gp);

        for (uint dof = 0; dof < 11; ++dof) {
            if (!almost_equal(u_proj(0, dof), u_proj_res(0, dof))) {
                error_found = true;

                std::cout << "Error found edge int L2Projection" << std::endl;
            }
        }
    }

    for (uint dof = 0; dof < 11; ++dof) {
        set_constant(u, 0.0);
        row(u, 0)[dof] = 1.0;

        for (uint n_int = 0; n_int < 3; n_int++) {
            // just checking orthogonality of Legendre 1D poly
            u_gp = edge_interfaces[n_int].ComputeUgp(u);

            double inner_product = edge_interfaces[n_int].IntegrationLambda(dof, u_gp)[0];

            if (!almost_equal(inner_product, 1.0 / (2 * dof + 1))) {
                error_found = true;

                std::cout << "Error found edge int int lambda" << std::endl;
            }

            set_constant(u_gp, 1.0);
            ;

            inner_product = edge_interfaces[n_int].IntegrationLambdaLambda(dof, dof, u_gp)[0];

            if (!almost_equal(inner_product, 1.0 / (2 * dof + 1))) {
                error_found = true;

                std::cout << "Error found edge int int lambda lambda" << std::endl;
            }
        }
    }

    DynMatrix<double> u_phi(1, 66);                           // ndof for p = 10
    DynMatrix<double> u_phi_gp(1, integration.GetNumGP(21));  // ngp for 2*p+1
    DynMatrix<double> unit(1, integration.GetNumGP(21));
    set_constant(unit, 1.0);

    for (uint n_bound = 0; n_bound < 1; ++n_bound) {
        for (uint dof = 0; dof < 66; ++dof) {
            set_constant(u_phi, 0.0);
            row(u_phi, 0)[dof] = 1.0;

            u_phi_gp = edge_boundaries[n_bound].boundary.ComputeUgp(u_phi);

            for (uint doff = 0; doff < 11; ++doff) {
                if (!almost_equal(edge_boundaries[n_bound].IntegrationPhiLambda(dof, doff, unit)[0],
                                  edge_boundaries[n_bound].IntegrationLambda(doff, u_phi_gp)[0],
                                  1.e+03)) {
                    error_found = true;

                    std::cout << std::setprecision(16)
                              << edge_boundaries[n_bound].IntegrationPhiLambda(dof, doff, unit)[0] << ' '
                              << edge_boundaries[n_bound].IntegrationLambda(doff, u_phi_gp)[0] << std::endl;

                    std::cerr << "Error found in boundary edge in IntegrationPhiLambda\n"
                              << "  dof: " << dof << " doff: " << doff << std::endl;
                }
            }
        }
    }

    for (uint n_int = 0; n_int < 3; n_int++) {
        for (uint dof = 0; dof < 66; ++dof) {
            set_constant(u_phi, 0.0);
            row(u_phi, 0)[dof] = 1.0;

            u_phi_gp = edge_interfaces[n_int].interface.ComputeUgpIN(u_phi);

            for (uint doff = 0; doff < 11; ++doff) {
                if (!almost_equal(edge_interfaces[n_int].IntegrationPhiLambdaIN(dof, doff, unit)[0],
                                  edge_interfaces[n_int].IntegrationLambda(doff, u_phi_gp)[0],
                                  1.e+03)) {
                    error_found = true;

                    std::cerr << "Error found in interface edge in IntegrationPhiLambdaIN" << std::endl;
                }
            }

            for (uint doff = 0; doff < 11; ++doff) {
                if (!almost_equal(std::abs(edge_interfaces[n_int].IntegrationPhiLambdaIN(dof, doff, unit)[0]),
                                  std::abs(edge_interfaces[n_int].IntegrationPhiLambdaEX(dof, doff, unit)[0]),
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