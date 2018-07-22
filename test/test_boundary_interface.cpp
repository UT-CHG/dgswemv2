#include "general_definitions.hpp"
#include "utilities/almost_equal.hpp"
#include "geometry/mesh_definitions.hpp"
#include "problem/SWE/discretization_RKDG/rkdg_swe_problem.hpp"

int main() {
    using Utilities::almost_equal;
    bool error_found = false;

    using MasterType  = Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>;
    using ShapeType   = Shape::StraightTriangle;
    using ElementType = Geometry::Element<2, MasterType, ShapeType, SWE::RKDG::Data>;

    using RawBoundaryType = Geometry::RawBoundary<1, SWE::RKDG::Data>;
    using BoundaryType    = Geometry::Boundary<1, Integration::GaussLegendre_1D, SWE::RKDG::Data, SWE::RKDG::BC::Land>;
    using InterfaceType =
        Geometry::Interface<1, Integration::GaussLegendre_1D, SWE::RKDG::Data, SWE::RKDG::IS::Internal>;

    // make an equilateral triangle
    DynVector<Point<3>> vrtxs(3);
    vrtxs[0] = {-0.5, 0., 0.};
    vrtxs[1] = {0.5, 0., 0.};
    vrtxs[2] = {0, std::sqrt(3.) / 2., 0.};

    MasterType master(10);
    ShapeType shape(std::move(vrtxs));

    vrtxs.resize(3);
    vrtxs[0] = {-0.5, 0., 0.};
    vrtxs[1] = {0.5, 0., 0.};
    vrtxs[2] = {0, std::sqrt(3.) / 2., 0.};

    ElementType triangle(0,
                         master,
                         std::move(vrtxs),
                         std::move(DynVector<uint>(3, 0)),
                         std::move(DynVector<uint>{DEFAULT_ID, DEFAULT_ID, DEFAULT_ID}),
                         std::move(DynVector<unsigned char>{
                             SWE::BoundaryTypes::land, SWE::BoundaryTypes::land, SWE::BoundaryTypes::land}));

    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>> raw_boundary;

    // generate boundaries
    triangle.CreateRawBoundaries(raw_boundary);

    std::vector<BoundaryType> boundaries;

    for (auto& rb : raw_boundary[SWE::BoundaryTypes::land]) {
        boundaries.emplace_back(BoundaryType(std::move(rb.second)));
    }

    // Check Integrations
    Integration::Dunavant_2D integ_2D;
    DynVector<Point<2>> gp_2D = integ_2D.GetRule(20).second;

    DynMatrix<double> x_node(1, 3);
    DynMatrix<double> y_node(1, 3);

    x_node(0, 0) = -0.5;
    x_node(0, 1) = 0.5;
    x_node(0, 2) = 0.0;

    y_node(0, 0) = 0.0;
    y_node(0, 1) = 0.0;
    y_node(0, 2) = std::sqrt(3.0) / 2.0;

    DynMatrix<double> x(1, triangle.data.get_ngp_internal());
    DynMatrix<double> y(1, triangle.data.get_ngp_internal());

    x = triangle.ComputeNodalUgp(x_node);
    y = triangle.ComputeNodalUgp(y_node);

    DynMatrix<double> F_vals_int(2, triangle.data.get_ngp_internal());
    DynMatrix<double> divF_vals_int(1, triangle.data.get_ngp_internal());

    for (uint gp = 0; gp < triangle.data.get_ngp_internal(); gp++) {
        F_vals_int(GlobalCoord::x, gp) = std::pow(x(0, gp), 2);
        F_vals_int(GlobalCoord::y, gp) = std::pow(y(0, gp), 2);

        divF_vals_int(0, gp) = 2 * x(0, gp) + 2 * y(0, gp);
    }

    Integration::GaussLegendre_1D integ_1D;
    DynVector<Point<1>> gp_1D = integ_1D.GetRule(21).second;
    DynVector<Point<2>> gp_bound;

    DynMatrix<double> F_vals_bound(2, triangle.data.get_ngp_boundary(0));
    DynMatrix<double> Fn_vals_bound(3, triangle.data.get_ngp_boundary(0));

    for (uint n_bound = 0; n_bound < 3; n_bound++) {
        gp_bound = master.BoundaryToMasterCoordinates(n_bound, gp_1D);

        x.resize(1, triangle.data.get_ngp_boundary(0));
        y.resize(1, triangle.data.get_ngp_boundary(0));

        x = boundaries[n_bound].ComputeNodalUgp(x_node);
        y = boundaries[n_bound].ComputeNodalUgp(y_node);

        for (uint gp = 0; gp < triangle.data.get_ngp_boundary(0); gp++) {
            F_vals_bound(GlobalCoord::x, gp) = std::pow(x(0, gp), 2);
            F_vals_bound(GlobalCoord::y, gp) = std::pow(y(0, gp), 2);

            Fn_vals_bound(n_bound, gp) =
                F_vals_bound(GlobalCoord::x, gp) * boundaries[n_bound].surface_normal[gp][GlobalCoord::x] +
                F_vals_bound(GlobalCoord::y, gp) * boundaries[n_bound].surface_normal[gp][GlobalCoord::y];
        }
    }

    double divergence;

    for (uint dof = 0; dof < 66; dof++) {
        divergence = triangle.IntegrationPhi(dof, row(divF_vals_int, 0)) +
                     triangle.IntegrationDPhi(GlobalCoord::x, dof, row(F_vals_int, GlobalCoord::x)) +
                     triangle.IntegrationDPhi(GlobalCoord::y, dof, row(F_vals_int, GlobalCoord::y)) -
                     boundaries[0].IntegrationPhi(dof, row(Fn_vals_bound, 0)) -
                     boundaries[1].IntegrationPhi(dof, row(Fn_vals_bound, 1)) -
                     boundaries[2].IntegrationPhi(dof, row(Fn_vals_bound, 2));

        if (!almost_equal(divergence, 0.0, 1.e+03)) {
            error_found = true;

            std::cerr << "Error found in boundary in IntegrationPhi" << std::endl;
        }
    }

    DynMatrix<double> mod_vals(1, triangle.data.get_ndof());
    DynMatrix<double> gp_vals(1, triangle.data.get_ngp_boundary(0));
    DynMatrix<double> unit(1, triangle.data.get_ngp_boundary(0), 1.0);

    // Check ComputeUgp
    for (uint dof = 0; dof < 66; dof++) {
        std::fill(row(mod_vals, 0).begin(), row(mod_vals, 0).end(), 0.0);
        row(mod_vals, 0)[dof] = 1.0;

        for (uint n_bound = 0; n_bound < 3; n_bound++) {
            gp_vals = boundaries[n_bound].ComputeUgp(mod_vals);
            if (!almost_equal(boundaries[n_bound].IntegrationPhi(dof, unit)[0],
                              boundaries[n_bound].IntegrationPhi(0, gp_vals)[0],
                              1.e+03)) {
                error_found = true;

                std::cerr << "Error found in boundary in ComputeUgp" << std::endl;
            }
        }
    }

    // Check ComputeNodalUgp, ComputeBoundaryNodalUgp and Integration
    DynMatrix<double> nodal_vals(1, 3);
    nodal_vals(0, 0) = 1.0;
    nodal_vals(0, 1) = 2.0;
    nodal_vals(0, 2) = 3.0;

    DynMatrix<double> bound_nodal_vals(1, 2);
    DynMatrix<double> nodal_vals_gp(1, gp_1D.size());
    DynMatrix<double> nodal_vals_gp_comp(1, gp_1D.size());

    for (uint n_bound = 0; n_bound < 3; n_bound++) {
        gp_bound = master.BoundaryToMasterCoordinates(n_bound, gp_1D);

        DynMatrix<double> psi_gp = shape.GetPsi(gp_bound);

        for (uint gp = 0; gp < gp_bound.size(); gp++) {
            nodal_vals_gp(0, gp) =
                psi_gp(0, gp) * nodal_vals(0, 0) + psi_gp(1, gp) * nodal_vals(0, 1) + psi_gp(2, gp) * nodal_vals(0, 2);
        }

        nodal_vals_gp_comp = boundaries[n_bound].ComputeNodalUgp(nodal_vals);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            if (!almost_equal(nodal_vals_gp(0, gp), nodal_vals_gp_comp(0, gp))) {
                error_found = true;

                std::cerr << "Error found in boundary in ComputeNodalUgp" << std::endl;
            }
        }

        bound_nodal_vals(0, 0) = nodal_vals(0, (n_bound + 1) % 3);
        bound_nodal_vals(0, 1) = nodal_vals(0, (n_bound + 2) % 3);

        DynMatrix<double> psi_bound_gp = shape.GetBoundaryPsi(n_bound, gp_1D);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            nodal_vals_gp(0, gp) =
                psi_bound_gp(0, gp) * bound_nodal_vals(0, 0) + psi_bound_gp(1, gp) * bound_nodal_vals(0, 1);
        }

        nodal_vals_gp_comp = boundaries[n_bound].ComputeBoundaryNodalUgp(bound_nodal_vals);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            if (!almost_equal(nodal_vals_gp(0, gp), nodal_vals_gp_comp(0, gp))) {
                error_found = true;

                std::cerr << "Error found in boundary in ComputeBoundaryNodalUgp" << std::endl;
            }
        }

        double exact_integration = (nodal_vals(0, (n_bound + 1) % 3) + nodal_vals(0, (n_bound + 2) % 3)) / 2.0;

        if (!almost_equal(exact_integration, boundaries[n_bound].Integration(nodal_vals_gp_comp)[0])) {
            error_found = true;

            std::cerr << "Error found in boundary in Integration" << std::endl;
        }
    }

    // Check integrations PhiPhi
    for (uint n_bound = 0; n_bound < 3; n_bound++) {
        for (uint dof = 0; dof < 66; dof++) {
            std::fill(row(mod_vals, 0).begin(), row(mod_vals, 0).end(), 0.0);
            row(mod_vals, 0)[dof] = 1.0;

            gp_vals = boundaries[n_bound].ComputeUgp(mod_vals);

            for (uint doff = 0; doff < 66; doff++) {
                if (!almost_equal(boundaries[n_bound].IntegrationPhi(doff, gp_vals)[0],
                                  boundaries[n_bound].IntegrationPhiPhi(dof, doff, unit)[0],
                                  1.e+03)) {
                    error_found = true;

                    std::cerr << "Error found in boundary in IntegrationPhiPhi" << std::endl;
                }
            }
        }
    }

    // generate Interfaces
    std::vector<InterfaceType> interfaces;

    for (auto& rb : raw_boundary[SWE::BoundaryTypes::land]) {
        interfaces.emplace_back(InterfaceType(std::move(rb.second), std::move(rb.second)));
    }

    // Check IntegrationPhiIN/EX
    for (uint dof = 0; dof < 66; dof++) {
        std::fill(row(mod_vals, 0).begin(), row(mod_vals, 0).end(), 0.0);
        row(mod_vals, 0)[dof] = 1.0;

        for (uint n_intface = 0; n_intface < 3; n_intface++) {
            if (!almost_equal(boundaries[n_intface].IntegrationPhi(dof, row(Fn_vals_bound, n_intface)),
                              interfaces[n_intface].IntegrationPhiIN(dof, row(Fn_vals_bound, n_intface)))) {
                error_found = true;

                std::cerr << "Error found in interface in IntegrationPhiIN" << std::endl;
            }

            if (!almost_equal(boundaries[n_intface].IntegrationPhi(dof, row(Fn_vals_bound, n_intface)),
                              interfaces[n_intface].IntegrationPhiEX(dof, row(Fn_vals_bound, n_intface)))) {
                error_found = true;

                std::cerr << "Error found in interface in IntegrationPhiEX" << std::endl;
            }
        }
    }

    // Check ComputeUgpIN/EX
    for (uint dof = 0; dof < 66; dof++) {
        std::fill(row(mod_vals, 0).begin(), row(mod_vals, 0).end(), 0.0);
        row(mod_vals, 0)[dof] = 1.0;

        for (uint n_intface = 0; n_intface < 3; n_intface++) {
            gp_vals = interfaces[n_intface].ComputeUgpIN(mod_vals);
            if (!almost_equal(interfaces[n_intface].IntegrationPhiIN(dof, unit)[0],
                              interfaces[n_intface].IntegrationPhiIN(0, gp_vals)[0],
                              1.e+03)) {
                error_found = true;

                std::cerr << "Error found in interface in ComputeUgpIN" << std::endl;
            }

            gp_vals = interfaces[n_intface].ComputeUgpEX(mod_vals);
            if (!almost_equal(interfaces[n_intface].IntegrationPhiEX(dof, unit)[0],
                              interfaces[n_intface].IntegrationPhiEX(0, gp_vals)[0],
                              1.e+03)) {
                error_found = true;

                std::cerr << "Error found in interface in ComputeUgpEX" << std::endl;
            }
        }
    }

    // Check ComputeNodalUgp, ComputeBoundaryNodalUgp and Integration
    for (uint n_intface = 0; n_intface < 3; n_intface++) {
        gp_bound = master.BoundaryToMasterCoordinates(n_intface, gp_1D);

        DynMatrix<double> psi_gp = shape.GetPsi(gp_bound);

        for (uint gp = 0; gp < gp_bound.size(); gp++) {
            nodal_vals_gp(0, gp) =
                psi_gp(0, gp) * nodal_vals(0, 0) + psi_gp(1, gp) * nodal_vals(0, 1) + psi_gp(2, gp) * nodal_vals(0, 2);
        }

        nodal_vals_gp_comp = interfaces[n_intface].ComputeNodalUgpIN(nodal_vals);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            if (!almost_equal(nodal_vals_gp(0, gp), nodal_vals_gp_comp(0, gp))) {
                error_found = true;

                std::cerr << "Error found in boundary in ComputeNodalUgpIN" << std::endl;
            }
        }

        bound_nodal_vals(0, 0) = nodal_vals(0, (n_intface + 1) % 3);
        bound_nodal_vals(0, 1) = nodal_vals(0, (n_intface + 2) % 3);

        DynMatrix<double> psi_bound_gp = shape.GetBoundaryPsi(n_intface, gp_1D);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            nodal_vals_gp(0, gp) =
                psi_bound_gp(0, gp) * bound_nodal_vals(0, 0) + psi_bound_gp(1, gp) * bound_nodal_vals(0, 1);
        }

        nodal_vals_gp_comp = interfaces[n_intface].ComputeBoundaryNodalUgpIN(bound_nodal_vals);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            if (!almost_equal(nodal_vals_gp(0, gp), nodal_vals_gp_comp(0, gp))) {
                error_found = true;

                std::cerr << "Error found in boundary in ComputeBoundaryNodalUgpIN" << std::endl;
            }
        }

        double exact_integration = (nodal_vals(0, (n_intface + 1) % 3) + nodal_vals(0, (n_intface + 2) % 3)) / 2.0;

        if (!almost_equal(exact_integration, interfaces[n_intface].IntegrationIN(nodal_vals_gp_comp)[0])) {
            error_found = true;

            std::cerr << "Error found in boundary in IntegrationIN" << std::endl;
        }

        nodal_vals_gp_comp = interfaces[n_intface].ComputeNodalUgpEX(nodal_vals);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            if (!almost_equal(nodal_vals_gp(0, gp), nodal_vals_gp_comp(0, gp))) {
                error_found = true;

                std::cerr << "Error found in boundary in ComputeNodalUgpEX" << std::endl;
            }
        }

        psi_bound_gp = shape.GetBoundaryPsi(n_intface, gp_1D);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            nodal_vals_gp(0, gp) =
                psi_bound_gp(0, gp) * bound_nodal_vals(0, 0) + psi_bound_gp(1, gp) * bound_nodal_vals(0, 1);
        }

        nodal_vals_gp_comp = interfaces[n_intface].ComputeBoundaryNodalUgpEX(bound_nodal_vals);

        for (uint gp = 0; gp < gp_1D.size(); gp++) {
            if (!almost_equal(nodal_vals_gp(0, gp), nodal_vals_gp_comp(0, gp))) {
                error_found = true;

                std::cerr << "Error found in boundary in ComputeBoundaryNodalUgpEX" << std::endl;
            }
        }

        if (!almost_equal(exact_integration, interfaces[n_intface].IntegrationEX(nodal_vals_gp_comp)[0])) {
            error_found = true;

            std::cerr << "Error found in boundary in IntegrationEX" << std::endl;
        }
    }

    // Check integrations PhiPhi
    for (uint n_intface = 0; n_intface < 3; n_intface++) {
        for (uint dof = 0; dof < 66; dof++) {
            std::fill(row(mod_vals, 0).begin(), row(mod_vals, 0).end(), 0.0);
            row(mod_vals, 0)[dof] = 1.0;

            gp_vals = interfaces[n_intface].ComputeUgpIN(mod_vals);

            for (uint doff = 0; doff < 66; doff++) {
                if (!almost_equal(interfaces[n_intface].IntegrationPhiIN(doff, gp_vals)[0],
                                  interfaces[n_intface].IntegrationPhiPhiIN(dof, doff, unit)[0],
                                  1.e+03)) {
                    error_found = true;

                    std::cerr << "Error found in interface in IntegrationPhiPhiIN" << std::endl;
                }
            }

            gp_vals = interfaces[n_intface].ComputeUgpEX(mod_vals);

            for (uint doff = 0; doff < 66; doff++) {
                if (!almost_equal(interfaces[n_intface].IntegrationPhiEX(doff, gp_vals)[0],
                                  interfaces[n_intface].IntegrationPhiPhiEX(dof, doff, unit)[0],
                                  1.e+03)) {
                    error_found = true;

                    std::cerr << "Error found in interface in IntegrationPhiPhiEX" << std::endl;
                }
            }
        }
    }

    if (error_found) {
        return 1;
    }

    return 0;
}