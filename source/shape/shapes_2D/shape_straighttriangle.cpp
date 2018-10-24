#include "../shapes_2D.hpp"

namespace Shape {
StraightTriangle::StraightTriangle(std::vector<Point<3>>&& nodal_coordinates) : Shape<2>(std::move(nodal_coordinates)) {
    // check if element nodes are ccw, swap if necessary
    // note that point selection doesn't matter since the Jacobian is
    // constant
    if (this->GetJdet(std::vector<Point<2>>(0))[0] < 0) {
        std::swap(this->nodal_coordinates[0], this->nodal_coordinates[2]);
    }
}

std::vector<uint> StraightTriangle::GetBoundaryNodeID(const uint bound_id, const std::vector<uint> node_ID) {
    std::vector<uint> bound_node_ID(2);

    bound_node_ID[0] = node_ID[(bound_id + 1) % 3];
    bound_node_ID[1] = node_ID[(bound_id + 2) % 3];

    return bound_node_ID;
}

Point<2> StraightTriangle::GetBarycentricCoordinates() {
    Point<2> baryctr_coord;

    baryctr_coord[GlobalCoord::x] =
        (this->nodal_coordinates[0][GlobalCoord::x] + this->nodal_coordinates[1][GlobalCoord::x] +
         this->nodal_coordinates[2][GlobalCoord::x]) /
        3.0;

    baryctr_coord[GlobalCoord::y] =
        (this->nodal_coordinates[0][GlobalCoord::y] + this->nodal_coordinates[1][GlobalCoord::y] +
         this->nodal_coordinates[2][GlobalCoord::y]) /
        3.0;

    return baryctr_coord;
}

std::vector<Point<2>> StraightTriangle::GetMidpointCoordinates() {
    std::vector<Point<2>> midpoint_coord(3);

    for (uint midpt = 0; midpt < 3; ++midpt) {
        midpoint_coord[midpt][GlobalCoord::x] = (this->nodal_coordinates[(midpt + 1) % 3][GlobalCoord::x] +
                                                 this->nodal_coordinates[(midpt + 2) % 3][GlobalCoord::x]) /
                                                2.0;

        midpoint_coord[midpt][GlobalCoord::y] = (this->nodal_coordinates[(midpt + 1) % 3][GlobalCoord::y] +
                                                 this->nodal_coordinates[(midpt + 2) % 3][GlobalCoord::y]) /
                                                2.0;
    }

    return midpoint_coord;
}

DynVector<double> StraightTriangle::GetJdet(const std::vector<Point<2>>& points) {
    DynVector<double> J_det(1);

    StatMatrix<double, 2, 2> J;

    J(0, 0) = (this->nodal_coordinates[1][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0;
    J(0, 1) = (this->nodal_coordinates[2][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0;
    J(1, 0) = (this->nodal_coordinates[1][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0;
    J(1, 1) = (this->nodal_coordinates[2][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0;

    J_det[0] = determinant(J);

    return J_det;
}

DynVector<double> StraightTriangle::GetSurfaceJ(const uint bound_id, const std::vector<Point<2>>& points) {
    DynVector<double> surface_J(1);

    surface_J[0] = sqrt(pow(this->nodal_coordinates[(bound_id + 2) % 3][GlobalCoord::x] -
                                this->nodal_coordinates[(bound_id + 1) % 3][GlobalCoord::x],
                            2.0) +
                        pow(this->nodal_coordinates[(bound_id + 2) % 3][GlobalCoord::y] -
                                this->nodal_coordinates[(bound_id + 1) % 3][GlobalCoord::y],
                            2.0)) /
                   2.0;  // half length for straight edge

    return surface_J;
}

AlignedVector<StatMatrix<double, 2, 2>>  StraightTriangle::GetJinv(const std::vector<Point<2>>& points) {
    AlignedVector<StatMatrix<double, 2, 2>> J_inv(1);

    StatMatrix<double, 2, 2> J;

    J(0, 0) = (this->nodal_coordinates[1][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0;
    J(0, 1) = (this->nodal_coordinates[2][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0;
    J(1, 0) = (this->nodal_coordinates[1][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0;
    J(1, 1) = (this->nodal_coordinates[2][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0;

    J_inv[0] = inverse(J);

    return J_inv;
}

AlignedVector<StatVector<double, 2>> StraightTriangle::GetSurfaceNormal(const uint bound_id,
									const std::vector<Point<2>>& points) {
    AlignedVector<StatVector<double, 2>> surface_normal(1);

    StatMatrix<double, 2, 2> J;

    J(0, 0) = (this->nodal_coordinates[1][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0;
    J(0, 1) = (this->nodal_coordinates[2][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0;
    J(1, 0) = (this->nodal_coordinates[1][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0;
    J(1, 1) = (this->nodal_coordinates[2][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0;

    double det_J = determinant(J);
    double cw    = det_J / std::abs(det_J);  // CW or CCW

    double length = std::hypot(this->nodal_coordinates[(bound_id + 2) % 3][GlobalCoord::x] -
			       this->nodal_coordinates[(bound_id + 1) % 3][GlobalCoord::x],
			       this->nodal_coordinates[(bound_id + 2) % 3][GlobalCoord::y] -
			       this->nodal_coordinates[(bound_id + 1) % 3][GlobalCoord::y]);

    surface_normal[0][GlobalCoord::x] = cw *
                                        (this->nodal_coordinates[(bound_id + 2) % 3][GlobalCoord::y] -
                                         this->nodal_coordinates[(bound_id + 1) % 3][GlobalCoord::y]) /
                                        length;
    surface_normal[0][GlobalCoord::y] = -cw *
                                        (this->nodal_coordinates[(bound_id + 2) % 3][GlobalCoord::x] -
                                         this->nodal_coordinates[(bound_id + 1) % 3][GlobalCoord::x]) /
                                        length;

    return surface_normal;
}

DynMatrix<double> StraightTriangle::GetPsi(const std::vector<Point<2>>& points) {
    uint ndof = 3;
    uint npt  = points.size();

    DynMatrix<double> psi(ndof, npt);

    for (uint pt = 0; pt < npt; ++pt) {
        psi(0, pt) = -(points[pt][LocalCoordTri::z1] + points[pt][LocalCoordTri::z2]) / 2;  // N1
        psi(1, pt) = (1 + points[pt][LocalCoordTri::z1]) / 2;                               // N2
        psi(2, pt) = (1 + points[pt][LocalCoordTri::z2]) / 2;                               // N3
    }

    return psi;
}

std::array<DynMatrix<double>, 2> StraightTriangle::GetDPsi(const std::vector<Point<2>>& points) {
    uint ndof = 3;
    uint npt  = points.size();

    std::array<DynMatrix<double>, 2> dpsi;

    DynMatrix<double> dpsi_dx(ndof, npt);
    DynMatrix<double> dpsi_dy(ndof, npt);

    StatMatrix<double, 2, 2> J_inv = this->GetJinv(points)[0];

    for (uint pt = 0; pt < points.size(); ++pt) {
        dpsi_dx(0, pt) =
            -0.5 * J_inv(LocalCoordTri::z1, GlobalCoord::x) - 0.5 * J_inv(LocalCoordTri::z2, GlobalCoord::x);
        dpsi_dy(0, pt) =
            -0.5 * J_inv(LocalCoordTri::z1, GlobalCoord::y) - 0.5 * J_inv(LocalCoordTri::z2, GlobalCoord::y);

        dpsi_dx(1, pt) = 0.5 * J_inv(LocalCoordTri::z1, GlobalCoord::x);
        dpsi_dy(1, pt) = 0.5 * J_inv(LocalCoordTri::z1, GlobalCoord::y);

        dpsi_dx(2, pt) = 0.5 * J_inv(LocalCoordTri::z2, GlobalCoord::x);
        dpsi_dy(2, pt) = 0.5 * J_inv(LocalCoordTri::z2, GlobalCoord::y);
    }

    dpsi[GlobalCoord::x] = dpsi_dx;
    dpsi[GlobalCoord::y] = dpsi_dy;

    return dpsi;
}

DynMatrix<double> StraightTriangle::GetBoundaryPsi(const uint bound_id, const std::vector<Point<1>>& points) {
    uint ndof = 2;
    uint npt  = points.size();

    DynMatrix<double> psi_bound(ndof, npt);

    for (uint pt = 0; pt < npt; ++pt) {
        psi_bound(0, pt) = (1 - points[pt][LocalCoordTri::z1]) / 2;  // N1
        psi_bound(1, pt) = (1 + points[pt][LocalCoordTri::z1]) / 2;  // N2
    }

    return psi_bound;
}

std::vector<Point<2>> StraightTriangle::LocalToGlobalCoordinates(const std::vector<Point<2>>& points) {
    uint npt = points.size();

    std::vector<Point<2>> global_coordinates(npt);

    DynMatrix<double> psi_pts = this->GetPsi(points);

    for (uint pt = 0; pt < npt; ++pt) {
        global_coordinates[pt][GlobalCoord::x] = this->nodal_coordinates[0][GlobalCoord::x] * psi_pts(0, pt) +
                                                 this->nodal_coordinates[1][GlobalCoord::x] * psi_pts(1, pt) +
                                                 this->nodal_coordinates[2][GlobalCoord::x] * psi_pts(2, pt);

        global_coordinates[pt][GlobalCoord::y] = this->nodal_coordinates[0][GlobalCoord::y] * psi_pts(0, pt) +
                                                 this->nodal_coordinates[1][GlobalCoord::y] * psi_pts(1, pt) +
                                                 this->nodal_coordinates[2][GlobalCoord::y] * psi_pts(2, pt);
    }

    return global_coordinates;
}

void StraightTriangle::GetVTK(std::vector<Point<3>>& points, Array2D<uint>& cells) {
    uint number_pt = points.size();

    double z1;
    double z2;
    double dz = 2.0 / N_DIV;

    for (uint i = 0; i <= N_DIV; ++i) {
        for (uint j = 0; j <= N_DIV - i; ++j) {
            points.push_back({0, 0, 0});

            z1 = -1.0 + dz * j;
            z2 = -1.0 + dz * i;

            points.back()[0] = this->nodal_coordinates[0][GlobalCoord::x] * (-(z1 + z2) / 2.0) +
                               this->nodal_coordinates[1][GlobalCoord::x] * ((1 + z1) / 2.0) +
                               this->nodal_coordinates[2][GlobalCoord::x] * ((1 + z2) / 2.0);

            points.back()[1] = this->nodal_coordinates[0][GlobalCoord::y] * (-(z1 + z2) / 2.0) +
                               this->nodal_coordinates[1][GlobalCoord::y] * ((1 + z1) / 2.0) +
                               this->nodal_coordinates[2][GlobalCoord::y] * ((1 + z2) / 2.0);

            points.back()[2] = 0;
        }
    }

    uint pt_ID;

    for (uint i = 0; i < N_DIV; ++i) {
        for (uint j = 0; j < N_DIV - i; ++j) {
            cells.push_back(std::vector<uint>(4));

            pt_ID = number_pt + (N_DIV + 1) * (N_DIV + 2) / 2 - (N_DIV - i + 1) * (N_DIV - i + 2) / 2 + j;

            cells.back()[0] = VTKElementTypes::straight_triangle;
            cells.back()[1] = pt_ID;
            cells.back()[2] = pt_ID + 1;
            cells.back()[3] = pt_ID + (N_DIV + 1 - i);
        }
    }

    for (uint i = 1; i < N_DIV; ++i) {
        for (uint j = 0; j < N_DIV - i; ++j) {
            cells.push_back(std::vector<uint>(4));

            pt_ID = number_pt + (N_DIV + 1) * (N_DIV + 2) / 2 - (N_DIV - i + 1) * (N_DIV - i + 2) / 2 + j;

            cells.back()[0] = VTKElementTypes::straight_triangle;
            cells.back()[1] = pt_ID;
            cells.back()[2] = pt_ID + 1;
            cells.back()[3] = pt_ID - (N_DIV + 2 - i) + 1;
        }
    }
}
}
