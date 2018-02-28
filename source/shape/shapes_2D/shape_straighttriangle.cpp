#include "../shapes_2D.hpp"

namespace Shape {
StraightTriangle::StraightTriangle(const std::vector<Point<2>>& nodal_coordinates) : Shape<2>(nodal_coordinates) {
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
    assert(this->nodal_coordinates.size() > 0);
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

    for (uint midpt = 0; midpt < 3; midpt++) {
        midpoint_coord[midpt][GlobalCoord::x] = (this->nodal_coordinates[(midpt + 1) % 3][GlobalCoord::x] +
                                                 this->nodal_coordinates[(midpt + 2) % 3][GlobalCoord::x]) /
                                                2.0;

        midpoint_coord[midpt][GlobalCoord::y] = (this->nodal_coordinates[(midpt + 1) % 3][GlobalCoord::y] +
                                                 this->nodal_coordinates[(midpt + 2) % 3][GlobalCoord::y]) /
                                                2.0;
    }

    return midpoint_coord;
}

std::vector<double> StraightTriangle::GetJdet(const std::vector<Point<2>>& points) {
    assert(this->nodal_coordinates.size() > 0);
    std::vector<double> J_det;

    Array2D<double> J(2);
    J[0].reserve(2);
    J[1].reserve(2);

    J[0].push_back((this->nodal_coordinates[1][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0);
    J[0].push_back((this->nodal_coordinates[2][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0);
    J[1].push_back((this->nodal_coordinates[1][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0);
    J[1].push_back((this->nodal_coordinates[2][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0);

    J_det.push_back(J[0][0] * J[1][1] - J[0][1] * J[1][0]);

    return J_det;
}

Array3D<double> StraightTriangle::GetJinv(const std::vector<Point<2>>& points) {
    assert(this->nodal_coordinates.size() > 0);
    Array3D<double> J_inv(2);
    J_inv[0].resize(2);
    J_inv[1].resize(2);

    Array2D<double> J(2);
    J[0].reserve(2);
    J[1].reserve(2);

    J[0].push_back((this->nodal_coordinates[1][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0);
    J[0].push_back((this->nodal_coordinates[2][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0);
    J[1].push_back((this->nodal_coordinates[1][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0);
    J[1].push_back((this->nodal_coordinates[2][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0);

    double det_J = (J[0][0] * J[1][1] - J[0][1] * J[1][0]);

    J_inv[0][0].push_back(J[1][1] / det_J);
    J_inv[0][1].push_back(-J[0][1] / det_J);
    J_inv[1][0].push_back(-J[1][0] / det_J);
    J_inv[1][1].push_back(J[0][0] / det_J);

    return J_inv;
}

std::vector<double> StraightTriangle::GetSurfaceJ(const uint bound_id, const std::vector<Point<2>>& points) {
    assert(this->nodal_coordinates.size() > 0);
    std::vector<double> surface_J;

    surface_J.push_back(sqrt(pow(this->nodal_coordinates[(bound_id + 2) % 3][GlobalCoord::x] -
                                     this->nodal_coordinates[(bound_id + 1) % 3][GlobalCoord::x],
                                 2.0) +
                             pow(this->nodal_coordinates[(bound_id + 2) % 3][GlobalCoord::y] -
                                     this->nodal_coordinates[(bound_id + 1) % 3][GlobalCoord::y],
                                 2.0)) /
                        2.0);  // half length for straight edge

    return surface_J;
}

Array2D<double> StraightTriangle::GetSurfaceNormal(const uint bound_id, const std::vector<Point<2>>& points) {
    assert(this->nodal_coordinates.size() > 0);
    Array2D<double> surface_normal(1);

    Array2D<double> J(2);
    J[0].reserve(2);
    J[1].reserve(2);

    J[0].push_back((this->nodal_coordinates[1][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0);
    J[0].push_back((this->nodal_coordinates[2][GlobalCoord::x] - this->nodal_coordinates[0][GlobalCoord::x]) / 2.0);
    J[1].push_back((this->nodal_coordinates[1][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0);
    J[1].push_back((this->nodal_coordinates[2][GlobalCoord::y] - this->nodal_coordinates[0][GlobalCoord::y]) / 2.0);

    double det_J = (J[0][0] * J[1][1] - J[0][1] * J[1][0]);
    double cw    = det_J / std::abs(det_J);  // CW or CCW

    double length = sqrt(pow(this->nodal_coordinates[(bound_id + 2) % 3][GlobalCoord::x] -
                                 this->nodal_coordinates[(bound_id + 1) % 3][GlobalCoord::x],
                             2.0) +
                         pow(this->nodal_coordinates[(bound_id + 2) % 3][GlobalCoord::y] -
                                 this->nodal_coordinates[(bound_id + 1) % 3][GlobalCoord::y],
                             2.0));

    surface_normal[0].push_back(cw *
                                (this->nodal_coordinates[(bound_id + 2) % 3][GlobalCoord::y] -
                                 this->nodal_coordinates[(bound_id + 1) % 3][GlobalCoord::y]) /
                                length);
    surface_normal[0].push_back(-cw *
                                (this->nodal_coordinates[(bound_id + 2) % 3][GlobalCoord::x] -
                                 this->nodal_coordinates[(bound_id + 1) % 3][GlobalCoord::x]) /
                                length);

    return surface_normal;
}

Array2D<double> StraightTriangle::GetPsi(const std::vector<Point<2>>& points) {
    Array2D<double> psi(3, std::vector<double>(points.size()));

    for (uint pt = 0; pt < points.size(); pt++) {
        psi[0][pt] = -(points[pt][LocalCoordTri::z1] + points[pt][LocalCoordTri::z2]) / 2;  // N1
        psi[1][pt] = (1 + points[pt][LocalCoordTri::z1]) / 2;                               // N2
        psi[2][pt] = (1 + points[pt][LocalCoordTri::z2]) / 2;                               // N3
    }

    return psi;
}

Array3D<double> StraightTriangle::GetDPsi(const std::vector<Point<2>>& points) {
    Array3D<double> dpsi(3, Array2D<double>(2, std::vector<double>(points.size())));

    Array3D<double> J_inv = this->GetJinv(points);

    for (uint pt = 0; pt < points.size(); pt++) {
        dpsi[0][GlobalCoord::x][pt] =
            -0.5 * J_inv[LocalCoordTri::z1][GlobalCoord::x][0] - 0.5 * J_inv[LocalCoordTri::z2][GlobalCoord::x][0];
        dpsi[0][GlobalCoord::y][pt] =
            -0.5 * J_inv[LocalCoordTri::z1][GlobalCoord::y][0] - 0.5 * J_inv[LocalCoordTri::z2][GlobalCoord::y][0];

        dpsi[1][GlobalCoord::x][pt] = 0.5 * J_inv[LocalCoordTri::z1][GlobalCoord::x][0];
        dpsi[1][GlobalCoord::y][pt] = 0.5 * J_inv[LocalCoordTri::z1][GlobalCoord::y][0];

        dpsi[2][GlobalCoord::x][pt] = 0.5 * J_inv[LocalCoordTri::z2][GlobalCoord::x][0];
        dpsi[2][GlobalCoord::y][pt] = 0.5 * J_inv[LocalCoordTri::z2][GlobalCoord::y][0];
    }

    return dpsi;
}

Array2D<double> StraightTriangle::GetBoundaryPsi(const uint bound_id, const std::vector<Point<1>>& points) {
    Array2D<double> psi_bound(2, std::vector<double>(points.size()));

    for (uint pt = 0; pt < points.size(); pt++) {
        psi_bound[0][pt] = (1 - points[pt][LocalCoordTri::z1]) / 2;  // N1
        psi_bound[1][pt] = (1 + points[pt][LocalCoordTri::z1]) / 2;  // N2
    }

    return psi_bound;
}

std::vector<Point<2>> StraightTriangle::LocalToGlobalCoordinates(const std::vector<Point<2>>& points) {
    assert(this->nodal_coordinates.size() > 0);
    std::vector<Point<2>> global_coordinates(points.size());

    Array2D<double> psi_pts = this->GetPsi(points);

    for (uint pt = 0; pt < points.size(); pt++) {
        global_coordinates[pt][GlobalCoord::x] = this->nodal_coordinates[0][GlobalCoord::x] * psi_pts[0][pt] +
                                                 this->nodal_coordinates[1][GlobalCoord::x] * psi_pts[1][pt] +
                                                 this->nodal_coordinates[2][GlobalCoord::x] * psi_pts[2][pt];

        global_coordinates[pt][GlobalCoord::y] = this->nodal_coordinates[0][GlobalCoord::y] * psi_pts[0][pt] +
                                                 this->nodal_coordinates[1][GlobalCoord::y] * psi_pts[1][pt] +
                                                 this->nodal_coordinates[2][GlobalCoord::y] * psi_pts[2][pt];
    }

    return global_coordinates;
}

void StraightTriangle::GetVTK(std::vector<Point<3>>& points, Array2D<uint>& cells) {
    assert(this->nodal_coordinates.size() > 0);
    uint number_pt = points.size();

    double z1;
    double z2;
    double dz = 2.0 / N_DIV;

    for (uint i = 0; i <= N_DIV; i++) {
        for (uint j = 0; j <= N_DIV - i; j++) {
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

    for (uint i = 0; i < N_DIV; i++) {
        for (uint j = 0; j < N_DIV - i; j++) {
            cells.push_back(std::vector<uint>(4));

            pt_ID = number_pt + (N_DIV + 1) * (N_DIV + 2) / 2 - (N_DIV - i + 1) * (N_DIV - i + 2) / 2 + j;

            cells.back()[0] = VTKElementTypes::straight_triangle;
            cells.back()[1] = pt_ID;
            cells.back()[2] = pt_ID + 1;
            cells.back()[3] = pt_ID + (N_DIV + 1 - i);
        }
    }

    for (uint i = 1; i < N_DIV; i++) {
        for (uint j = 0; j < N_DIV - i; j++) {
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
