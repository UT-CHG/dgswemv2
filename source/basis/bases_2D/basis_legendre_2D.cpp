#include "../bases_2D.hpp"

namespace Basis {
Array2D<double> Legendre_2D::GetPhi(const uint p, const std::vector<Point<2>>& points) {
    uint n_pts = points.size();

    Array2D<double> phi((p + 1) * (p + 2) / 2);

    std::vector<double> n1(n_pts);
    std::vector<double> n2(n_pts);
    std::vector<double> Px(n_pts);
    std::vector<double> Py(n_pts);
    std::vector<double> ans(n_pts);



    for (uint pt = 0; pt < n_pts; pt++) {
        n1[pt] = points[pt][LocalCoordTri::z1];
        n2[pt] = points[pt][LocalCoordTri::z2];
    }

    
    for (uint dof = 0; dof < phi.size(); dof++) {
        uint tri_num_indx = (uint)std::ceil((-3. + std::sqrt(1. + 8. * (dof + 1))) / 2.);
        uint lower_tri_num = (tri_num_indx + 1) * tri_num_indx / 2;

        uint p = dof - lower_tri_num;
        uint q = tri_num_indx - p;

        
        Px=jacobi_polynomial(p,0,0,n1);
        Py=jacobi_polynomial(q,0,0,n2);

        for (uint count = 0; count < phi.size(); count++) {

            ans[count]=Px[count]*Py[count];

        }

        phi[dof]=ans;

    }

    return  phi;
}

Array3D<double> Legendre_2D::GetDPhi(const uint p, const std::vector<Point<2>>& points) {
    uint n_pts = points.size();

    Array3D<double> dphi((p + 1) * (p + 2) / 2);

    std::vector<double> n1(n_pts);
    std::vector<double> n2(n_pts);
    std::vector<double> Px(n_pts);
    std::vector<double> Pdx(n_pts);
    std::vector<double> Py(n_pts);
    std::vector<double> Pdy(n_pts);
    std::vector<double> Ux(n_pts);
    std::vector<double> Uy(n_pts);


    for (uint pt = 0; pt < n_pts; pt++) {
        n1[pt] = points[pt][LocalCoordTri::z1];
        n2[pt] = points[pt][LocalCoordTri::z2];
    }

    
    for (uint dof = 0; dof < dphi.size(); dof++) {
        uint tri_num_indx = (uint)std::ceil((-3. + std::sqrt(1. + 8. * (dof + 1))) / 2.);
        uint lower_tri_num = (tri_num_indx + 1) * tri_num_indx / 2;

        uint p = dof - lower_tri_num;
        uint q = tri_num_indx - p;

        
        Px=jacobi_polynomial(p,0,0,n1);
        Pdx=jacobi_polynomial_derivative(p,0,0,n1);
        Py=jacobi_polynomial(q,0,0,n2);
        Pdy=jacobi_polynomial_derivative(q,0,0,n2);

        for (uint count = 0; count < dphi.size(); count++) {

            Ux[count]=Px[count]*Pdy[count];
            Uy[count]=Pdx[count]*Py[count];

        }

        dphi[dof][0]=Ux;
        dphi[dof][1]=Uy;

    }

    return dphi;
}

std::pair<bool, Array2D<double>> Legendre_2D::GetMinv(const uint p) {
    return std::pair<bool, Array2D<double>>{};
}
}