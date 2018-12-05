#include "test_element_triangle.hpp"

int main() {
    // make an equilateral triangle
    std::vector<Point<3>> vrtxs(3);
    vrtxs[0] = {-0.5, 0., 0.};
    vrtxs[1] = {0.5, 0., 0.};
    vrtxs[2] = {0, std::sqrt(3.) / 2., 0.};

    SWE::RKDG::SoAContainer data_holder(3 /*ndofs*/, 1 /*stage*/, 1 /*element*/);

    MasterType master(10);

    ElementType triangle(0,
                         master,
                         std::move(data_holder.at(0)),
                         std::move(vrtxs),
                         std::move(std::vector<uint>(0)),
                         std::move(std::vector<uint>(0)),
                         std::move(std::vector<unsigned char>(0)));

    Integration::Dunavant_2D integ;
    std::vector<Point<2>> gp = integ.GetRule(20).second;

    std::size_t ngp = gp.size();

    DynMatrix<double> x_node(1, 3);
    DynMatrix<double> y_node(1, 3);

    x_node(0, 0) = -0.5;
    x_node(0, 1) = 0.5;
    x_node(0, 2) = 0.0;

    y_node(0, 0) = 0.0;
    y_node(0, 1) = 0.0;
    y_node(0, 2) = std::sqrt(3.0) / 2.0;

    DynMatrix<double> x(1, ngp);
    DynMatrix<double> y(1, ngp);

    x = triangle.ComputeNodalUgp(x_node);
    y = triangle.ComputeNodalUgp(y_node);

    DynMatrix<double> f_vals(1, ngp);
    for (uint gp = 0; gp < ngp; ++gp) {
        f_vals(0, gp) = std::pow(x(0, gp) + 1., 2) + std::pow(y(0, gp) - 1., 2);
    }

    bool error_found = check_for_error(triangle, f_vals);

    if (error_found) {
        return 1;
    }

    return 0;
}