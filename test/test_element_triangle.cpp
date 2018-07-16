#include "test_element_triangle.hpp"

int main() {
    // make an equilateral triangle
    DynVector<Point<3>> vrtxs(3);
    vrtxs[0] = {-0.5, 0., 0.};
    vrtxs[1] = {0.5, 0., 0.};
    vrtxs[2] = {0, std::sqrt(3.) / 2., 0.};

    MasterType master(10);

    ElementType triangle(0,
                         master,
                         std::move(vrtxs),
                         std::move(std::vector<uint>(0)),
                         std::move(std::vector<uint>(0)),
                         std::move(std::vector<unsigned char>(0)));

    Integration::Dunavant_2D integ;
    DynVector<Point<2>> gp = integ.GetRule(20).second;

    std::vector<double> x(gp.size());
    std::vector<double> y(gp.size());

    triangle.ComputeNodalUgp({-0.5, 0.5, 0}, x);
    triangle.ComputeNodalUgp({0, 0, std::sqrt(3.) / 2.}, y);

    std::size_t ngp = x.size();
    std::vector<double> f_vals(ngp);

    for (uint gp = 0; gp < ngp; gp++) {
        f_vals[gp] = std::pow(x[gp] + 1., 2) + std::pow(y[gp] - 1., 2);
    }

    bool error_found = check_for_error(triangle, f_vals);

    if (error_found) {
        return 1;
    }

    return 0;
}