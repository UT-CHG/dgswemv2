#include "test_element_triangle.hpp"

int main() {
    // make an equilateral triangle
    std::vector<Point<3>> vrtxs(3);
    vrtxs[0] = {-0.5, 0., 0.};
    vrtxs[1] = {0.5, 0., 0.};
    vrtxs[2] = {0, std::sqrt(3.) / 2., 0.};

    MasterType master(10);

    ElementType o_triangle(0,
                           master,
                           std::move(vrtxs),
                           std::move(std::vector<uint>(0)),
                           std::move(std::vector<uint>(0)),
                           std::move(std::vector<unsigned char>(0)));

    Integration::Dunavant_2D integ;
    std::vector<Point<2>> gp = integ.GetRule(20).second;

    std::vector<double> x(gp.size());
    std::vector<double> y(gp.size());

    o_triangle.ComputeNodalUgp({-0.5, 0.5, 0}, x);
    o_triangle.ComputeNodalUgp({0, 0, std::sqrt(3.) / 2.}, y);

    std::size_t ngp = x.size();
    std::vector<double> f_vals(ngp);

    for (uint gp = 0; gp < ngp; gp++) {
        f_vals[gp] = std::pow(x[gp] + 1., 2) + std::pow(y[gp] - 1., 2);
    }

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_triangle;

    hpx::serialization::input_archive i_archive(buffer);
    ElementType i_triangle;
    i_archive >> i_triangle;
    i_triangle.SetMaster(master);

    bool error_found = check_for_error(i_triangle, f_vals);

    if (error_found) {
        return 1;
    }
    return 0;
}