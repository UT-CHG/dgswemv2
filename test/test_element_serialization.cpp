#include "test_element_triangle.hpp"

int main() {
    // make an equilateral triangle
    AlignedVector<Point<3>> vrtxs(3);
    vrtxs[0] = {-0.5, 0., 0.};
    vrtxs[1] = {0.5, 0., 0.};
    vrtxs[2] = {0, std::sqrt(3.) / 2., 0.};

    MasterType master(10);

    ElementType o_triangle(
        0, master, std::move(vrtxs), std::vector<uint>(0), std::vector<uint>(0), std::vector<uchar>(0));

    Integration::Dunavant_2D integ;
    AlignedVector<Point<2>> gp = integ.GetRule(20).second;

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

    x = o_triangle.ComputeNodalUgp(x_node);
    y = o_triangle.ComputeNodalUgp(y_node);

    DynMatrix<double> f_vals(1, ngp);
    for (uint gp = 0; gp < ngp; ++gp) {
        f_vals(0, gp) = std::pow(x(0, gp) + 1., 2) + std::pow(y(0, gp) - 1., 2);
    }

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_triangle;

    hpx::serialization::input_archive i_archive(buffer);
    ElementType i_triangle;
    i_archive >> i_triangle;
    i_triangle.SetMaster(master);
    i_triangle.Initialize();

    bool error_found = check_for_error(i_triangle, f_vals);

    if (error_found) {
        return 1;
    }
    return 0;
}