<<<<<<< b6e108c4d075e33e1ea0e6068557db30e44f4337


int main() {
    using Utilities::almost_equal;
    bool error_found = false;


=======
#include "test_element_triangle.hpp"

int main() {
>>>>>>> Adding serialization for `Geometry::Element`
    // make an equilateral triangle
    std::vector<Point<2>> vrtxs(3);
    vrtxs[0] = {-0.5, 0.};
    vrtxs[1] = {0.5, 0.};
    vrtxs[2] = {0, std::sqrt(3.) / 2.};

    MasterType master(10);
    ShapeType shape(vrtxs);

    ElementType triangle(0, master, vrtxs, std::vector<uint>(0), std::vector<uint>(0), std::vector<unsigned char>(0));

    Integration::Dunavant_2D integ;
    std::vector<Point<2>> gp = integ.GetRule(20).second;

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

    if ( error_found ) {
        return 1;
    }
    return 0;
}