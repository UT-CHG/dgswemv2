#include "general_definitions.hpp"
#include "utilities/almost_equal.hpp"
#include "geometry/mesh_definitions.hpp"
#include "problem/SWE/swe_problem.hpp"

const std::vector<double> IntegrationPhi_true = {
    6.881941874331419e-01, -8.89156081756484e-02, 7.216878364870322e-02, 1.20281306081172e-02, 0.00000000000000e+00,
    4.811252243246882e-03, 0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00, 0.00000000000000e+00,
    0.00000000000000e+00};

const std::vector<double> IntegrationDPhiDX_true = {
    0.00000000000000e+00,  0.00000000000000e+00,  1.376388374866284,     0.00000000000000e+00,  6.212068893253612e-01,
    4.330127018922193e-01, 0.00000000000000e+00,  4.541823898251592e-01, 1.732050807568877e-01, 8.09844822527753e-01,
    0.00000000000000e+00,  3.593266739736606e-01, 8.66025403784439e-02,  4.752989329117831e-01, 2.886751345948129e-01,
    0.00000000000000e+00,  2.976844657418943e-01, 4.948716593053934e-02, 3.477531724792759e-01, 1.649572197684645e-01,
    5.677806551742285e-01, 0.00000000000000e+00,  2.542695064773368e-01, 3.092947870658709e-02, 2.799779589772887e-01,
    1.030982623552903e-01, 3.839882003961648e-01, 2.165063509461096e-01, 0.00000000000000e+00,  2.219874007710344e-01,
    2.061965247105806e-02, 2.368352966075515e-01, 6.873217490352688e-02, 2.924290225727635e-01, 1.443375672974064e-01,
    4.382557445972097e-01, 0.00000000000000e+00,  1.970208683922866e-01, 1.443375672974064e-02, 2.063521709566122e-01,
    4.811252243246881e-02, 2.388944320048427e-01, 1.010362971081845e-01, 3.216406511354539e-01, 1.732050807568877e-01,
    0.00000000000000e+00,  1.771259393312325e-01, 1.049727762162956e-02, 1.833727084180284e-01, 3.499092540543186e-02,
    2.038020106931869e-01, 7.348094335140692e-02, 2.541597625891523e-01, 1.259673314595547e-01, 3.571900656194553e-01,
    0.00000000000000e+00,  1.6089440460229e-01,   7.87295821622217e-03,  1.652850866122932e-01, 2.624319405407389e-02,
    1.788406944880994e-01, 5.511070751355518e-02, 2.112510229573667e-01, 9.4475498594666e-02,   2.765124421822678e-01,
    1.443375672974064e-01};

const std::vector<double> IntegrationDPhiDY_true = {
    0.00000000000000e+00,  2.383974596215561,      0.00000000000000e+00,  -1.744016935856293,    4.166666666666668e-01,
    5.639957660359269e-01, 2.337820631441114,      -3.00000000000000e-01, 3.443590455313351e-01, 2.00000000000000e-01,
    -1.879743494847109,    3.50000000000000e-01,   1.913828550551446e-01, 1.00000000000000e-01,  3.851221159381804e-01,
    2.252991532071854,     -3.238095238095238e-01, 2.131929837166789e-01, 5.714285714285716e-02, 2.453297042212707e-01,
    1.428571428571429e-01, -1.937912020128888,     3.392857142857143e-01, 1.178190458792872e-01, 3.571428571428573e-02,
    1.83043385083546e-01,  8.92857142857143e-02,   2.855387752616669e-01, 2.210576982387224,     -3.293650793650794e-01,
    1.641331427947897e-01, 2.380952380952382e-02,  1.483567763383696e-01, 5.952380952380954e-02, 2.021367719068524e-01,
    1.111111111111111e-01, -1.970227867507653,     3.361111111111112e-01, 8.21387896549881e-02,  1.666666666666667e-02,
    1.26044828084989e-01,  4.166666666666668e-02,  1.569447932368951e-01, 7.777777777777779e-02, 2.272196374916589e-01,
    2.185128252576446,     -3.313131313131313e-01, 1.367871463596913e-01, 1.212121212121213e-02, 1.102706931333926e-01,
    3.030303030303032e-02, 1.293162666349232e-01,  5.656565656565658e-02, 1.71701644642061e-01,  9.09090909090909e-02,
    -1.990792497657777,    3.348484848484849e-01,  6.041767889566378e-02, 9.0909090909091e-03,   9.83877593333601e-02,
    2.272727272727273e-02, 1.1080228133985e-01,    4.242424242424243e-02, 1.378321850860735e-01, 6.818181818181821e-02,
    1.8880147989604e-01};

int main() {
    using Utilities::almost_equal;
    bool error_found = false;

    using MasterType = Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>;
    using ShapeType = Shape::StraightTriangle;
    using ElementType = Geometry::Element<2, MasterType, ShapeType, SWE::Data>;

    // make an equilateral triangle
    std::vector<Point<2>> vrtxs(3);
    vrtxs[0] = {-0.5, 0.};
    vrtxs[1] = {0.5, 0.};
    vrtxs[2] = {0, std::sqrt(3.) / 2.};

    MasterType master(10);
    ShapeType shape(vrtxs);

    ElementType triangle(0, master, vrtxs, std::vector<uint>(0), std::vector<unsigned char>(0));

    // Check integrations
    Integration::Dunavant_2D integ;
    std::vector<Point<2>> gp = integ.GetRule(20).second;

    std::vector<double> x = shape.InterpolateNodalValues({-0.5, 0.5, 0}, gp);
    std::vector<double> y = shape.InterpolateNodalValues({0, 0, std::sqrt(3.) / 2.}, gp);

    std::vector<double> f_vals(triangle.data.get_ngp_internal());

    for (uint gp = 0; gp < triangle.data.get_ngp_internal(); gp++) {
        f_vals[gp] = std::pow(x[gp] + 1., 2) + std::pow(y[gp] - 1., 2);
    }

    if (!almost_equal(IntegrationPhi_true[0], triangle.Integration(f_vals), 1.e+04)) {
        error_found = true;
        std::cerr << "Error found in Triangle element in Integration" << std::endl;
    }

    for (uint dof = 0; dof < 66; dof++) {
        if (!almost_equal(IntegrationPhi_true[dof], triangle.IntegrationPhi(dof, f_vals), 1.e+04)) {
            error_found = true;
            std::cerr << "Error found in Triangle element in IntegrationPhi" << std::endl;
        }
    }

    for (uint dof = 0; dof < 66; dof++) {
        if (!almost_equal(IntegrationDPhiDX_true[dof], triangle.IntegrationDPhi(GlobalCoord::x, dof, f_vals), 1.e+04)) {
            error_found = true;

            std::cerr << "Error found in Triangle element in IntegrationDPhi "
                         "in x direction" << std::endl;
        }
    }
    // Add 7 more modes
    for (uint dof = 0; dof < 66; dof++) {
        if (!almost_equal(IntegrationDPhiDY_true[dof], triangle.IntegrationDPhi(GlobalCoord::y, dof, f_vals), 1.e+04)) {
            error_found = true;

            std::cerr << "Error found in Triangle element in IntegrationDPhi "
                         "in y direction" << std::endl;
        }
    }

    // Check ComputeUgp and SolveLSE
    std::vector<double> mod_vals(triangle.data.get_ndof());
    std::vector<double> gp_vals(triangle.data.get_ngp_internal());

    for (uint dof = 0; dof < 66; dof++) {
        std::fill(mod_vals.begin(), mod_vals.end(), 0.0);
        mod_vals[dof] = 1.0;

        triangle.ComputeUgp(mod_vals, gp_vals);

        for (uint doff = 0; doff < 66; doff++) {
            if (dof == doff) {
                if (!almost_equal(
                         (1. / triangle.IntegrationPhi(doff, gp_vals)), triangle.SolveLSE(mod_vals)[dof], 1.e+03)) {
                    error_found = true;

                    std::cerr << "Error found in Triangle element in SolveLSE" << std::endl;
                }
            } else {
                if (!almost_equal(triangle.IntegrationPhi(doff, gp_vals), 0.0)) {
                    error_found = true;

                    std::cerr << "Error found in Triangle element in ComputeUgp" << std::endl;
                }
            }
        }
    }

    // Check ComputeDUgp
    std::vector<double> gp_dvals(triangle.data.get_ngp_internal());
    std::fill(gp_vals.begin(), gp_vals.end(), 1.0);

    for (uint dof = 0; dof < 66; dof++) {
        std::fill(mod_vals.begin(), mod_vals.end(), 0.0);
        mod_vals[dof] = 1.0;

        triangle.ComputeDUgp(GlobalCoord::x, mod_vals, gp_dvals);

        if (!almost_equal(triangle.IntegrationPhi(0, gp_dvals),
                          triangle.IntegrationDPhi(GlobalCoord::x, dof, gp_vals))) {
            error_found = true;

            std::cerr << "Error found in Triangle element in ComputeDUgp in x "
                         "direction" << std::endl;
        }

        triangle.ComputeDUgp(GlobalCoord::y, mod_vals, gp_dvals);

        if (!almost_equal(triangle.IntegrationPhi(0, gp_dvals),
                          triangle.IntegrationDPhi(GlobalCoord::y, dof, gp_vals))) {
            error_found = true;

            std::cerr << "Error found in Triangle element in ComputeDUgp in y "
                         "direction" << std::endl;
        }
    }

    // Check L2 projection
    std::vector<double> nodal_vals{1.0, 2.0, 3.0};
    std::vector<double> modal_vals_true{2.0, 0.5, 0.5};

    std::vector<double> modal_vals_computed = triangle.L2Projection(nodal_vals);

    for (uint i = 0; i < 3; i++) {
        if (!almost_equal(modal_vals_computed[i], modal_vals_true[i])) {
            error_found = true;

            std::cerr << "Error found in Triangle element in L2 projection" << std::endl;
        }
    }

    for (uint i = 3; i < 66; i++) {
        if (!almost_equal(modal_vals_computed[i], 0.0, 1.e+04)) {
            error_found = true;

            std::cerr << "Error found in Triangle element in L2 projection" << std::endl;
        }
    }

    if (error_found) {
        return 1;
    }

    return 0;
}