#include "integration/integrations_2D.hpp"
#include "utilities/almost_equal.hpp"

int main() {
    using Utilities::almost_equal;

    bool any_error = false;
    Integration::Dunavant_2D dunavant;
    std::pair<DynVector<double>, std::vector<Point<2>>> rule;

    for (uint p = 1; p < 20; ++p) {
        double exact_integration =
            1 / ((double)p + 1) *
            ((1 - pow(-1.0, p)) / ((double)p + 2) + 2 * pow(-1.0, p));  // S(x^p)dxdy over triangle

        rule = dunavant.GetRule(p);

        uint num_gp = dunavant.GetNumGP(p);

        double num_integration = 0;
        for (uint gp = 0; gp < rule.first.size(); ++gp) {
            num_integration += pow(rule.second[gp][GlobalCoord::x], p) * rule.first[gp];
        }

        if (!almost_equal(num_integration, exact_integration, 1.e+03)) {
            any_error = true;
            std::cerr << "Error found in Dunavant 2D at " << std::to_string(p)
                      << " - integration true value: " << exact_integration
                      << ", integration computed value: " << num_integration << std::endl;
        }

        if (num_gp != rule.first.size()) {
            any_error = true;
            std::cerr << "Error found in Dunavant 2D at " << std::to_string(p) << " gp_vector has size "
                      << rule.first.size() << " and return num_gp value " << num_gp << std::endl;
        }
    }

    if (any_error) {
        return 1;
    }

    return 0;
}