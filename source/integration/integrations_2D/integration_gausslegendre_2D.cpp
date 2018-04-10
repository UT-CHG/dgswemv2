#include "../integrations_2D.hpp"
#include "../integrations_1D.hpp"

namespace Integration {
std::pair<std::vector<double>, std::vector<Point<2>>> GaussLegendre_2D::GetRule(const uint p) {
    
    GaussLegendre_1D Rule1D;
    std::pair<std::vector<double>, std::vector<Point<1>>> Rule1Ddata=Rule1D.GetRule(p);
    std::pair<std::vector<double>, std::vector<Point<2>>>  Rule2Ddata;

    for (uint i=0;i<Rule1D.GetNumGP(p);++i) {

        for (uint j=0;j<Rule1D.GetNumGP(p);++j) {

            Rule2Ddata.first.push_back(Rule1Ddata.first[i]*Rule1Ddata.first[j]);
            Point<2>PTE;
            PTE[0]=Rule1Ddata.second[i][0];
            PTE[1]=Rule1Ddata.second[j][0];
            Rule2Ddata.second.push_back(PTE);

        }

    }
    
    return Rule2Ddata;
}

uint GaussLegendre_2D::GetNumGP(const uint p) {
    /* returns number of gauss points*/
    GaussLegendre_1D Rule1D;
    uint NumGP1D=Rule1D.GetNumGP(p);

    return NumGP1D*NumGP1D;
}
}