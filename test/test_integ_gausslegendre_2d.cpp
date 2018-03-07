#include "integration/integrations_2D.hpp"
#include "utilities/almost_equal.hpp"
#include <iostream>
#include <stdio.h>

int main() {


    Integration::GaussLegendre_2D gausslegendre;
    
    for (uint p=1;p<65;++p){

    std::pair<std::vector<double>, std::vector<Point<2>>> Rule2Ddata=gausslegendre.GetRule(p);
    double integration = 0;

    for (uint i=0;i<gausslegendre.GetNumGP(p);++i){

            double x=Rule2Ddata.second[i][0];
            double y=Rule2Ddata.second[i][1];
            double w=Rule2Ddata.first[i];
            integration+=(std::pow((x+1),p)+std::pow((y+1),p))*w;

    }

    if (!Utilities::almost_equal(integration,4*std::pow(2.0,p+1)/(p+1))) {

        return 1;

    }
    
}
    return 0;
}