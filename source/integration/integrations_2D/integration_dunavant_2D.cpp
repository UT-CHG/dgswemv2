#include "../integrations_2D.hpp"

namespace Integration {
std::pair<std::vector<double>, std::vector<Point<2>>> Dunavant_2D::GetRule(const uint p) {
    if (p < 0 || p > 20) {
        printf("\n");
        printf("DUNAVANT 2D - Fatal error!\n");
        printf("Illegal P = %d\n", p);
        exit(1);
    }

    std::vector<uint>                                     permutation = this->PermutationData(p);
    std::pair<std::vector<double>, std::vector<Point<3>>> gp          = this->GPData(p);

    if (permutation.size() != gp.first.size()) {
        printf("\n");
        printf("DUNAVANT 2D - Fatal error!\n");
        printf("Permutation size != GP size\n");
        exit(1);
    }

    uint number_gp = 0;
    for (uint i = 0; i < permutation.size(); i++)
        number_gp += permutation[i];

    std::pair<std::vector<double>, std::vector<Point<2>>> rule;
    rule.first.reserve(number_gp);
    rule.second.reserve(number_gp);

    for (uint i = 0; i < gp.first.size(); i++) {
        if (permutation[i] == 1) {
            rule.first.push_back(2 * gp.first[i]);
            rule.second.push_back({2 * gp.second[i][0] - 1, 2 * gp.second[i][1] - 1});
        } else if (permutation[i] == 3) {
            rule.first.push_back(2 * gp.first[i]);
            rule.second.push_back({2 * gp.second[i][0] - 1, 2 * gp.second[i][1] - 1});

            rule.first.push_back(2 * gp.first[i]);
            rule.second.push_back({2 * gp.second[i][1] - 1, 2 * gp.second[i][2] - 1});

            rule.first.push_back(2 * gp.first[i]);
            rule.second.push_back({2 * gp.second[i][2] - 1, 2 * gp.second[i][0] - 1});
        } else if (permutation[i] == 6) {
            rule.first.push_back(2 * gp.first[i]);
            rule.second.push_back({2 * gp.second[i][0] - 1, 2 * gp.second[i][1] - 1});

            rule.first.push_back(2 * gp.first[i]);
            rule.second.push_back({2 * gp.second[i][1] - 1, 2 * gp.second[i][2] - 1});

            rule.first.push_back(2 * gp.first[i]);
            rule.second.push_back({2 * gp.second[i][2] - 1, 2 * gp.second[i][0] - 1});

            rule.first.push_back(2 * gp.first[i]);
            rule.second.push_back({2 * gp.second[i][1] - 1, 2 * gp.second[i][0] - 1});

            rule.first.push_back(2 * gp.first[i]);
            rule.second.push_back({2 * gp.second[i][2] - 1, 2 * gp.second[i][1] - 1});

            rule.first.push_back(2 * gp.first[i]);
            rule.second.push_back({2 * gp.second[i][0] - 1, 2 * gp.second[i][2] - 1});
        }
    }

    return rule;
}

uint Dunavant_2D::GetNumGP(const uint p) {
    uint num_gp{0};
    switch (p) {
        case (1):
            num_gp = 1;
            break;
        case (2):
            num_gp = 3;
            break;
        case (3):
            num_gp = 4;
            break;
        case (4):
            num_gp = 6;
            break;
        case (5):
            num_gp = 7;
            break;
        case (6):
            num_gp = 12;
            break;
        case (7):
            num_gp = 13;
            break;
        case (8):
            num_gp = 16;
            break;
        case (9):
            num_gp = 19;
            break;
        case (10):
            num_gp = 25;
            break;
        case (11):
            num_gp = 27;
            break;
        case (12):
            num_gp = 33;
            break;
        case (13):
            num_gp = 37;
            break;
        case (14):
            num_gp = 42;
            break;
        case (15):
            num_gp = 48;
            break;
        case (16):
            num_gp = 52;
            break;
        case (17):
            num_gp = 61;
            break;
        case (18):
            num_gp = 70;
            break;
        case (19):
            num_gp = 73;
            break;
        case (20):
            num_gp = 79;
            break;
    }
    return num_gp;
}

std::vector<uint> Dunavant_2D::PermutationData(const uint p) {
    std::vector<uint> permutations;

    if (p == 1) {
        permutations.reserve(1);

        permutations.push_back(1);
    } else if (p == 2) {
        permutations.reserve(1);

        permutations.push_back(3);
    } else if (p == 3) {
        permutations.reserve(2);

        permutations.push_back(1);
        permutations.push_back(3);
    } else if (p == 4) {
        permutations.reserve(2);

        permutations.push_back(3);
        permutations.push_back(3);
    } else if (p == 5) {
        permutations.reserve(3);

        permutations.push_back(1);
        permutations.push_back(3);
        permutations.push_back(3);
    } else if (p == 6) {
        permutations.reserve(3);

        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
    } else if (p == 7) {
        permutations.reserve(4);

        permutations.push_back(1);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
    } else if (p == 8) {
        permutations.reserve(5);

        permutations.push_back(1);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
    } else if (p == 9) {
        permutations.reserve(6);

        permutations.push_back(1);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
    } else if (p == 10) {
        permutations.reserve(6);

        permutations.push_back(1);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
    } else if (p == 11) {
        permutations.reserve(7);

        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
        permutations.push_back(6);
    } else if (p == 12) {
        permutations.reserve(8);

        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
    } else if (p == 13) {
        permutations.reserve(10);

        permutations.push_back(1);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
    } else if (p == 14) {
        permutations.reserve(10);

        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
    } else if (p == 15) {
        permutations.reserve(11);

        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
    } else if (p == 16) {
        permutations.reserve(13);

        permutations.push_back(1);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
    } else if (p == 17) {
        permutations.reserve(15);

        permutations.push_back(1);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
    } else if (p == 18) {
        permutations.reserve(17);

        permutations.push_back(1);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
    } else if (p == 19) {
        permutations.reserve(17);

        permutations.push_back(1);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
    } else if (p == 20) {
        permutations.reserve(19);

        permutations.push_back(1);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(3);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
        permutations.push_back(6);
    }

    return permutations;
}

std::pair<std::vector<double>, std::vector<Point<3>>> Dunavant_2D::GPData(const uint p) {
    std::pair<std::vector<double>, std::vector<Point<3>>> gp;

    if (p == 1) {
        gp.first.reserve(1);
        gp.second.reserve(1);

        gp.first.push_back(1.000000000000000);
        gp.second.push_back({0.333333333333333, 0.333333333333333, 0.333333333333333});
    } else if (p == 2) {
        gp.first.reserve(1);
        gp.second.reserve(1);

        gp.first.push_back(0.333333333333333);
        gp.second.push_back({0.666666666666667, 0.166666666666667, 0.166666666666667});
    } else if (p == 3) {
        gp.first.reserve(2);
        gp.second.reserve(2);

        gp.first.push_back(-0.562500000000000);
        gp.second.push_back({0.333333333333333, 0.333333333333333, 0.333333333333333});

        gp.first.push_back(0.520833333333333);
        gp.second.push_back({0.600000000000000, 0.200000000000000, 0.200000000000000});
    } else if (p == 4) {
        gp.first.reserve(2);
        gp.second.reserve(2);

        gp.first.push_back(0.223381589678011);
        gp.second.push_back({0.108103018168070, 0.445948490915965, 0.445948490915965});

        gp.first.push_back(0.109951743655322);
        gp.second.push_back({0.816847572980459, 0.091576213509771, 0.091576213509771});
    } else if (p == 5) {
        gp.first.reserve(3);
        gp.second.reserve(3);

        gp.first.push_back(0.225000000000000);
        gp.second.push_back({0.333333333333333, 0.333333333333333, 0.333333333333333});

        gp.first.push_back(0.132394152788506);
        gp.second.push_back({0.059715871789770, 0.470142064105115, 0.470142064105115});

        gp.first.push_back(0.125939180544827);
        gp.second.push_back({0.797426985353087, 0.101286507323456, 0.101286507323456});
    } else if (p == 6) {
        gp.first.reserve(3);
        gp.second.reserve(3);

        gp.first.push_back(0.116786275726379);
        gp.second.push_back({0.501426509658179, 0.249286745170910, 0.249286745170910});

        gp.first.push_back(0.050844906370207);
        gp.second.push_back({0.873821971016996, 0.063089014491502, 0.063089014491502});

        gp.first.push_back(0.082851075618374);
        gp.second.push_back({0.053145049844817, 0.310352451033784, 0.636502499121399});
    } else if (p == 7) {
        gp.first.reserve(4);
        gp.second.reserve(4);

        gp.first.push_back(-0.149570044467682);
        gp.second.push_back({0.333333333333333, 0.333333333333333, 0.333333333333333});

        gp.first.push_back(0.175615257433208);
        gp.second.push_back({0.479308067841920, 0.260345966079040, 0.260345966079040});

        gp.first.push_back(0.053347235608838);
        gp.second.push_back({0.869739794195568, 0.065130102902216, 0.065130102902216});

        gp.first.push_back(0.077113760890257);
        gp.second.push_back({0.048690315425316, 0.312865496004874, 0.638444188569810});
    } else if (p == 8) {
        gp.first.reserve(5);
        gp.second.reserve(5);

        gp.first.push_back(0.144315607677787);
        gp.second.push_back({0.333333333333333, 0.333333333333333, 0.333333333333333});

        gp.first.push_back(0.095091634267285);
        gp.second.push_back({0.081414823414554, 0.459292588292723, 0.459292588292723});

        gp.first.push_back(0.103217370534718);
        gp.second.push_back({0.658861384496480, 0.170569307751760, 0.170569307751760});

        gp.first.push_back(0.032458497623198);
        gp.second.push_back({0.898905543365938, 0.050547228317031, 0.050547228317031});

        gp.first.push_back(0.027230314174435);
        gp.second.push_back({0.008394777409958, 0.263112829634638, 0.728492392955404});
    } else if (p == 9) {
        gp.first.reserve(6);
        gp.second.reserve(6);

        gp.first.push_back(0.097135796282799);
        gp.second.push_back({0.333333333333333, 0.333333333333333, 0.333333333333333});

        gp.first.push_back(0.031334700227139);
        gp.second.push_back({0.020634961602525, 0.489682519198738, 0.489682519198738});

        gp.first.push_back(0.077827541004774);
        gp.second.push_back({0.125820817014127, 0.437089591492937, 0.437089591492937});

        gp.first.push_back(0.079647738927210);
        gp.second.push_back({0.623592928761935, 0.188203535619033, 0.188203535619033});

        gp.first.push_back(0.025577675658698);
        gp.second.push_back({0.910540973211095, 0.044729513394453, 0.044729513394453});

        gp.first.push_back(0.043283539377289);
        gp.second.push_back({0.036838412054736, 0.221962989160766, 0.741198598784498});
    } else if (p == 10) {
        gp.first.reserve(6);
        gp.second.reserve(6);

        gp.first.push_back(0.090817990382754);
        gp.second.push_back({0.333333333333333, 0.333333333333333, 0.333333333333333});

        gp.first.push_back(0.036725957756467);
        gp.second.push_back({0.028844733232685, 0.485577633383657, 0.485577633383657});

        gp.first.push_back(0.045321059435528);
        gp.second.push_back({0.781036849029926, 0.109481575485037, 0.109481575485037});

        gp.first.push_back(0.072757916845420);
        gp.second.push_back({0.141707219414880, 0.307939838764121, 0.550352941820999});

        gp.first.push_back(0.028327242531057);
        gp.second.push_back({0.025003534762686, 0.246672560639903, 0.728323904597411});

        gp.first.push_back(0.009421666963733);
        gp.second.push_back({0.009540815400299, 0.066803251012200, 0.923655933587500});
    } else if (p == 11) {
        gp.first.reserve(7);
        gp.second.reserve(7);

        gp.first.push_back(0.000927006328961);
        gp.second.push_back({-0.069222096541517, 0.534611048270758, 0.534611048270758});

        gp.first.push_back(0.077149534914813);
        gp.second.push_back({0.202061394068290, 0.398969302965855, 0.398969302965855});

        gp.first.push_back(0.059322977380774);
        gp.second.push_back({0.593380199137435, 0.203309900431282, 0.203309900431282});

        gp.first.push_back(0.036184540503418);
        gp.second.push_back({0.761298175434837, 0.119350912282581, 0.119350912282581});

        gp.first.push_back(0.013659731002678);
        gp.second.push_back({0.935270103777448, 0.032364948111276, 0.032364948111276});

        gp.first.push_back(0.052337111962204);
        gp.second.push_back({0.050178138310495, 0.356620648261293, 0.593201213428213});

        gp.first.push_back(0.020707659639141);
        gp.second.push_back({0.021022016536166, 0.171488980304042, 0.807489003159792});
    } else if (p == 12) {
        gp.first.reserve(8);
        gp.second.reserve(8);

        gp.first.push_back(0.025731066440455);
        gp.second.push_back({0.02356522045239, 0.488217389773805, 0.488217389773805});

        gp.first.push_back(0.043692544538038);
        gp.second.push_back({0.120551215411079, 0.43972439229446, 0.43972439229446});

        gp.first.push_back(0.062858224217885);
        gp.second.push_back({0.457579229975768, 0.271210385012116, 0.271210385012116});

        gp.first.push_back(0.034796112930709);
        gp.second.push_back({0.744847708916828, 0.127576145541586, 0.127576145541586});

        gp.first.push_back(0.006166261051559);
        gp.second.push_back({0.957365299093579, 0.02131735045321, 0.02131735045321});

        gp.first.push_back(0.040371557766381);
        gp.second.push_back({0.115343494534698, 0.275713269685514, 0.608943235779788});

        gp.first.push_back(0.022356773202303);
        gp.second.push_back({0.022838332222257, 0.28132558098994, 0.695836086787803});

        gp.first.push_back(0.017316231108659);
        gp.second.push_back({0.02573405054833, 0.116251915907597, 0.858014033544073});
    } else if (p == 13) {
        gp.first.reserve(10);
        gp.second.reserve(10);

        gp.first.push_back(0.052520923400802);
        gp.second.push_back({0.333333333333333, 0.333333333333333, 0.333333333333333});

        gp.first.push_back(0.01128014520933);
        gp.second.push_back({0.009903630120591, 0.495048184939705, 0.495048184939705});

        gp.first.push_back(0.031423518362454);
        gp.second.push_back({0.062566729780852, 0.468716635109574, 0.468716635109574});

        gp.first.push_back(0.047072502504194);
        gp.second.push_back({0.170957326397447, 0.414521336801277, 0.414521336801277});

        gp.first.push_back(0.047363586536355);
        gp.second.push_back({0.541200855914337, 0.229399572042831, 0.229399572042831});

        gp.first.push_back(0.031167529045794);
        gp.second.push_back({0.77115100960734, 0.11442449519633, 0.11442449519633});

        gp.first.push_back(0.007975771465074);
        gp.second.push_back({0.950377217273082, 0.024811391363459, 0.024811391363459});

        gp.first.push_back(0.036848402728732);
        gp.second.push_back({0.094853828379579, 0.268794997058761, 0.63635117456166});

        gp.first.push_back(0.017401463303822);
        gp.second.push_back({0.018100773278807, 0.291730066734288, 0.690169159986905});

        gp.first.push_back(0.015521786839045);
        gp.second.push_back({0.02223307667409, 0.126357385491669, 0.851409537834241});
    } else if (p == 14) {
        gp.first.reserve(10);
        gp.second.reserve(10);

        gp.first.push_back(0.021883581369429);
        gp.second.push_back({0.022072179275643, 0.488963910362179, 0.488963910362179});

        gp.first.push_back(0.032788353544125);
        gp.second.push_back({0.164710561319092, 0.417644719340454, 0.417644719340454});

        gp.first.push_back(0.051774104507292);
        gp.second.push_back({0.453044943382323, 0.273477528308839, 0.273477528308839});

        gp.first.push_back(0.042162588736993);
        gp.second.push_back({0.645588935174913, 0.177205532412543, 0.177205532412543});

        gp.first.push_back(0.014433699669777);
        gp.second.push_back({0.876400233818255, 0.061799883090873, 0.061799883090873});

        gp.first.push_back(0.0049234036024);
        gp.second.push_back({0.961218077502598, 0.019390961248701, 0.019390961248701});

        gp.first.push_back(0.024665753212564);
        gp.second.push_back({0.057124757403648, 0.172266687821356, 0.770608554774996});

        gp.first.push_back(0.038571510787061);
        gp.second.push_back({0.092916249356972, 0.336861459796345, 0.570222290846683});

        gp.first.push_back(0.014436308113534);
        gp.second.push_back({0.014646950055654, 0.298372882136258, 0.686980167808088});

        gp.first.push_back(0.005010228838501);
        gp.second.push_back({0.001268330932872, 0.118974497696957, 0.879757171370171});
    } else if (p == 15) {
        gp.first.reserve(11);
        gp.second.reserve(11);

        gp.first.push_back(0.001916875642849);
        gp.second.push_back({-0.013945833716486, 0.506972916858243, 0.506972916858243});

        gp.first.push_back(0.044249027271145);
        gp.second.push_back({0.137187291433955, 0.431406354283023, 0.431406354283023});

        gp.first.push_back(0.051186548718852);
        gp.second.push_back({0.444612710305711, 0.277693644847144, 0.277693644847144});

        gp.first.push_back(0.023687735870688);
        gp.second.push_back({0.747070217917492, 0.126464891041254, 0.126464891041254});

        gp.first.push_back(0.013289775690021);
        gp.second.push_back({0.858383228050628, 0.070808385974686, 0.070808385974686});

        gp.first.push_back(0.004748916608192);
        gp.second.push_back({0.962069659517853, 0.018965170241073, 0.018965170241073});

        gp.first.push_back(0.038550072599593);
        gp.second.push_back({0.133734161966621, 0.261311371140087, 0.604954466893291});

        gp.first.push_back(0.027215814320624);
        gp.second.push_back({0.036366677396917, 0.388046767090269, 0.575586555512814});

        gp.first.push_back(0.002182077366797);
        gp.second.push_back({-0.010174883126571, 0.285712220049916, 0.724462663076655});

        gp.first.push_back(0.021505319847731);
        gp.second.push_back({0.036843869875878, 0.215599664072284, 0.747556466051838});

        gp.first.push_back(0.007673942631049);
        gp.second.push_back({0.012459809331199, 0.103575616576386, 0.883964574092416});
    } else if (p == 16) {
        gp.first.reserve(13);
        gp.second.reserve(13);

        gp.first.push_back(0.046875697427642);
        gp.second.push_back({0.333333333333333, 0.333333333333333, 0.333333333333333});

        gp.first.push_back(0.006405878578585);
        gp.second.push_back({0.005238916103123, 0.497380541948438, 0.497380541948438});

        gp.first.push_back(0.041710296739387);
        gp.second.push_back({0.173061122901295, 0.413469438549352, 0.413469438549352});

        gp.first.push_back(0.026891484250064);
        gp.second.push_back({0.059082801866017, 0.470458599066991, 0.470458599066991});

        gp.first.push_back(0.04213252276165);
        gp.second.push_back({0.518892500060958, 0.240553749969521, 0.240553749969521});

        gp.first.push_back(0.030000266842773);
        gp.second.push_back({0.704068411554854, 0.147965794222573, 0.147965794222573});

        gp.first.push_back(0.014200098925024);
        gp.second.push_back({0.849069624685052, 0.075465187657474, 0.075465187657474});

        gp.first.push_back(0.003582462351273);
        gp.second.push_back({0.96680719475395, 0.016596402623025, 0.016596402623025});

        gp.first.push_back(0.032773147460627);
        gp.second.push_back({0.103575692245252, 0.296555596579887, 0.599868711174861});

        gp.first.push_back(0.015298306248441);
        gp.second.push_back({0.020083411655416, 0.337723063403079, 0.642193524941505});

        gp.first.push_back(0.002386244192839);
        gp.second.push_back({-0.004341002614139, 0.204748281642812, 0.799592720971327});

        gp.first.push_back(0.019084792755899);
        gp.second.push_back({0.04194178646801, 0.189358492130623, 0.768699721401368});

        gp.first.push_back(0.006850054546542);
        gp.second.push_back({0.014317320230681, 0.085283615682657, 0.900399064086661});
    } else if (p == 17) {
        gp.first.reserve(15);
        gp.second.reserve(15);

        gp.first.push_back(0.033437199290803);
        gp.second.push_back({0.333333333333333, 0.333333333333333, 0.333333333333333});

        gp.first.push_back(0.005093415440507);
        gp.second.push_back({0.005658918886452, 0.497170540556774, 0.497170540556774});

        gp.first.push_back(0.014670864527638);
        gp.second.push_back({0.035647354750751, 0.482176322624625, 0.482176322624625});

        gp.first.push_back(0.024350878353672);
        gp.second.push_back({0.099520061958437, 0.450239969020782, 0.450239969020782});

        gp.first.push_back(0.031107550868969);
        gp.second.push_back({0.199467521245206, 0.400266239377397, 0.400266239377397});

        gp.first.push_back(0.03125711121862);
        gp.second.push_back({0.495717464058095, 0.252141267970953, 0.252141267970953});

        gp.first.push_back(0.024815654339665);
        gp.second.push_back({0.675905990683077, 0.162047004658461, 0.162047004658461});

        gp.first.push_back(0.014056073070557);
        gp.second.push_back({0.848248235478508, 0.075875882260746, 0.075875882260746});

        gp.first.push_back(0.003194676173779);
        gp.second.push_back({0.968690546064356, 0.015654726967822, 0.015654726967822});

        gp.first.push_back(0.008119655318993);
        gp.second.push_back({0.010186928826919, 0.334319867363658, 0.655493203809423});

        gp.first.push_back(0.026805742283163);
        gp.second.push_back({0.135440871671036, 0.292221537796944, 0.57233759053202});

        gp.first.push_back(0.018459993210822);
        gp.second.push_back({0.054423924290583, 0.31957488542319, 0.626001190286228});

        gp.first.push_back(0.008476868534328);
        gp.second.push_back({0.012868560833637, 0.190704224192292, 0.796427214974071});

        gp.first.push_back(0.018292796770025);
        gp.second.push_back({0.067165782413524, 0.180483211648746, 0.752351005937729});

        gp.first.push_back(0.006665632004165);
        gp.second.push_back({0.014663182224828, 0.080711313679564, 0.904625504095608});
    } else if (p == 18) {
        gp.first.reserve(17);
        gp.second.reserve(17);

        gp.first.push_back(0.030809939937647);
        gp.second.push_back({0.333333333333333, 0.333333333333333, 0.333333333333333});

        gp.first.push_back(0.009072436679404);
        gp.second.push_back({0.013310382738157, 0.493344808630921, 0.493344808630921});

        gp.first.push_back(0.018761316939594);
        gp.second.push_back({0.061578811516086, 0.469210594241957, 0.469210594241957});

        gp.first.push_back(0.019441097985477);
        gp.second.push_back({0.127437208225989, 0.436281395887006, 0.436281395887006});

        gp.first.push_back(0.02775394861081);
        gp.second.push_back({0.210307658653168, 0.394846170673416, 0.394846170673416});

        gp.first.push_back(0.032256225351457);
        gp.second.push_back({0.500410862393686, 0.249794568803157, 0.249794568803157});

        gp.first.push_back(0.025074032616922);
        gp.second.push_back({0.677135612512315, 0.161432193743843, 0.161432193743843});

        gp.first.push_back(0.015271927971832);
        gp.second.push_back({0.846803545029257, 0.076598227485371, 0.076598227485371});

        gp.first.push_back(0.006793922022963);
        gp.second.push_back({0.9514951212931, 0.02425243935345, 0.02425243935345});

        gp.first.push_back(-0.00222309872992);
        gp.second.push_back({0.913707265566071, 0.043146367216965, 0.043146367216965});

        gp.first.push_back(0.006331914076406);
        gp.second.push_back({0.00843053620242, 0.358911494940944, 0.632657968856636});

        gp.first.push_back(0.027257538049138);
        gp.second.push_back({0.131186551737188, 0.294402476751957, 0.574410971510855});

        gp.first.push_back(0.017676785649465);
        gp.second.push_back({0.050203151565675, 0.325017801641814, 0.624779046792512});

        gp.first.push_back(0.01837948463807);
        gp.second.push_back({0.066329263810916, 0.184737559666046, 0.748933176523037});

        gp.first.push_back(0.008104732808192);
        gp.second.push_back({0.011996194566236, 0.218796800013321, 0.769207005420443});

        gp.first.push_back(0.007634129070725);
        gp.second.push_back({0.014858100590125, 0.101179597136408, 0.883962302273467});

        gp.first.push_back(0.000046187660794);
        gp.second.push_back({-0.035222015287949, 0.020874755282586, 1.01434726000536});
    } else if (p == 19) {
        gp.first.reserve(17);
        gp.second.reserve(17);

        gp.first.push_back(0.032906331388919);
        gp.second.push_back({0.333333333333333, 0.333333333333333, 0.333333333333333});

        gp.first.push_back(0.010330731891272);
        gp.second.push_back({0.020780025853987, 0.489609987073006, 0.489609987073006});

        gp.first.push_back(0.022387247263016);
        gp.second.push_back({0.090926214604215, 0.454536892697893, 0.454536892697893});

        gp.first.push_back(0.030266125869468);
        gp.second.push_back({0.197166638701138, 0.401416680649431, 0.401416680649431});

        gp.first.push_back(0.030490967802198);
        gp.second.push_back({0.488896691193805, 0.255551654403098, 0.255551654403098});

        gp.first.push_back(0.024159212741641);
        gp.second.push_back({0.645844115695741, 0.17707794215213, 0.17707794215213});

        gp.first.push_back(0.016050803586801);
        gp.second.push_back({0.779877893544096, 0.110061053227952, 0.110061053227952});

        gp.first.push_back(0.008084580261784);
        gp.second.push_back({0.888942751496321, 0.05552862425184, 0.05552862425184});

        gp.first.push_back(0.002079362027485);
        gp.second.push_back({0.974756272445543, 0.012621863777229, 0.012621863777229});

        gp.first.push_back(0.003884876904981);
        gp.second.push_back({0.003611417848412, 0.395754787356943, 0.600633794794645});

        gp.first.push_back(0.025574160612022);
        gp.second.push_back({0.13446675453078, 0.307929983880436, 0.557603261588784});

        gp.first.push_back(0.008880903573338);
        gp.second.push_back({0.014446025776115, 0.26456694840652, 0.720987025817365});

        gp.first.push_back(0.016124546761731);
        gp.second.push_back({0.046933578838178, 0.358539352205951, 0.594527068955871});

        gp.first.push_back(0.002491941817491);
        gp.second.push_back({0.002861120350567, 0.157807405968595, 0.839331473680839});

        gp.first.push_back(0.018242840118951);
        gp.second.push_back({0.223861424097916, 0.075050596975911, 0.701087978926173});

        gp.first.push_back(0.010258563736199);
        gp.second.push_back({0.03464707481676, 0.142421601113383, 0.822931324069857});

        gp.first.push_back(0.003799928855302);
        gp.second.push_back({0.010161119296278, 0.065494628082938, 0.924344252620784});
    } else if (p == 20) {
        gp.first.reserve(19);
        gp.second.reserve(19);

        gp.first.push_back(0.033057055541624);
        gp.second.push_back({0.333333333333333, 0.333333333333333, 0.333333333333333});

        gp.first.push_back(0.000867019185663);
        gp.second.push_back({-0.0019009287044, 0.5009504643522, 0.5009504643522});

        gp.first.push_back(0.011660052716448);
        gp.second.push_back({0.023574084130543, 0.488212957934729, 0.488212957934729});

        gp.first.push_back(0.022876936356421);
        gp.second.push_back({0.089726636099435, 0.455136681950283, 0.455136681950283});

        gp.first.push_back(0.030448982673938);
        gp.second.push_back({0.196007481363421, 0.401996259318289, 0.401996259318289});

        gp.first.push_back(0.030624891725355);
        gp.second.push_back({0.488214180481157, 0.255892909759421, 0.255892909759421});

        gp.first.push_back(0.0243680576768);
        gp.second.push_back({0.647023488009788, 0.176488255995106, 0.176488255995106});

        gp.first.push_back(0.015997432032024);
        gp.second.push_back({0.791658289326483, 0.104170855336758, 0.104170855336758});

        gp.first.push_back(0.007698301815602);
        gp.second.push_back({0.89386207231814, 0.05306896384093, 0.05306896384093});

        gp.first.push_back(-0.000632060497488);
        gp.second.push_back({0.916762569607942, 0.041618715196029, 0.041618715196029});

        gp.first.push_back(0.001751134301193);
        gp.second.push_back({0.976836157186356, 0.011581921406822, 0.011581921406822});

        gp.first.push_back(0.016465839189576);
        gp.second.push_back({0.048741583664839, 0.344855770229001, 0.60640264610616});

        gp.first.push_back(0.004839033540485);
        gp.second.push_back({0.006314115948605, 0.377843269594854, 0.615842614456541});

        gp.first.push_back(0.02580490653465);
        gp.second.push_back({0.134316520547348, 0.306635479062357, 0.559048000390295});

        gp.first.push_back(0.008471091054441);
        gp.second.push_back({0.013973893962392, 0.249419362774742, 0.736606743262866});

        gp.first.push_back(0.01835491410628);
        gp.second.push_back({0.075549132909764, 0.212775724802802, 0.711675142287434});

        gp.first.push_back(0.000704404677908);
        gp.second.push_back({-0.008368153208227, 0.146965436053239, 0.861402717154987});

        gp.first.push_back(0.010112684927462);
        gp.second.push_back({0.026686063258714, 0.137726978828923, 0.835586957912363});

        gp.first.push_back(0.00357390938595);
        gp.second.push_back({0.010547719294141, 0.059696109149007, 0.929756171556853});
    }

    return gp;
}
}