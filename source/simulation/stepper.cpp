#include "stepper.hpp"

Stepper::Stepper(const StepperInput& stepper_input)
    : order(stepper_input.order),
      nstages(stepper_input.nstages),
      dt(stepper_input.dt),
      step(0),
      stage(0),
      timestamp(0),
      t(0.),
      ramp(0.) {
    // Allocate the time stepping arrays
    this->ark.reserve(this->nstages);
    this->brk.reserve(this->nstages);
    this->crk.reserve(this->nstages);
    this->drk = std::vector<double>(this->nstages, 0);

    for (uint step = 1; step <= this->nstages; ++step) {
        this->ark.emplace_back(std::vector<double>(step, 0));
        this->brk.emplace_back(std::vector<double>(step, 0));
        this->crk.emplace_back(std::vector<double>(step, 0));
    }

    // The forward Euler method
    if ((this->nstages == 1) && (this->order == 1)) {
        this->ark[0][0] = 1.;
        this->brk[0][0] = 1.;

        // SSP(s,2) schemes
    } else if ((this->nstages == 2) && (this->order == 2)) {
        this->ark[0][0] = 1;
        this->ark[1][0] = 0.5;
        this->ark[1][1] = 0.5;

        this->brk[0][0] = 1.;
        this->brk[1][1] = 0.5;

        // SSP(3,3) scheme
    } else if ((this->nstages == 3) && (this->order == 3)) {
        this->ark[0][0] = 1.;
        this->ark[1][0] = 3. / 4.;
        this->ark[1][1] = 1. / 4.;
        this->ark[2][0] = 1. / 3.;
        this->ark[2][2] = 2. / 3.;

        this->brk[0][0] = 1.;
        this->brk[1][1] = 1. / 4.;
        this->brk[2][2] = 2. / 3.;

        // SP(4,3) scheme
    } else if ((this->nstages == 4) && (this->order == 3)) {
        this->ark[0][0] = 1.;
        this->ark[1][1] = 1.;
        this->ark[2][0] = 2. / 3.;
        this->ark[2][2] = 1. / 3.;
        this->ark[3][3] = 1.;

        this->brk[0][0] = 1. / 2.;
        this->brk[1][1] = 1. / 2.;
        this->brk[2][2] = 1. / 6.;
        this->brk[3][3] = 1. / 2.;

        // SSP(5,3) scheme
    } else if ((this->nstages == 5) && (this->order == 3)) {
        this->ark[0][0] = 1.;
        this->ark[1][1] = 1.;
        this->ark[2][0] = 0.355909775063327;
        this->ark[2][2] = 0.644090224936674;
        this->ark[3][0] = 0.367933791638137;
        this->ark[3][3] = 0.632066208361863;
        this->ark[4][2] = 0.237593836598569;
        this->ark[4][4] = 0.762406163401431;

        this->brk[0][0] = 0.377268915331368;
        this->brk[1][1] = 0.377268915331368;
        this->brk[2][2] = 0.242995220537396;
        this->brk[3][3] = 0.238458932846290;
        this->brk[4][4] = 0.287632146308408;

        // SSP(6,3) scheme
    } else if ((this->nstages == 6) && (this->order == 3)) {
        this->ark[0][0] = 1.;
        this->ark[1][1] = 1.;
        this->ark[2][2] = 1.;
        this->ark[3][0] = 0.476769811285196;
        this->ark[3][1] = 0.098511733286064;
        this->ark[3][3] = 0.424718455428740;
        this->ark[4][4] = 1.;
        this->ark[5][2] = 0.155221702560091;
        this->ark[5][5] = 0.844778297439909;

        this->brk[0][0] = 0.284220721334261;
        this->brk[1][1] = 0.284220721334261;
        this->brk[2][2] = 0.284220721334261;
        this->brk[3][3] = 0.120713785765930;
        this->brk[4][4] = 0.284220721334261;
        this->brk[5][5] = 0.240103497065900;

        // SSP(7,3) scheme
    } else if ((this->nstages == 7) && (this->order == 3)) {
        this->ark[0][0] = 1.;
        this->ark[1][1] = 1.;
        this->ark[2][2] = 1.;
        this->ark[3][0] = 0.184962588071072;
        this->ark[3][3] = 0.815037411928928;
        this->ark[4][0] = 0.180718656570380;
        this->ark[4][1] = 0.314831034403793;
        this->ark[4][4] = 0.504450309025826;
        this->ark[5][5] = 1.;
        this->ark[6][3] = 0.120199000000000;
        this->ark[6][6] = 0.879801000000000;

        this->brk[0][0] = 0.233213863663009;
        this->brk[1][1] = 0.233213863663009;
        this->brk[2][2] = 0.233213863663009;
        this->brk[3][3] = 0.190078023865845;
        this->brk[4][4] = 0.117644805593912;
        this->brk[5][5] = 0.233213863663009;
        this->brk[6][6] = 0.205181790464579;

        // SSP(8,3) scheme
    } else if ((this->nstages == 8) && (this->order == 3)) {
        this->ark[0][0] = 1.;
        this->ark[1][1] = 1.;
        this->ark[2][2] = 1.;
        this->ark[3][3] = 1.;
        this->ark[4][0] = 0.421366967085359;
        this->ark[4][1] = 0.005949401107575;
        this->ark[4][4] = 0.572683631807067;
        this->ark[5][1] = 0.004254010666365;
        this->ark[5][5] = 0.995745989333635;
        this->ark[6][2] = 0.104380143093325;
        this->ark[6][3] = 0.243265240906726;
        this->ark[6][6] = 0.652354615999950;
        this->ark[7][7] = 1.;

        this->brk[0][0] = 0.195804015330143;
        this->brk[1][1] = 0.195804015330143;
        this->brk[2][2] = 0.195804015330143;
        this->brk[3][3] = 0.195804015330143;
        this->brk[4][4] = 0.112133754621673;
        this->brk[5][5] = 0.194971062960412;
        this->brk[6][6] = 0.127733653231944;
        this->brk[7][7] = 0.195804015330143;

        // SSP(5,4) scheme
    } else if ((this->nstages == 5) && (this->order == 4)) {
        this->ark[0][0] = 1.;
        this->ark[1][0] = 0.44437049406734;
        this->ark[1][1] = 0.55562950593266;
        this->ark[2][0] = 0.62010185138540;
        this->ark[2][2] = 0.37989814861460;
        this->ark[3][0] = 0.17807995410773;
        this->ark[3][3] = 0.82192004589227;
        this->ark[4][0] = 0.00683325884039;
        this->ark[4][2] = 0.51723167208978;
        this->ark[4][3] = 0.12759831133288;
        this->ark[4][4] = 0.34833675773694;

        this->brk[0][0] = 0.39175222700392;
        this->brk[1][1] = 0.36841059262959;
        this->brk[2][2] = 0.25189177424738;
        this->brk[3][3] = 0.54497475021237;
        this->brk[4][3] = 0.08460416338212;
        this->brk[4][4] = 0.22600748319395;

        // SSP(6,4) scheme
    } else if ((this->nstages == 6) && (this->order == 4)) {
        this->ark[0][0] = 1.00000000000000;
        this->ark[1][0] = 0.30948026455053;
        this->ark[1][1] = 0.69051973544947;
        this->ark[2][0] = 0.54205244285557;
        this->ark[2][2] = 0.45794755714443;
        this->ark[3][0] = 0.35984960863377;
        this->ark[3][3] = 0.64015039136623;
        this->ark[4][4] = 1.00000000000000;
        this->ark[5][0] = 0.05776282890116;
        this->ark[5][2] = 0.44216432622405;
        this->ark[5][4] = 0.10115567086469;
        this->ark[5][5] = 0.39891717401009;

        this->brk[0][0] = 0.39270746575722;
        this->brk[1][1] = 0.30154043149172;
        this->brk[2][2] = 0.19997937335132;
        this->brk[3][3] = 0.27954483459696;
        this->brk[4][4] = 0.43668618869443;
        this->brk[5][2] = 0.09150931531680;
        this->brk[5][4] = 0.04417328437472;
        this->brk[5][5] = 0.14911300530736;

        // SSP(7,4) scheme
    } else if ((this->nstages == 7) && (this->order == 4)) {
        this->ark[0][0] = 1.;
        this->ark[1][0] = 0.20161507213829;
        this->ark[1][1] = 0.79838492786171;
        this->ark[2][0] = 0.19469598207921;
        this->ark[2][2] = 0.80530401792079;
        this->ark[3][0] = 0.58143386885601;
        this->ark[3][3] = 0.41856613114399;
        this->ark[4][0] = 0.01934367892154;
        this->ark[4][4] = 0.98065632107846;
        this->ark[5][5] = 1.;
        this->ark[6][0] = 0.06006304558847;
        this->ark[6][2] = 0.30152730794242;
        this->ark[6][3] = 0.10518998496676;
        this->ark[6][4] = 0.01483791154585;
        this->ark[6][6] = 0.51838174995650;

        this->brk[0][0] = 0.30111872706068;
        this->brk[1][1] = 0.24040865318216;
        this->brk[2][2] = 0.24249212077315;
        this->brk[3][3] = 0.12603810060080;
        this->brk[4][4] = 0.29529398308716;
        this->brk[5][5] = 0.30111872706068;
        this->brk[6][2] = 0.09079551914158;
        this->brk[6][3] = 0.02888359354880;
        this->brk[6][6] = 0.15609445267839;

        // SSP(8,4) scheme
    } else if ((this->nstages == 8) && (this->order == 4)) {
        this->ark[0][0] = 1.;
        this->ark[1][0] = 0.10645325745007;
        this->ark[1][1] = 0.89354674254993;
        this->ark[2][2] = 1.;
        this->ark[3][0] = 0.57175518477257;
        this->ark[3][3] = 0.42824481522743;
        this->ark[4][0] = 0.19161667219044;
        this->ark[4][4] = 0.80838332780956;
        this->ark[5][5] = 1.;
        this->ark[6][6] = 1.;
        this->ark[7][0] = 0.02580435327923;
        this->ark[7][2] = 0.03629901341774;
        this->ark[7][3] = 0.31859181340256;
        this->ark[7][4] = 0.05186768980103;
        this->ark[7][5] = 0.03944076217320;
        this->ark[7][6] = 0.00511633747411;
        this->ark[7][7] = 0.52288003045213;

        this->brk[0][0] = 0.24120020561311;
        this->brk[1][1] = 0.21552365802797;
        this->brk[2][2] = 0.24120020561311;
        this->brk[3][3] = 0.10329273748560;
        this->brk[4][4] = 0.19498222488188;
        this->brk[5][5] = 0.24120020561311;
        this->brk[6][6] = 0.24120020561311;
        this->brk[7][2] = 0.00875532949991;
        this->brk[7][3] = 0.06195575835101;
        this->brk[7][5] = 0.00951311994571;
        this->brk[7][7] = 0.12611877085604;

    } else {
        throw std::logic_error("Fatal Error: invalid Runge-Kutta method entered!");
    }

    // Compute the time dependent parameters

    for (uint i = 0; i < this->nstages; ++i) {
        for (uint k = 0; k < i; ++k) {
            double casum = 0.;
            for (uint l = k + 1; l < i; ++l) {
                casum += this->crk.at(l).at(k) * this->ark.at(i).at(l);
            }
            this->crk.at(i).at(k) = this->brk.at(i).at(k) + casum;
        }
    }

    for (uint k = 1; k < this->nstages; ++k) {
        for (uint l = 0; l < k; ++l) {
            this->drk[k] = this->drk[k] + this->crk[k - 1][l];
        }
    }

    // Compute the maximum beta over alpha ratio at each stage

    /*for ( uint stage = 0; stage < this->nstages; ++stage ) {
        double max_boa = 0.;
        for ( uint i = 0; i <= stage; ++stage ) {
        if (this->ark(stage,i) != 0) {
            if (max_boa < this->brk[stage][i]/this->ark[stage][i]) { max_boa =
       this->brk[stage][i]/this->ark[stage][i]; }
        }
        }
        max_boa_dt[stage] = max_boa_dt*dt;
        }*/
}