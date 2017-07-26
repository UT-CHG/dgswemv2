#ifndef SWE_IC_SRC_FUNCTIONS_HPP
#define SWE_IC_SRC_FUNCTIONS_HPP

namespace SWE {
	//This sets up manufactured solution (hopefuly)

	double ic_ze(double t, Point<2>& pt) {
		double x1 = 40000.;
		double x2 = 148000.;
		double y1 = 10000.;
		double y2 = 53200.;

		double Ho = 2.;
		double zo = 0.25;

		double w = 2 * PI / 43200.;
		double tau = 0;

		return 2 * zo*cos(w*(pt[GlobalCoord::x] - x1)) * cos(w*(pt[GlobalCoord::y] - y1)) *
			cos(w*(t + tau)) / (cos(w*(x2 - x1)) * cos(w*(y2 - y1)));
	}

	double ic_qx(double t, Point<2>& pt) {
		double x1 = 40000.;
		double x2 = 148000.;
		double y1 = 10000.;
		double y2 = 53200.;

		double Ho = 2.;
		double zo = 0.25;

		double w = 2 * PI / 43200.;
		double tau = 0;

		return zo*sin(w*(pt[GlobalCoord::x] - x1)) * cos(w*(pt[GlobalCoord::y] - y1)) *
			sin(w*(t + tau)) / (cos(w*(x2 - x1)) * cos(w*(y2 - y1)));
	}
	
	double ic_qy(double t, Point<2>& pt) {
		double x1 = 40000.;
		double x2 = 148000.;
		double y1 = 10000.;
		double y2 = 53200.;

		double Ho = 2.;
		double zo = 0.25;

		double w = 2 * PI / 43200.;
		double tau = 0;

		return zo*cos(w*(pt[GlobalCoord::x] - x1)) * sin(w*(pt[GlobalCoord::y] - y1)) *
			sin(w*(t + tau)) / (cos(w*(x2 - x1)) * cos(w*(y2 - y1)));
	}

	double source_ze(double t, Point<2>& pt) {
		return 0;
	}

	double source_qx(double t, Point<2>& pt) {
		double x1 = 40000.;
		double x2 = 148000.;
		double y1 = 10000.;
		double y2 = 53200.;

		double Ho = 2.;
		double zo = 0.25;

		double w = 2 * PI / 43200.;
		double tau = 0;

		double x = pt[GlobalCoord::x];
		double y = pt[GlobalCoord::y];

		return w*zo*cos((t + tau)*w)*cos(w*(y - y1))*(1. / cos(w*(-x1 + x2)))*(1. / cos(w*(-y1 + y2)))*
			sin(w*(x - x1)) - 2 * SWE::Global::g*Ho*w*zo*cos((t + tau)*w)*cos(w*(y - y1))*
			(1. / cos(w*(-x1 + x2)))*(1. / cos(w*(-y1 + y2)))*sin(w*(x - x1)) -
			4.*SWE::Global::g*w*pow(zo, 2)*pow(cos((t + tau)*w), 2)*cos(w*(x - x1))*
			pow(cos(w*(y - y1)), 2)*pow((1. / cos(w*(-x1 + x2))), 2)*pow((1. / cos(w*(-y1 + y2))), 2)*
			sin(w*(x - x1)) + (3 * w*pow(zo, 2)*cos(w*(x - x1))*pow(cos(w*(y - y1)), 2)*
				pow((1. / cos(w*(-x1 + x2))), 2)*pow((1. / cos(w*(-y1 + y2))), 2)*
				pow(sin((t + tau)*w), 2)*sin(w*(x - x1))) /
				(Ho + 2 * zo*cos((t + tau)*w)*cos(w*(x - x1))*cos(w*(y - y1))*(1. / cos(w*(-x1 + x2)))*
			(1. / cos(w*(-y1 + y2)))) + (2 * w*pow(zo, 3)*cos((t + tau)*w)*
				pow(cos(w*(y - y1)), 3)*pow((1. / cos(w*(-x1 + x2))), 3)*
				pow((1. / cos(w*(-y1 + y2))), 3)*pow(sin((t + tau)*w), 2)*pow(sin(w*(x - x1)), 3))
			/ pow(Ho + 2 * zo*cos((t + tau)*w)*cos(w*(x - x1))*cos(w*(y - y1))*
			(1. / cos(w*(-x1 + x2)))*(1. / cos(w*(-y1 + y2))), 2) +
				(2 * w*pow(zo, 3)*cos((t + tau)*w)*pow(cos(w*(x - x1)), 2)*cos(w*(y - y1))*
					pow((1. / cos(w*(-x1 + x2))), 3)*pow((1. / cos(w*(-y1 + y2))), 3)*
					pow(sin((t + tau)*w), 2)*sin(w*(x - x1))*pow(sin(w*(y - y1)), 2)) /
			pow(Ho + 2 * zo*cos((t + tau)*w)*cos(w*(x - x1))*cos(w*(y - y1))*
			(1. / cos(w*(-x1 + x2)))*(1. / cos(w*(-y1 + y2))), 2) -
				(w*pow(zo, 2)*cos(w*(x - x1))*pow((1. / cos(w*(-x1 + x2))), 2)*
					pow((1. / cos(w*(-y1 + y2))), 2)*pow(sin((t + tau)*w), 2)*sin(w*(x - x1))*
					pow(sin(w*(y - y1)), 2)) /
					(Ho + 2 * zo*cos((t + tau)*w)*cos(w*(x - x1))*cos(w*(y - y1))*(1. / cos(w*(-x1 + x2)))*
			(1. / cos(w*(-y1 + y2))));
	}
	
	double source_qy(double t, Point<2>& pt) {
		double x1 = 40000.;
		double x2 = 148000.;
		double y1 = 10000.;
		double y2 = 53200.;

		double Ho = 2.;
		double zo = 0.25;

		double w = 2 * PI / 43200.;
		double tau = 0;

		double x = pt[GlobalCoord::x];
		double y = pt[GlobalCoord::y];
		
		return w*zo*cos((t + tau)*w)*cos(w*(x - x1))*(1. / cos(w*(-x1 + x2)))*(1. / cos(w*(-y1 + y2)))*
			sin(w*(y - y1)) - 2 * SWE::Global::g*Ho*w*zo*cos((t + tau)*w)*cos(w*(x - x1))*
			(1. / cos(w*(-x1 + x2)))*(1. / cos(w*(-y1 + y2)))*sin(w*(y - y1)) -
			4.*SWE::Global::g*w*pow(zo, 2)*pow(cos((t + tau)*w), 2)*pow(cos(w*(x - x1)), 2)*
			cos(w*(y - y1))*pow((1. / cos(w*(-x1 + x2))), 2)*pow((1. / cos(w*(-y1 + y2))), 2)*
			sin(w*(y - y1)) + (3 * w*pow(zo, 2)*pow(cos(w*(x - x1)), 2)*cos(w*(y - y1))*
				pow((1. / cos(w*(-x1 + x2))), 2)*pow((1. / cos(w*(-y1 + y2))), 2)*
				pow(sin((t + tau)*w), 2)*sin(w*(y - y1))) /
				(Ho + 2 * zo*cos((t + tau)*w)*cos(w*(x - x1))*cos(w*(y - y1))*(1. / cos(w*(-x1 + x2)))*
			(1. / cos(w*(-y1 + y2)))) + (2 * w*pow(zo, 3)*cos((t + tau)*w)*cos(w*(x - x1))*
				pow(cos(w*(y - y1)), 2)*pow((1. / cos(w*(-x1 + x2))), 3)*
				pow((1. / cos(w*(-y1 + y2))), 3)*pow(sin((t + tau)*w), 2)*pow(sin(w*(x - x1)), 2)*
				sin(w*(y - y1))) /
			pow(Ho + 2 * zo*cos((t + tau)*w)*cos(w*(x - x1))*cos(w*(y - y1))*
			(1. / cos(w*(-x1 + x2)))*(1. / cos(w*(-y1 + y2))), 2) -
				(w*pow(zo, 2)*cos(w*(y - y1))*pow((1. / cos(w*(-x1 + x2))), 2)*
					pow((1. / cos(w*(-y1 + y2))), 2)*pow(sin((t + tau)*w), 2)*pow(sin(w*(x - x1)), 2)*
					sin(w*(y - y1))) /
					(Ho + 2 * zo*cos((t + tau)*w)*cos(w*(x - x1))*cos(w*(y - y1))*(1. / cos(w*(-x1 + x2)))*
			(1. / cos(w*(-y1 + y2)))) + (2 * w*pow(zo, 3)*cos((t + tau)*w)*
				pow(cos(w*(x - x1)), 3)*pow((1. / cos(w*(-x1 + x2))), 3)*
				pow((1. / cos(w*(-y1 + y2))), 3)*pow(sin((t + tau)*w), 2)*pow(sin(w*(y - y1)), 3))
			/ pow(Ho + 2 * zo*cos((t + tau)*w)*cos(w*(x - x1))*cos(w*(y - y1))*
			(1. / cos(w*(-x1 + x2)))*(1. / cos(w*(-y1 + y2))), 2);
	}
}

#endif
