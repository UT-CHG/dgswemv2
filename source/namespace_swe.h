#ifndef SWE_H
#define SWE_H
#include<cmath>
//todo: remove global variable

namespace SWE {
  const double g =9.81;

  auto F00 = [](double ze, double qx, double qy, double b)->double
  {
    return qx;
  };

  auto F01 = [](double ze, double qx, double qy, double b)->double
  {
    return qx*qx/(ze+b) + g *(0.5*ze*ze + ze* b);
  };

  auto F02 = [](double ze, double qx, double qy, double b)->double
  {
    return qy*qx/(ze+b);
  };

  auto F10 = [](double ze, double qx, double qy, double b)->double
  {
    return qy;
  };

  auto F11 = [](double ze, double qx, double qy, double b)->double
  {
    return qy*qx/(ze+b);
  };

  auto F12 = [](double ze, double qx, double qy, double b)->double
  {
    return qy*qy/(ze+b) + g * (0.5*ze*ze + ze*b);
  };

  auto S0 = [](double ze, double dbdx)->double
  {
    return 0;
  };

  auto S1 = [](double ze, double dbdx)->double
  {
    return g*ze*dbdx;
  };

  auto S2 = [](double ze, double dbdy)->double
  {
    return g*ze*dbdy;
  };

  auto eigmax = [](double ze, double qx, double qy, double b)->double
  {
    return std::sqrt(g*(ze+b)) + std::hypot(qx/(ze+b),qy/(ze+b));
  };
}
#endif
