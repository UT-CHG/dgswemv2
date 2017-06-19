#ifndef SWE_BOUNDARY_CONDITIONS_HPP
#define SWE_BOUNDARY_CONDITIONS_HPP

#include <array>
#include <vector>

namespace SWE {
enum class BoundaryConditions
{
  _land,
  _tidal
};

class BoundaryBase
{
public:
  BoundaryBase(uint n_gp)
    : ze_at_gp(n_gp), qx_at_gp(n_gp), qy_at_gp(n_gp), bath_at_gp(n_gp),
      ze_numerical_flux_at_gp(n_gp), qx_numerical_flux_at_gp(n_gp),
      qy_numerical_flux_at_gp(n_gp)
  {}

  virtual void get_ex(const double t, const double ze_in,const double qx_in,
                      const double qy_in, const double bath,
                      const std::array<double,2>& normal,
                      double& ze_ex,  double& qx_ex, double& qy_ex
                      ) = 0;

  std::vector<double> ze_at_gp;
  std::vector<double> qx_at_gp;
  std::vector<double> qy_at_gp;
  std::vector<double> bath_at_gp;

  std::vector<double> ze_numerical_flux_at_gp;
  std::vector<double> qx_numerical_flux_at_gp;
  std::vector<double> qy_numerical_flux_at_gp;
};

class LandBoundary final : public BoundaryBase
{
public:
  LandBoundary(uint n_gp)
    : BoundaryBase(n_gp)
  {}

  void get_ex(const double t, const double ze_in,const double qx_in,const double qy_in,
              const double bath,
              const std::array<double,2>& normal,
              double& ze_ex,  double& qx_ex, double& qy_ex)
  {
    ze_ex = ze_in;

    double qn = normal[0]*qx_in + normal[1]*qy_in; //normal component of momentum
    double qt =-normal[1]*qx_in + normal[0]*qy_in; //tangential component of momentum

    //observe that the matrix created is orthogonal so it can be inverted by taking the transpose
    qx_ex = normal[0]*-qn - normal[1]*qt;
    qy_ex = normal[1]*-qn + normal[0]*qt;
  }

};

class TidalBoundary final : public BoundaryBase
{
public:
  TidalBoundary(uint n_gp, double amplitude, std::vector<double>& frequency, double phase,
                double t_ramp)
    : BoundaryBase(n_gp), _amplitude(amplitude), _frequency(std::move(frequency)),
      _phase(phase), _t_ramp(t_ramp)
  {}

  void get_ex(const double t, const double ze_in,const double qx_in,const double qy_in,
              const double bath,
              const std::array<double,2>& normal,
              double& ze_ex,  double& qx_ex, double& qy_ex)
  {
    double ramp_frac = t/_t_ramp;

    double frequency_eval = 0;
    for ( auto& frequency_component : _frequency ) {
      frequency_eval += std::sin( frequency_component * (t + _phase) );
    }

    ze_ex = std::min(1., ramp_frac)* _amplitude* frequency_eval;

    qx_ex = qx_in;
    qy_ex = qy_in;
  }

private:
  const double _amplitude;
  const std::vector<double> _frequency;
  const double _phase;
  const double _t_ramp;
};
}
#endif
