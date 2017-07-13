#ifndef DUBINER_HPP
#define DUBINER_HPP

#include <cassert>

#include "basis_base.hpp"

//These are two helper classes designed to assist in computing the value of the gradient
//at the derivative
struct DxSingularityHelper
{
  DxSingularityHelper(int q_max)
  {
    assert( q_max >= 0 );
    if ( q_max == 0 ) {

      a_mem = { 3. };
      eval_mem = { 1. };

    } else {

      DxSingularityHelper temp(q_max - 1);

      a_mem = std::move(temp.a_mem);
      eval_mem = std::move(temp.eval_mem);

      eval_mem.push_back( a_mem.back() + eval_mem.back() );

      double a_new = a_mem.back() + a_mem.size() + 2;
      a_mem.push_back(a_new);
    }
  }

  double get(int p, int q)
  {
    if ( p == 1 ) {
      return eval_mem[q];
    }

    return 0;
  }

  //memoizations of calls for q <= q_max
  std::vector<double> a_mem;
  std::vector<double> eval_mem;
};

struct DySingularityHelper
{
  DySingularityHelper(int q_max)
  {
    assert( q_max >= 0 );
    if ( q_max == 0 ) {

      a_mem = { 0.5 };
      eval_mem = { 0 };

    } else {

      DySingularityHelper temp(q_max - 1);

      a_mem = std::move(temp.a_mem);
      eval_mem = std::move(temp.eval_mem);

      eval_mem.push_back( a_mem.back() + eval_mem.back() );

      double a_new = a_mem.back() + 0.5* (a_mem.size() + 1);
      a_mem.push_back(a_new);

    }

  }

  double get(int p, int q)
  {
    if ( p > 1 ) {
      return 0;
    } else if ( p == 0 ) {
      return 3*eval_mem.at(q);
    }
    // p == 1
    return eval_mem.at(q+1);
  }
  //eval_mem stores the values corresponding to p = 1
  std::vector<double> a_mem;
  std::vector<double> eval_mem;
};


namespace Geometry {
namespace Basis {
//CRTP
class Dubiner : public Base<Dubiner,2>
{
public:
  Dubiner(int p) : _p(p), _grad_x(p), _grad_y(p+1) {}

  constexpr static int num_faces {3};

  using Base<Dubiner,2>::Point;

  std::vector<double> get_impl(const int id,
                               const std::vector<Point>& point);

  std::vector<std::array<double,2> > get_grad_impl(const int id,
                                                   const std::vector<Point>& point);

private:
  int _p;
  DxSingularityHelper _grad_x;
  DySingularityHelper _grad_y;

};
}
}
#endif
