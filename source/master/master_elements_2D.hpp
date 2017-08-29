#ifndef CLASS_MASTER_ELEMENT_HPP
#define CLASS_MASTER_ELEMENT_HPP

#include "../general_definitions.hpp"

namespace Master {
template <class BasisType, class IntegrationType>
class Triangle : public Master<2> {
  public:
    BasisType basis;
    IntegrationType integration;

  public:
    Triangle(uint);

    std::vector<Point<2>> BoundaryToMasterCoordinates(uint, const std::vector<Point<1>>&);

  private:
    std::vector<Point<2>> VTKPostCell();
    std::vector<Point<2>> VTKPostPoint();
};
}

#include "elements_2D/master_triangle.tpp"

#endif