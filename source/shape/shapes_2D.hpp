#ifndef SHAPES_2D_HPP
#define SHAPES_2D_HPP

#include "../general_definitions.hpp"

namespace Shape {
class StraightTriangle : public Shape<2> {
  public:
    StraightTriangle()=default;
    StraightTriangle(const std::vector<Point<2>>& nodal_coordinates) : Shape<2>(nodal_coordinates) {}

    bool CheckJacobianPositive(const Point<2>& point);

    Point<2> GetBarycentricCoordinates();
    std::vector<Point<2>> GetMidpointCoordinates();

    std::vector<double> GetJdet(const std::vector<Point<2>>& points);
    Array3D<double> GetJinv(const std::vector<Point<2>>& points);
    std::vector<double> GetSurfaceJ(const uint bound_id, const std::vector<Point<2>>& points);
    Array2D<double> GetSurfaceNormal(const uint bound_id, const std::vector<Point<2>>& points);

    std::vector<double> InterpolateNodalValues(const std::vector<double>& nodal_values,
                                               const std::vector<Point<2>>& points);
    std::vector<Point<2>> LocalToGlobalCoordinates(const std::vector<Point<2>>& points);

    void GetVTK(std::vector<Point<3>>& points, Array2D<uint>& cells);

#ifdef HAS_HPX
  public:
    template <typename Archive>
    void serialize(Archive& ar, unsigned);
    HPX_SERIALIZATION_POLYMORPHIC_TEMPLATE(StraightTriangle);
#endif
};

#ifdef HAS_HPX
template <typename Archive>
void StraightTriangle::serialize(Archive& ar, unsigned) {
    ar & hpx::serialization::base_object<Shape<2>>(*this);
}
#endif
}
#endif
