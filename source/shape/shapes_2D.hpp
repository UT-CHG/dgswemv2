#ifndef SHAPES_2D_HPP
#define SHAPES_2D_HPP

#include "general_definitions.hpp"

namespace Shape {
class StraightTriangle : public Shape<2> {
  public:
    StraightTriangle() = default;
    StraightTriangle(std::vector<Point<3>>&& nodal_coordinates);

    std::vector<uint> GetBoundaryNodeID(const uint bound_id, const std::vector<uint> node_ID);

    Point<2> GetBarycentricCoordinates();
    std::vector<Point<2>> GetMidpointCoordinates();

    std::vector<double> GetJdet(const std::vector<Point<2>>& points);
    Array3D<double> GetJinv(const std::vector<Point<2>>& points);
    std::vector<double> GetSurfaceJ(const uint bound_id, const std::vector<Point<2>>& points);
    Array2D<double> GetSurfaceNormal(const uint bound_id, const std::vector<Point<2>>& points);

    Array2D<double> GetPsi(const std::vector<Point<2>>& points);
    Array3D<double> GetDPsi(const std::vector<Point<2>>& points);

    Array2D<double> GetBoundaryPsi(const uint bound_id, const std::vector<Point<1>>& points);

    std::vector<Point<2>> LocalToGlobalCoordinates(const std::vector<Point<2>>& points);

    void GetVTK(std::vector<Point<3>>& points, Array2D<uint>& cells);

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & hpx::serialization::base_object<Shape<2>>(*this);
        // clang-format on
    }
    HPX_SERIALIZATION_POLYMORPHIC_TEMPLATE(StraightTriangle);
#endif
};
}
#endif
