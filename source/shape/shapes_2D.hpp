#ifndef SHAPES_2D_HPP
#define SHAPES_2D_HPP

#include "general_definitions.hpp"

namespace Shape {
class StraightTriangle : public Shape<2> {
  public:
    StraightTriangle() = default;
    StraightTriangle(std::vector<Point<3>>&& nodal_coordinates);

    std::vector<uint> GetBoundaryNodeID(const uint bound_id, const std::vector<uint>& node_ID) const;

    Point<2> GetBarycentricCoordinates() const;
    std::vector<Point<2>> GetMidpointCoordinates() const;

    DynVector<double> GetJdet(const std::vector<Point<2>>& points) const;
    DynVector<double> GetSurfaceJ(const uint bound_id, const std::vector<Point<2>>& points) const;
    AlignedVector<StatMatrix<double, 2, 2>> GetJinv(const std::vector<Point<2>>& points) const;
    AlignedVector<StatVector<double, 2>> GetSurfaceNormal(const uint bound_id,
                                                          const std::vector<Point<2>>& points) const;

    DynMatrix<double> GetPsi(const std::vector<Point<2>>& points) const;
    std::array<DynMatrix<double>, 2> GetDPsi(const std::vector<Point<2>>& points) const;
    DynMatrix<double> GetBoundaryPsi(const uint bound_id, const std::vector<Point<1>>& points) const;

    std::vector<Point<2>> LocalToGlobalCoordinates(const std::vector<Point<2>>& points) const;
    std::vector<Point<2>> GlobalToLocalCoordinates(const std::vector<Point<2>>& points) const;
    bool ContainsPoint(const Point<2>& point) const;

    void GetVTK(std::vector<Point<3>>& points, Array2D<uint>& cells) const;

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
