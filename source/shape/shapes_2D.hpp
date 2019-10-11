#ifndef SHAPES_2D_HPP
#define SHAPES_2D_HPP

#include "general_definitions.hpp"

namespace Shape {
class StraightTriangle : public Shape<2> {
  public:
    StraightTriangle() = default;
    StraightTriangle(AlignedVector<Point<3>>&& nodal_coordinates);

    std::vector<uint> GetBoundaryNodeID(const uint bound_id, const std::vector<uint>& node_ID) const override;

    double GetArea() const override;
    Point<2> GetBarycentricCoordinates() const override;
    AlignedVector<Point<2>> GetMidpointCoordinates() const override;

    DynVector<double> GetJdet(const AlignedVector<Point<2>>& points) const override;
    DynVector<double> GetSurfaceJ(const uint bound_id, const AlignedVector<Point<2>>& points) const override;
    AlignedVector<StatMatrix<double, 2, 2>> GetJinv(const AlignedVector<Point<2>>& points) const override;
    AlignedVector<StatVector<double, 2>> GetSurfaceNormal(const uint bound_id,
                                                          const AlignedVector<Point<2>>& points) const override;

    DynMatrix<double> GetPsi(const AlignedVector<Point<2>>& points) const override;
    std::array<DynMatrix<double>, 2> GetDPsi(const AlignedVector<Point<2>>& points) const override;
    DynMatrix<double> GetBoundaryPsi(const uint bound_id, const AlignedVector<Point<1>>& points) const override;

    AlignedVector<Point<2>> LocalToGlobalCoordinates(const AlignedVector<Point<2>>& points) const override;
    AlignedVector<Point<2>> GlobalToLocalCoordinates(const AlignedVector<Point<2>>& points) const override;
    bool ContainsPoint(const Point<2>& point) const override;

    void GetVTK(AlignedVector<Point<3>>& points, Array2D<uint>& cells) const override;

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
