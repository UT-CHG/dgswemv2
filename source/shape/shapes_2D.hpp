#ifndef SHAPES_2D_HPP
#define SHAPES_2D_HPP

#include "general_definitions.hpp"

namespace Shape {
class StraightTriangle : public Shape<2> {
  public:
    StraightTriangle() = default;
    StraightTriangle(DynVector<Point<3>>&& nodal_coordinates);

    DynVector<uint> GetBoundaryNodeID(const uint bound_id, const DynVector<uint> node_ID);

    Point<2> GetBarycentricCoordinates();
    DynVector<Point<2>> GetMidpointCoordinates();

    DynVector<double> GetJdet(const DynVector<Point<2>>& points);
    DynVector<StatMatrix<double, 2, 2>> GetJinv(const DynVector<Point<2>>& points);
    DynVector<double> GetSurfaceJ(const uint bound_id, const DynVector<Point<2>>& points);
    DynVector<StatVector<double, 2>> GetSurfaceNormal(const uint bound_id, const DynVector<Point<2>>& points);

    DynMatrix<double> GetPsi(const DynVector<Point<2>>& points);
    StatVector<DynMatrix<double>, 2> GetDPsi(const DynVector<Point<2>>& points);

    DynMatrix<double> GetBoundaryPsi(const uint bound_id, const DynVector<Point<1>>& points);

    DynVector<Point<2>> LocalToGlobalCoordinates(const DynVector<Point<2>>& points);

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
