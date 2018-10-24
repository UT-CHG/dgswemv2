#include "shape/shapes_2D.hpp"
#include "utilities/almost_equal.hpp"

bool is_equal(Shape::StraightTriangle& o_tri, Shape::StraightTriangle& i_tri) {
    using Utilities::almost_equal;

    std::vector<Point<2>> test_points{{0.5, 0.5}, {-0.5, 0.5}, {0, 0}};

    bool is_equal{true};

    for (uint bd_id = 0; bd_id < 3; ++bd_id) {
        std::vector<Point<2>> test_pt{test_points[bd_id]};

	using Vec = StatVector<double,2>;
        std::vector<Vec,AlignedAllocator<Vec>> o_normal = o_tri.GetSurfaceNormal(bd_id, test_pt);
        std::vector<Vec,AlignedAllocator<Vec>> i_normal = i_tri.GetSurfaceNormal(bd_id, test_pt);

        is_equal &= almost_equal(o_normal[0][0], i_normal[0][0]);
        is_equal &= almost_equal(o_normal[0][1], i_normal[0][1]);
    }

    return is_equal;
}

int main() {
    std::vector<Point<3>> nodal_coord(3);
    nodal_coord[0] = {-1, 0, 0};
    nodal_coord[1] = {1, 0, 0};
    nodal_coord[2] = {0, 1, 0};

    Shape::StraightTriangle o_tri(std::move(nodal_coord));
    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_tri;

    hpx::serialization::input_archive i_archive(buffer);
    Shape::StraightTriangle i_tri;
    i_archive >> i_tri;

    if (!is_equal(o_tri, i_tri)) {
        return 1;
    }
    return 0;
}
