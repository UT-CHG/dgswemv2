#ifndef SWE_PRE_CREATE_BOUND_HPP
#define SWE_PRE_CREATE_BOUND_HPP

namespace SWE {
template <typename RawBoundaryType>
void Problem::create_boundaries_kernel(ProblemMeshType& mesh,
                                       std::map<uchar, std::vector<RawBoundaryType>>& pre_boundaries,
                                       Writer<SWE::Problem>& writer) {
    uint n_bound_old_land = 0;
    uint n_bound_old_tidal = 0;
    uint n_bound_old_flow = 0;

    using BoundaryTypes = Geometry::BoundaryTypeTuple<SWE::Data, SWE::Land, SWE::Tidal, SWE::Flow>;

    for (auto it = pre_boundaries.begin(); it != pre_boundaries.end(); it++) {
        switch (it->first) {
            case SWE::land:
                using BoundaryTypeLand = typename std::tuple_element<0, BoundaryTypes>::type;

                n_bound_old_land = mesh.GetNumberBoundaries();

                for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
                    mesh.template CreateBoundary<BoundaryTypeLand>(*itt);
                }

                if (writer.WritingLog()) {
                    writer.GetLogFile() << "Number of land boundaries: " << mesh.GetNumberBoundaries() -
                                                                                n_bound_old_land << std::endl;
                }

                break;
            case SWE::tidal:
                using BoundaryTypeTidal = typename std::tuple_element<1, BoundaryTypes>::type;

                n_bound_old_tidal = mesh.GetNumberBoundaries();

                for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
                    mesh.template CreateBoundary<BoundaryTypeTidal>(*itt);
                }

                if (writer.WritingLog()) {
                    writer.GetLogFile() << "Number of tidal boundaries: " << mesh.GetNumberBoundaries() -
                                                                                 n_bound_old_tidal << std::endl;
                }

                break;
            case SWE::flow:
                using BoundaryTypeFlow = typename std::tuple_element<2, BoundaryTypes>::type;

                n_bound_old_flow = mesh.GetNumberBoundaries();

                for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
                    mesh.template CreateBoundary<BoundaryTypeFlow>(*itt);
                }

                if (writer.WritingLog()) {
                    writer.GetLogFile() << "Number of flow boundaries: " << mesh.GetNumberBoundaries() -
                                                                                n_bound_old_flow << std::endl;
                }

                break;
        }
    }
}
}

#endif