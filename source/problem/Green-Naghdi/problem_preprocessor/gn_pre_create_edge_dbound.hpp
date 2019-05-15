#ifndef EHDG_GN_PRE_CREATE_EDGE_DBOUND_HPP
#define EHDG_GN_PRE_CREATE_EDGE_DBOUND_HPP

namespace GN {
template <typename ProblemType>
void create_edge_distributeds(typename ProblemType::ProblemMeshType& mesh,
                              typename ProblemType::ProblemMeshSkeletonType& mesh_skeleton,
                              typename ProblemType::ProblemWriterType& writer) {
    using DBTypeDistributed =
        typename std::tuple_element<0, typename ProblemType::ProblemDistributedBoundaryTypes>::type;
    using DBTypeDistributedLevee =
        typename std::tuple_element<1, typename ProblemType::ProblemDistributedBoundaryTypes>::type;

    using EDTypeDistributed = typename std::tuple_element<0, typename ProblemType::ProblemEdgeDistributedTypes>::type;
    using EDTypeDistributedLevee =
        typename std::tuple_element<1, typename ProblemType::ProblemEdgeDistributedTypes>::type;

    uint old_edge_distributed = mesh_skeleton.GetNumberEdgeDistributeds();

    mesh.template CallForEachDistributedBoundaryOfType<DBTypeDistributed>([&mesh_skeleton](auto& dbound) {
        uint locality_in = dbound.boundary_condition.exchanger.locality_in;
        uint submesh_in  = dbound.boundary_condition.exchanger.submesh_in;

        uint locality_ex = dbound.boundary_condition.exchanger.locality_ex;
        uint submesh_ex  = dbound.boundary_condition.exchanger.submesh_ex;

        bool ccw;

        if (locality_in < locality_ex || (locality_in == locality_ex && submesh_in < submesh_ex)) {
            ccw = true;  // trace basis is defined counterclockwise around element
        } else if (locality_in > locality_ex || (locality_in == locality_ex && submesh_in > submesh_ex)) {
            ccw = false;  // trace basis is defined clockwise around element
        }

        mesh_skeleton.template CreateEdgeDistributed<EDTypeDistributed>(dbound, ccw);
    });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of distributed edges: "
                            << mesh_skeleton.GetNumberEdgeDistributeds() - old_edge_distributed << std::endl;
    }

    uint old_edge_distributed_levee = mesh_skeleton.GetNumberEdgeDistributeds();

    mesh.template CallForEachDistributedBoundaryOfType<DBTypeDistributedLevee>([&mesh_skeleton](auto& dbound) {
        uint locality_in = dbound.boundary_condition.exchanger.locality_in;
        uint submesh_in  = dbound.boundary_condition.exchanger.submesh_in;

        uint locality_ex = dbound.boundary_condition.exchanger.locality_ex;
        uint submesh_ex  = dbound.boundary_condition.exchanger.submesh_ex;

        bool ccw;

        if (locality_in < locality_ex || (locality_in == locality_ex && submesh_in < submesh_ex)) {
            ccw = true;  // trace basis is defined counterclockwise around element
        } else if (locality_in > locality_ex || (locality_in == locality_ex && submesh_in > submesh_ex)) {
            ccw = false;  // trace basis is defined clockwise around element
        }

        mesh_skeleton.template CreateEdgeDistributed<EDTypeDistributedLevee>(dbound, ccw);
    });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of distributed levee edges: "
                            << mesh_skeleton.GetNumberEdgeDistributeds() - old_edge_distributed_levee << std::endl;
    }
}
}

#endif
