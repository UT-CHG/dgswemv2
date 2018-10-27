#ifndef INITIALIZE_MESH_SKELETON_HPP
#define INITIALIZE_MESH_SKELETON_HPP

template <typename ProblemType>
void initialize_mesh_skeleton(typename ProblemType::ProblemMeshType& mesh,
                              typename ProblemType::ProblemMeshSkeletonType& mesh_skeleton,
                              typename ProblemType::ProblemWriterType& writer) {
    ProblemType::create_edge_boundaries(mesh, mesh_skeleton, writer);
    ProblemType::create_edge_interfaces(mesh, mesh_skeleton, writer);
    ProblemType::create_edge_distributeds(mesh, mesh_skeleton, writer);
}

#endif