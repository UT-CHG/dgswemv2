#ifndef HPX_MESH_HPP
#define HPX_MESH_HPP

template <typename ProblemType>
class hpx_mesh : public hpx::components::simple_component_base<hpx_mesh<ProblemType>> {
  public:
    typename ProblemType::mesh_type mesh;

    hpx_mesh() = default;  // default constructor to register hpx component
    hpx_mesh(uint p, const MeshMetaData& mesh_data) : mesh(p, mesh_data._mesh_name) {
        initialize_mesh<ProblemType>(this->mesh, mesh_data);
    };
};

#endif