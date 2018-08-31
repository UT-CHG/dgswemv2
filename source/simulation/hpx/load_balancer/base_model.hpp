#ifndef LOAD_BALANCER_BASE_MODEL_HPP
#define LOAD_BALANCER_BASE_MODEL_HPP

namespace LoadBalancer {
class SubmeshModel {
  public:
    SubmeshModel() = default;
    SubmeshModel(uint locality_id, uint submesh_id) : locality_id(locality_id), submesh_id(submesh_id) {}

    virtual void InStep(uint64_t compute_cost, uint64_t memory_cost) { std::cout << "Firing base_model::Instep()\n"; }

    virtual ~SubmeshModel() = default;

    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & locality_id
            & submesh_id;
        // clang-format on
    }
    HPX_SERIALIZATION_POLYMORPHIC(SubmeshModel);

    virtual std::pair<uint, uint> get_tag() { return std::make_pair(locality_id, submesh_id); }

  protected:
    uint locality_id, submesh_id;
};
}

#endif