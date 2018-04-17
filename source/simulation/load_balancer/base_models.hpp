#ifndef LOAD_BALANCER_BASE_MODELS_HPP
#define LOAD_BALANCER_BASE_MODELS_HPP

namespace LoadBalancer {
class SubmeshModel {
  public:
    SubmeshModel(uint locality_id, uint submesh_id)
        : locality_id(locality_id), submesh_id(submesh_id) {}

    virtual void InStep(uint64_t compute_cost, uint64_t memory_cost) {};

  protected:
    uint locality_id, submesh_id;
};
}
#endif