#ifndef ABSTRACT_LOAD_BALANCER_FACTORY_IMPL_HPP
#define ABSTRACT_LOAD_BALANCER_FACTORY_IMPL_HPP

#include "abstract_load_balancer_factory.hpp"
#include "random.hpp"


namespace LoadBalancer {
template <typename ProblemType>
hpx::future<void> AbstractFactory::initialize_locality_and_world_models(const uint locality_id) {
    return Random<ProblemType>::initialize_locality_and_world_models(locality_id);
}

template <typename ProblemType>
void AbstractFactory::reset_locality_and_world_models() {
    return Random<ProblemType>::reset_locality_and_world_models();
}

template <typename ProblemType>
std::unique_ptr<SubmeshModel> AbstractFactory::create_submesh_model(uint locality_id, uint submesh_id) {
    return Random<ProblemType>::create_submesh_model(locality_id, submesh_id);
}
}
#endif