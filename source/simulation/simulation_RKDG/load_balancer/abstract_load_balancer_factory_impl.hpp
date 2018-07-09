#ifndef ABSTRACT_LOAD_BALANCER_FACTORY_IMPL_HPP
#define ABSTRACT_LOAD_BALANCER_FACTORY_IMPL_HPP

#include "abstract_load_balancer_factory.hpp"
#include "random.hpp"

namespace RKDG {
namespace LoadBalancer {
template <typename ProblemType>
hpx::future<void> AbstractFactory::initialize_locality_and_world_models(const uint locality_id,
                                                                        const std::string& input_string) {
    return Random<ProblemType>::initialize_locality_and_world_models(locality_id, input_string);
}

template <typename ProblemType>
void AbstractFactory::reset_locality_and_world_models() {
    return Random<ProblemType>::reset_locality_and_world_models();
}

template <typename ProblemType>
std::unique_ptr<SubmeshModel> AbstractFactory::create_submesh_model(uint locality_id,
                                                                    uint submesh_id,
                                                                    const LoadBalancerInput& load_balancer_input) {
    if (load_balancer_input.use_load_balancer) {
        if (load_balancer_input.name == "random") {
            return Random<ProblemType>::create_submesh_model(
                locality_id, submesh_id, load_balancer_input.rebalance_frequency);
        } else {
            std::string err_msg{"Error: Unknown RKDG load balancer of type: " + load_balancer_input.name};
            throw std::logic_error(err_msg);
        }
    } else {
        return nullptr;
    }
}
}
}
#endif