#ifndef ABSTRACT_LOAD_BALANCER_FACTORY_HPP
#define ABSTRACT_LOAD_BALANCER_FACTORY_HPP

#include "base_model.hpp"

namespace LoadBalancer {
class AbstractFactory {
public:
    template <typename ProblemType>
    static hpx::future<void> initialize_locality_and_world_models(const uint locality_id,
                                                                  const std::string& input_string);

    template <typename ProblemType>
    static void reset_locality_and_world_models();

    template <typename ProblemType>
    static std::unique_ptr<SubmeshModel> create_submesh_model(uint locality_id, uint submesh_id);
};
}
#endif