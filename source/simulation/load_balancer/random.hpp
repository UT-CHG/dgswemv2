#ifndef LOAD_BALANCER_RANDOM_HPP
#define LOAD_BALANCER_RANDOM_HPP
#include "../../utilities/heartbeat.hpp"
#include "base_models.hpp"
#include "../hpx_simulation_unit.hpp"

#include <cstdlib>

namespace LoadBalancer {
namespace detail_random {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename ProblemType>
class WorldModel : public hpx::components::simple_component_base<WorldModel<ProblemType>> {
  public:
    using client_t = HPXSimulationUnitClient<ProblemType>;

    static constexpr const char* GetBasename() { return "LOAD_BALANCER_RANDOMIZED_WORLD_MODEL"; }

    WorldModel(/*meshrelated*/);
    void MigrateOneSubmesh();
    HPX_DEFINE_COMPONENT_ACTION(WorldModel, MigrateOneSubmesh, MigrateOneSubmeshAction);

  private:
    std::vector<client_t> simulation_unit_clients;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename ProblemType>
class WorldModelClient final :
        public hpx::components::client_base<WorldModelClient<ProblemType>, WorldModel<ProblemType>> {
  private:
    using BaseType = hpx::components::client_base<WorldModelClient<ProblemType>, WorldModel<ProblemType>>;

  public:
    WorldModelClient()=default;
    WorldModelClient(hpx::future<hpx::id_type>&&id) : BaseType(std::move(id)) {}

    void MigrateOneSubmesh() {
        using ActionType = typename WorldModel<ProblemType>::MigrateOneSubmesh;
        hpx::apply<ActionType>(this->get_id());
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename ProblemType>
class SubmeshModel final : public LoadBalancer::SubmeshModel {
  public:
    using BaseType = LoadBalancer::SubmeshModel;

    SubmeshModel(const std::chrono::duration<double>& rebalance_period,
                 uint locality_id, uint submesh_id);
    void InStep(uint64_t compute_cost, uint64_t memory_cost);

  private:
    Utilities::HeartBeat beat;
};
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename ProblemType>
struct Random {
    using SubmeshModel = detail_random::SubmeshModel<ProblemType>;
    using WorldModel = detail_random::WorldModel<ProblemType>;
    using WorldModelClient = detail_random::WorldModelClient<ProblemType>;

    static WorldModelClient world_model_client;
    static hpx::future<void> initialize_locality_and_world_models(const uint locality_id);
    static void reset_locality_and_world_models();
};
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define DGSWEMV2_REGISTER_RANDOM_LOAD_BALANCER(ProblemType)                                                           \
    using rndm_lb_world_client_ = LoadBalancer::Random<ProblemType>::WorldModelClient;                                \
    template<> rndm_lb_world_client_ LoadBalancer::Random<ProblemType>::world_model_client = rndm_lb_world_client_(); \
    using rndm_lb_world_model_ = LoadBalancer::Random<ProblemType>::WorldModel;                                       \
    using rndm_lb_world_model_component_ = hpx::components::simple_component<rndm_lb_world_model_>;                   \
    HPX_REGISTER_COMPONENT(rndm_lb_world_model_component_,rndm_lb_world_model_);
/**/
namespace LoadBalancer { namespace detail_random {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// WorldModel Implementation
template <typename ProblemType>
WorldModel<ProblemType>::WorldModel() {

}

template <typename ProblemType>
void WorldModel<ProblemType>::MigrateOneSubmesh() {
    std::cout << "Trying to steal one tile" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SubmeshModel Implementation
template <typename ProblemType>
SubmeshModel<ProblemType>::SubmeshModel(const std::chrono::duration<double>& rebalance_period,
                           uint locality_id, uint submesh_id)
    : BaseType(locality_id, submesh_id), beat(rebalance_period) {}

template <typename ProblemType>
void SubmeshModel<ProblemType>::InStep(uint64_t, uint64_t) {
    if ( locality_id == 0 && submesh_id == 0 && beat.Thump() ) {
        Random<ProblemType>::world_model_client.MigrateOneSubmesh();
    }
}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename ProblemType>
hpx::future<void> Random<ProblemType>::initialize_locality_and_world_models(const uint locality_id) {
    if ( locality_id == 0 ) {
        std::cout << "Initializing The world model\n";

        using WorldModelComponent = hpx::components::simple_component<Random<ProblemType>::WorldModel>;
        hpx::future<hpx::id_type> world_model_id = hpx::new_<WorldModelComponent>(hpx::find_here());

        Random<ProblemType>::world_model_client = Random<ProblemType>::WorldModelClient(std::move(world_model_id));

        hpx::future<void> registration_future = Random<ProblemType>::world_model_client.register_as(
            Random<ProblemType>::WorldModel::GetBasename()
            );
        std::cout << "All calls made, just need to wait on future\n";
        return registration_future;
    }
    return hpx::make_ready_future();

}

template <typename ProblemType>
void Random<ProblemType>::reset_locality_and_world_models() {
    Random<ProblemType>::world_model_client = WorldModelClient();

}
}
#endif