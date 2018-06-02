#ifndef LOAD_BALANCER_RANDOM_HPP
#define LOAD_BALANCER_RANDOM_HPP
#include "../../utilities/heartbeat.hpp"
#include "base_model.hpp"
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

    WorldModel(const std::string& input_string);
    void MigrateOneSubmesh();
    HPX_DEFINE_COMPONENT_ACTION(WorldModel, MigrateOneSubmesh, MigrateOneSubmeshAction);

  private:
    bool tried_moving_one_tile = false;
    using client_locality_id_pair = std::pair<client_t,uint>;
    std::vector<client_locality_id_pair> simulation_unit_clients;
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
        using ActionType = typename WorldModel<ProblemType>::MigrateOneSubmeshAction;
        hpx::apply<ActionType>(this->get_id());
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename ProblemType>
class SubmeshModel : public LoadBalancer::SubmeshModel {
  public:
    using BaseType = LoadBalancer::SubmeshModel;

    SubmeshModel() = default;
    SubmeshModel(const std::chrono::duration<double>& rebalance_period,
                 uint locality_id, uint submesh_id);

    void InStep(uint64_t compute_cost, uint64_t memory_cost);

    template <typename Archive>
    void serialize(Archive& ar, unsigned);
    HPX_SERIALIZATION_POLYMORPHIC(SubmeshModel<ProblemType>);

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
    static hpx::future<void> initialize_locality_and_world_models(const uint locality_id,
                                                                  const std::string& input_string);
    static void reset_locality_and_world_models();
    static std::unique_ptr<LoadBalancer::SubmeshModel> create_submesh_model(uint locality_id, uint submesh_id);
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
WorldModel<ProblemType>::WorldModel(const std::string& input_string) {
    InputParameters<typename ProblemType::ProblemInputType> input(input_string);

    for ( uint locality_id = 0; locality_id < hpx::get_initial_num_localities(); ++locality_id ) {
        std::string submesh_file_prefix = input.mesh_file_name.substr(0, input.mesh_file_name.find_last_of('.')) + "_" +
                                          std::to_string(locality_id) + '_';
        std::string submesh_file_postfix =
            input.mesh_file_name.substr(input.mesh_file_name.find_last_of('.'), input.mesh_file_name.size());

        uint submesh_id = 0;
        while ( Utilities::file_exists(submesh_file_prefix + std::to_string(submesh_id) + submesh_file_postfix) ) {
            this->simulation_unit_clients.push_back(
                std::make_pair( client_t(), locality_id )
                );

            this->simulation_unit_clients.back().first.connect_to(std::string{client_t::GetBasename()}+
                                                            std::to_string(locality_id)+'_'+
                                                            std::to_string(submesh_id));
            std::cout << "Adding client for locality: " << locality_id << " ad submesh: " << submesh_id << '\n';
            ++submesh_id;
        }
    }
}

template <typename ProblemType>
void WorldModel<ProblemType>::MigrateOneSubmesh() {
    if ( !this->tried_moving_one_tile ) {
        const std::vector<hpx::naming::id_type> localities = hpx::find_all_localities();

        client_locality_id_pair& curr_target = this->simulation_unit_clients[ 0xdeadbeef % this->simulation_unit_clients.size() ];
        uint target_locality = ( curr_target.second + 1 ) % localities.size();

        std::cout << "Stealing tile from " << curr_target.second << " and moving it to " << target_locality << std::endl;

        curr_target.first = hpx::components::migrate(curr_target.first,localities[target_locality]);
        std::cout << "Made it past steal" << std::endl;
        curr_target.second = target_locality;
        this->tried_moving_one_tile = true;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SubmeshModel Implementation
template <typename ProblemType>
SubmeshModel<ProblemType>::SubmeshModel(const std::chrono::duration<double>& rebalance_period,
                           uint locality_id, uint submesh_id)
    : BaseType(locality_id, submesh_id), beat(rebalance_period) {}

template <typename ProblemType>
void SubmeshModel<ProblemType>::InStep(uint64_t, uint64_t) {
    /*if ( this->locality_id == 0 && this->submesh_id == 0 && this->beat.Thump() ) {
        std::cout << "In submeshmode->Instep()\n";
        assert(Random<ProblemType>::world_model_client);
        Random<ProblemType>::world_model_client.MigrateOneSubmesh();
        }*/
}

template <typename ProblemType>
template <typename Archive>
void SubmeshModel<ProblemType>::serialize(Archive& ar, unsigned) {
    ar & hpx::serialization::base_object<LoadBalancer::SubmeshModel>(*this);
    ar & beat;
}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename ProblemType>
hpx::future<void> Random<ProblemType>::initialize_locality_and_world_models(const uint locality_id,
                                                                            const std::string& input_string) {
    if ( locality_id == 0 ) {
        std::cout << "Initializing The world model\n";

        using WorldModelComponent = hpx::components::simple_component<Random<ProblemType>::WorldModel>;
        hpx::future<hpx::id_type> world_model_id = hpx::new_<WorldModelComponent>(hpx::find_here(), input_string);

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

template <typename ProblemType>
std::unique_ptr<LoadBalancer::SubmeshModel> Random<ProblemType>::create_submesh_model(uint locality_id, uint submesh_id) {
    return std::make_unique<Random<ProblemType>::SubmeshModel>(std::chrono::duration<double>(2.), locality_id, submesh_id);
}
}
#endif