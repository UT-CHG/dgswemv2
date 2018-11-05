#ifndef SIM_UNIT_HPX_BASE_HPP
#define SIM_UNIT_HPX_BASE_HPP

#include "general_definitions.hpp"
#include "utilities/is_defined.hpp"

struct HPXSimulationUnitBase
//    : public hpx::components::migration_support<hpx::components::component_base<HPXSimulationUnit<ProblemType>>> {
    : public hpx::components::abstract_managed_component_base<HPXSimulationUnitBase> {

    virtual ~HPXSimulationUnitBase() = default;

    virtual hpx::future<void> Preprocessor() = 0;
    hpx::future<void> Preprocessor_() { return Preprocessor(); }
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnitBase, Preprocessor_, PreprocessorAction);

    virtual void Launch() = 0;
    void Launch_() { Launch(); }
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnitBase, Launch_, LaunchAction);

    virtual hpx::future<void> Step() = 0;
    hpx::future<void> Step_() { return Step(); }
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnitBase, Step_, StepAction);

    virtual double ResidualL2() = 0;
    double ResidualL2_() { return ResidualL2(); }
    HPX_DEFINE_COMPONENT_ACTION(HPXSimulationUnitBase, ResidualL2_, ResidualL2Action);
};

class HPXSimulationUnitClient
    : public hpx::components::client_base<HPXSimulationUnitClient, HPXSimulationUnitBase> {
  private:
    using BaseType = hpx::components::client_base<HPXSimulationUnitClient, HPXSimulationUnitBase>;

  public:
    HPXSimulationUnitClient() = default;
    HPXSimulationUnitClient(hpx::future<hpx::id_type>&& id) : BaseType(std::move(id)) {}
    HPXSimulationUnitClient(hpx::id_type&& id) : BaseType(std::move(id)) {}

    static constexpr const char* GetBasename() { return "Simulation_Unit_Client_"; }

    hpx::future<void> Preprocessor() {
        using ActionType = typename HPXSimulationUnitBase::PreprocessorAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Launch() {
        using ActionType = typename HPXSimulationUnitBase::LaunchAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<void> Step() {
        using ActionType = typename HPXSimulationUnitBase::StepAction;
        return hpx::async<ActionType>(this->get_id());
    }

    hpx::future<double> ResidualL2() {
        using ActionType = typename HPXSimulationUnitBase::ResidualL2Action;
        return hpx::async<ActionType>(this->get_id());
    }
};

template <typename ProblemType>
class HPXSimulationUnit;

struct HPXSimulationUnitFactory {
    static HPXSimulationUnitClient Create(const hpx::naming::id_type& here,
                                          const std::string& input_string,
                                          const uint locality_id,
                                          const uint submesh_id);
private:
    template <typename ProblemType>
    static HPXSimulationUnitClient CreateSimulationUnit(const hpx::naming::id_type& here,
                                                        const std::string& input_string,
                                                        const uint locality_id,
                                                        const uint submesh_id) {
        return HPXSimulationUnitFactory::CreateSimulationUnitImpl<ProblemType>(here, input_string, locality_id, submesh_id, Utilities::is_defined<ProblemType>{});
    }

    template<typename ProblemType>
    static HPXSimulationUnitClient CreateSimulationUnitImpl(const hpx::naming::id_type& here,
                                                            const std::string& input_string,
                                                            const uint locality_id,
                                                            const uint submesh_id,
                                                            std::true_type) {
        return hpx::components::new_<HPXSimulationUnit<ProblemType>>(here, input_string, locality_id, submesh_id);
    }

    template<typename ProblemType>
    static HPXSimulationUnitClient CreateSimulationUnitImpl(const hpx::naming::id_type& here,
                                                            const std::string& input_string,
                                                            const uint locality_id,
                                                            const uint submesh_id,
                                                            std::false_type) {
        throw std::runtime_error( "Problem class not supported, please check cmake configuration for proper support\n");
        return HPXSimulationUnitClient();
}
};


#endif