#ifndef EHDG_GN_DBC_DISTRIBUTED_HPP
#define EHDG_GN_DBC_DISTRIBUTED_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "communication/db_data_exchanger.hpp"

namespace GN {
namespace EHDG {
namespace DBC {
class Distributed : public SWE::EHDG::DBC::Distributed {
  public:
    Distributed() = default;
    Distributed(const DBDataExchanger& exchanger) : SWE::EHDG::DBC::Distributed(exchanger) {}
};
}
}
}

#endif