#include "general_definitions.hpp"

template<typename ElementType>
void scrutinize_solution(const Stepper&, ElementType&);

template<typename MeshType>
void run_simulation(double time_end, Stepper& rk_stepper, MeshType& mesh)
{
  //we write these gross looking wrapper functions to append the rk_stepper in a way that allows us to keep the
  //the nice std::for_each notation without having to define rk_stepper within each element

  auto volume_kernel = [&rk_stepper](auto& elt) {
    SWE::volume_kernel(rk_stepper,elt);
  };

  auto source_kernel = [&rk_stepper](auto& elt) {
    SWE::source_kernel(rk_stepper,elt);
  };

  auto interface_kernel = [&rk_stepper](auto& intface) {
    SWE::interface_kernel(rk_stepper,intface);
  };

  auto boundary_kernel = [&rk_stepper](auto& bound) {
    SWE::boundary_kernel(rk_stepper,bound);
  };

  auto update_kernel = [&rk_stepper](auto& elt) {
    SWE::update_kernel(rk_stepper, elt);
  };

  auto swap_states = [&rk_stepper](auto& elt) {
    SWE::swap_states(rk_stepper, elt);
  };

  /*auto scru_sol = [&rk_stepper](auto& elt) {
    scrutinize_solution(rk_stepper, elt);
  };*/

  uint nsteps = std::ceil(time_end/rk_stepper.get_dt());
  uint n_stages = rk_stepper.get_num_stages();

  auto resize_data_container = [n_stages](auto& elt) { 
    elt.data.resize(n_stages + 1);
  };
 
  mesh.CallForEachElement(resize_data_container);
  
  SWE::write_VTK_data(rk_stepper, mesh);

  for ( uint step = 1; step <= nsteps; ++step ) {
    for ( uint stage = 0; stage < rk_stepper.get_num_stages(); ++stage ) {

      mesh.CallForEachElement(volume_kernel);

      mesh.CallForEachInterface(interface_kernel);

      mesh.CallForEachBoundary(boundary_kernel);
      
      mesh.CallForEachElement(update_kernel);

      ++rk_stepper;
    }

    mesh.CallForEachElement(swap_states);

    if ( step % 300 == 0 ) {
      std::cout << "Step: " << step << "\n";
      SWE::write_VTK_data(rk_stepper, mesh);
  }
}
}

template<typename ElementType>
void scrutinize_solution(const Stepper& rk_stepper, ElementType& elt)
{
  uint rk_stage = rk_stepper.get_stage();

  auto& state = elt.data->state[rk_stage];

  for ( auto& ze_mode : state.ze ) {
    if ( isnan(ze_mode) ) {
      std::cerr << "Error: found isnan ze at Element " << elt.get_id();
      std::cerr << "       At rk_stage: " << rk_stage << "\n";
    }
  }

  for ( auto& qx_mode : state.qx ) {
    if ( isnan(qx_mode) ) {
      std::cerr << "Error: found isnan qx at Element " << elt.get_id();
      std::cerr << "       At rk_stage: " << rk_stage << "\n";
    }
  }

  for ( auto& qy_mode : state.qy ) {
    if ( isnan(qy_mode) ) {
      std::cerr << "Error: found isnan qy at Element " << elt.get_id();
      std::cerr << "       At rk_stage: " << rk_stage << "\n";
    }
  }

  for ( auto& rhs_ze_mode : state.rhs_ze ) {
    if ( isnan(rhs_ze_mode) ) {
      std::cerr << "Error: found isnan rhs_ze at Element " << elt.get_id();
      std::cerr << "       At rk_stage: " << rk_stage << "\n";
    }
  }

  for ( auto& rhs_qx_mode : state.rhs_qx ) {
    if ( isnan(rhs_qx_mode) ) {
      std::cerr << "Error: found isnan rhs_qx at Element " << elt.get_id();
      std::cerr << "       At rk_stage: " << rk_stage << "\n";
    }
  }

  for ( auto& rhs_qy_mode : state.rhs_qy ) {
    if ( isnan(rhs_qy_mode) ) {
      std::cerr << "Error: found isnan rhs_qy at Element " << elt.get_id();
      std::cerr << "       At rk_stage: " << rk_stage << "\n";
    }
  }
}