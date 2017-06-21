#include "stepper.hpp"
#include "problem/SWE/swe_kernels.hpp"

#include "general_definitions.h"

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


template<typename MeshType>
void run_simulation(double time_end, Stepper& rk_stepper, MeshType& mesh)
{
  //we write these gross looking wrapper functions to append the rk_stepper in a way that allows us to keep the
  //the nice std::for_each notation without having to define rk_stepper within each element

  auto volume_kernel = [&rk_stepper](auto& elt) {
    SWE::volume_kernel(rk_stepper,elt);
  };

  auto edge_kernel = [&rk_stepper](auto& edg) {
    SWE::edge_kernel(rk_stepper,edg);
  };

  auto boundary_kernel = [&rk_stepper](auto& edg) {
    SWE::boundary_kernel(rk_stepper,edg);
  };

  auto update_kernel = [&rk_stepper](auto& elt) {
    SWE::update_kernel(rk_stepper,elt);
  };

  auto swap_states = [&rk_stepper](auto& elt) {
    SWE::swap_states(rk_stepper, elt);
  };

  auto scru_sol = [&rk_stepper](auto& elt) {
    scrutinize_solution(rk_stepper, elt);
  };

  uint nsteps = std::ceil(time_end/rk_stepper.get_dt());

  Writer::ParaviewWriter pv_writer(1);
  std::string fname = "DG.63";
  Writer::DG63Writer writer(fname);

  for ( uint step = 0; step < nsteps; ++step ) {
    for ( uint stage = 0; stage < rk_stepper.get_num_stages(); ++stage ) {

      mesh.call_for_each_element(volume_kernel);
      //      std::cerr << "After volume kernel\n";
      //mesh.call_for_each_element(scru_sol);

      mesh.call_for_each_interior_edge(edge_kernel);
      //std::cerr << "After edge kernel\n";
      //mesh.call_for_each_element(scru_sol);

      mesh.call_for_each_boundary_edge(boundary_kernel);
      //std::cerr << "After boundary kernel\n";
      //mesh.call_for_each_element(scru_sol);

      mesh.call_for_each_element(update_kernel);
      //std::cerr << "After update kernel\n";
      //mesh.call_for_each_element(scru_sol);

      ++rk_stepper;
    }

    mesh.call_for_each_element(swap_states);
    //mesh.call_for_each_element(scru_sol);

    if ( step % 360 == 0 ) {
      std::cout << "Step: " << step << "\n";
      std::string fname = "output/P" + std::to_string(step) + ".vtk"; std::string var_name = "Ze";

      writer.write_output(mesh, rk_stepper.get_t_at_curr_stage(), step);
      pv_writer.write_output(fname, mesh, var_name,
                          [](double ze, double qx, double qy, double bath)->double {
                            return ze;
                            });
    }
  }
}
