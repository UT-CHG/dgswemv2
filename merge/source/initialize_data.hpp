#ifndef INITIALIZE_DATA_HPP
#define INITIALIZE_DATA_HPP

#include "geometry/reader/adcirc_format.hpp"
#include "geometry/shape_function/linear_interpolant.hpp"
#include "geometry/element/l2_projection.hpp"

#include "problem/swe/boundary_types.hpp"
#include "problem/swe/swe_boundary_conditions.hpp"

template<typename MeshType>
void initialize_data(MeshType& mesh, Geometry::AdcircFormat& mesh_file, Stepper& stepper)
{
  using Point = std::array<double,2>;

  //begin by reading in the bathymetry from the fort.14
  using BathymetryInfo = std::pair< std::array<Point,3>, std::array<double,3> >;
  std::unordered_map<int, BathymetryInfo > bathymetry_vrtxs;

  for ( auto elt_info : mesh_file.elements ) {
    bathymetry_vrtxs.insert({elt_info.first, BathymetryInfo() });
    for ( uint node = 0; node < 3; ++node ) {
      std::array<double,3>& node_detail = mesh_file.nodes.at(elt_info.second[node + 1]);

      bathymetry_vrtxs[elt_info.first].first[node] = { node_detail[0], node_detail[1] };
      bathymetry_vrtxs[elt_info.first].second[node]= node_detail[2];
    }
  }


  //Use RAII to initialize all the data as well as set the appropriate boundary conditions
  mesh.call_for_each_element([&stepper, &bathymetry_vrtxs](auto& elt) {
      std::size_t p = elt.get_p();
      std::size_t ndof = (p + 1)*(p+2)/2;

      elt.data = std::make_unique<std::decay_t<decltype(*(elt.data))> >(ndof, elt.get_num_area_gauss_points(),
                                                      stepper.get_num_stages() );


      int id = elt.get_id();
      if ( !bathymetry_vrtxs.count(id) ) {
        throw std::logic_error("Error: could not find bathymetry for element with id: " + id);
      }

      std::function<double(Point)> bath_func = Geometry::Shape::linear_interpolant( bathymetry_vrtxs[id].first,
                                                                                    bathymetry_vrtxs[id].second);

      std::vector<double> bath_modes = Geometry::Element::L2_projection( elt, bath_func );

      elt.fill_area_gauss_points( bath_modes, elt.data->bath_at_gp );

      {
        elt.data->bath_deriv_wrt_x_at_gp.reserve(elt.get_num_area_gauss_points() );
        elt.data->bath_deriv_wrt_y_at_gp.reserve(elt.get_num_area_gauss_points() );

        std::vector<std::array<double,2> > Dbath = elt.get_gradient_at_gp(bath_modes);
        for ( auto& grad_eval : Dbath ) {
          elt.data->bath_deriv_wrt_x_at_gp.push_back(grad_eval[0]);
          elt.data->bath_deriv_wrt_y_at_gp.push_back(grad_eval[1]);
        }
      }


      //this is where we would set other initial conditions
      std::fill( elt.data->state[0].ze.begin(), elt.data->state[0].ze.end(), 0);
      std::fill( elt.data->state[0].qx.begin(), elt.data->state[0].qx.end(), 0);
      std::fill( elt.data->state[0].qy.begin(), elt.data->state[0].qy.end(), 0);
    });


  //for the edges we can use the bathymetry at the area Gauss points to compute the gauss points at the edge
  mesh.call_for_each_interior_edge([](auto& edg) {

      edg.data = std::make_unique<std::decay_t<decltype(*(edg.data))> >(edg.get_num_surface_gauss_points() );

      std::vector<double>& bath_at_area_gp = edg._elts.first->data->bath_at_gp;

      std::vector<double> bath_modes = Geometry::Element::L2_projection( *(edg._elts.first), bath_at_area_gp );

      edg.fill_surface_gauss_points( 0, bath_modes, edg.data->bath_at_gp );

    });

  mesh.call_for_each_boundary_edge([&mesh_file](auto& edg) {
      int elt_id = edg._elt->get_id();
      int local_fid = edg.get_local_fid();

      std::array<int,2> node_pair { mesh_file.elements.at(elt_id).at(local_fid+1),
          mesh_file.elements.at(elt_id).at((local_fid+1)%3 + 1) };

      SWE::BOUNDARY_TYPE ibtype = mesh_file.get_ibtype(node_pair);
      switch (ibtype)
        {
        case SWE::BOUNDARY_TYPE::_land_boundary:
          edg.data = std::make_unique<SWE::LandBoundary>(edg.get_num_surface_gauss_points() );
          break;
        case SWE::BOUNDARY_TYPE::_tidal_boundary:
          std::vector<double> freq_comp { 0.0001405257046694307, 0.0002810514093388614,
              0.0004215771140082921, 0.0005621028186777228, 0.0007026285233471535,
              0.0008431542280165842, 0.0009836799326860149, 0.00112420563736,
              0.00126473134202};

          edg.data = std::make_unique<SWE::TidalBoundary>(edg.get_num_surface_gauss_points(),
                                                          0.2763, freq_comp, 0., 21600.);
          break;
        }
      //compute the gauss points at the edge
      std::vector<double>& bath_at_area_gp = edg._elt->data->bath_at_gp;

      std::vector<double> bath_modes = Geometry::Element::L2_projection( *(edg._elt), bath_at_area_gp );

      edg.fill_surface_gauss_points( 0, bath_modes, edg.data->bath_at_gp );
    });
}
#endif
