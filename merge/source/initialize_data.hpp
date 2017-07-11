#ifndef INITIALIZE_DATA_HPP
#define INITIALIZE_DATA_HPP

template <typename MeshType>
void initialize_data(MeshType& mesh, AdcircFormat& mesh_file) {
    mesh.CallForEachElement(
      [](auto& elt) {
        elt.data.initialize();
      }
    );

    std::unordered_map<uint, std::vector<double> > bathymetry;

    for (auto elt : mesh_file.elements) {
        bathymetry.insert({ elt.first, std::vector<double>(3) });
        for (uint node = 0; node < 3; node++) {
            bathymetry[elt.first][node] = mesh_file.nodes.at(elt.second[node + 1])[2];
        }
    }

    mesh.CallForEachElement([&bathymetry](auto& elt) {
        uint id = elt.GetID();

        if (!bathymetry.count(id)) {
            throw std::logic_error("Error: could not find bathymetry for element with id: " + id);
        }

        elt.data.state[0].bath = elt.L2Projection(bathymetry[id]);

        elt.ComputeUgp(elt.data.state[0].bath, elt.data.internal.bath_at_gp);

        elt.ComputeDUgp(GlobalCoord::x, elt.data.state[0].bath, elt.data.internal.bath_deriv_wrt_x_at_gp);
        elt.ComputeDUgp(GlobalCoord::y, elt.data.state[0].bath, elt.data.internal.bath_deriv_wrt_y_at_gp);

        //this is where we would set other initial conditions
        std::fill(elt.data.state[0].ze.begin(), elt.data.state[0].ze.end(), 0);
        std::fill(elt.data.state[0].qx.begin(), elt.data.state[0].qx.end(), 0);
        std::fill(elt.data.state[0].qy.begin(), elt.data.state[0].qy.end(), 0);
    });

    mesh.CallForEachInterface([](auto& intface) {
        intface.ComputeUgpIN(intface.data_in.state[0].bath, intface.data_in.boundary.bath_at_gp);
        intface.ComputeUgpEX(intface.data_ex.state[0].bath, intface.data_ex.boundary.bath_at_gp);
    });

    mesh.CallForEachBoundary([](auto& bound) {
        bound.ComputeUgp(bound.data.state[0].bath, bound.data.boundary.bath_at_gp);
        /*
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
        }*/
    });
}

#endif