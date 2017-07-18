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

	mesh.CallForEachElement(
		[&bathymetry](auto& elt) {
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
		}
	);

	mesh.CallForEachInterface(
		[](auto& intface) {
			intface.ComputeUgpIN(intface.data_in.state[0].bath, intface.data_in.boundary.bath_at_gp);
			intface.ComputeUgpEX(intface.data_ex.state[0].bath, intface.data_ex.boundary.bath_at_gp);
		}
	);

	mesh.CallForEachBoundary(
		[](auto& bound) {
			bound.ComputeUgp(bound.data.state[0].bath, bound.data.boundary.bath_at_gp);
		}
	);
}

#endif