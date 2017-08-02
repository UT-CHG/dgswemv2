#ifndef INITIALIZE_DATA_HPP
#define INITIALIZE_DATA_HPP

#include "problem/SWE/swe_ic_src_functions.hpp"

template <typename MeshType>
void initialize_data(MeshType& mesh, const MeshMetaData& mesh_data) {
	mesh.CallForEachElement(
		[](auto& elt) {
			elt.data.initialize();
		}
	);

	std::unordered_map<uint, std::vector<double> > bathymetry;

	for (const auto& elt : mesh_data._elements) {
	        bathymetry.insert({ elt.first, mesh_data.GetBathymetry(elt.first) });
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

			auto ze_init = [](Point<2>& pt) {
				return SWE::ic_ze(0, pt);
			};

			elt.data.state[0].ze = elt.L2Projection(ze_init);

			auto qx_init = [](Point<2>& pt){
				return SWE::ic_qx(0, pt);
			};

			elt.data.state[0].qx = elt.L2Projection(qx_init);

			auto qy_init = [](Point<2>& pt){
				return SWE::ic_qy(0, pt);
			};

			elt.data.state[0].qy = elt.L2Projection(qy_init);
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
