#ifndef SWE_KERNELS_HPP
#define SWE_KERNELS_HPP

#include "swe_LLF_flux.hpp"

namespace SWE {
	template<typename RawBoundaryType>
	void Problem::create_boundaries_kernel(mesh_type& mesh, std::map<unsigned char, std::vector<RawBoundaryType>>& pre_boundaries) {
		uint n_bound_old_land = 0;
		uint n_bound_old_tidal = 0;

		for (auto it = pre_boundaries.begin(); it != pre_boundaries.end(); it++) {
			switch (it->first) {
			case SWE::land:
				using BoundaryTypeLand = Geometry::Boundary<1, Integration::GaussLegendre_1D, SWE::Data, SWE::Land>;

				n_bound_old_land = mesh.GetNumberBoundaries();

				for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
					mesh.template CreateBoundary<BoundaryTypeLand>(*itt);
				}

				std::cout << "Number of land boundaries: " << mesh.GetNumberBoundaries() - n_bound_old_land << "\n";

				break;
			case SWE::tidal:
				using BoundaryTypeTidal = Geometry::Boundary<1, Integration::GaussLegendre_1D, SWE::Data, SWE::Tidal>;

				n_bound_old_tidal = mesh.GetNumberBoundaries();

				for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
					mesh.template CreateBoundary<BoundaryTypeTidal>(*itt);
				}

				std::cout << "Number of tidal boundaries: " << mesh.GetNumberBoundaries() - n_bound_old_tidal << "\n";

				break;
			}
		}
	}

	template<typename ElementType>
	void Problem::volume_kernel(const Stepper& stepper, ElementType& elt) {
		const uint stage = stepper.get_stage();

		auto& state = elt.data.state[stage];
		auto& internal = elt.data.internal;

		//get state at Gauss points
		elt.ComputeUgp(state.ze, internal.ze_at_gp);
		elt.ComputeUgp(state.qx, internal.qx_at_gp);
		elt.ComputeUgp(state.qy, internal.qy_at_gp);

		//assemble flux
		for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
			internal.water_column_hgt_at_gp[gp] = internal.ze_at_gp[gp] + internal.bath_at_gp[gp];

			internal.ze_flux_at_gp[GlobalCoord::x][gp] = internal.qx_at_gp[gp];
			internal.ze_flux_at_gp[GlobalCoord::y][gp] = internal.qy_at_gp[gp];

			internal.qx_flux_at_gp[GlobalCoord::x][gp] = std::pow(internal.qx_at_gp[gp], 2) / internal.water_column_hgt_at_gp[gp] +
				Global::g * (0.5 * std::pow(internal.ze_at_gp[gp], 2) + internal.ze_at_gp[gp] * internal.bath_at_gp[gp]);
			internal.qx_flux_at_gp[GlobalCoord::y][gp] = internal.qx_at_gp[gp] * internal.qy_at_gp[gp] / internal.water_column_hgt_at_gp[gp];

			internal.qy_flux_at_gp[GlobalCoord::x][gp] = internal.qx_at_gp[gp] * internal.qy_at_gp[gp] / internal.water_column_hgt_at_gp[gp];
			internal.qy_flux_at_gp[GlobalCoord::y][gp] = std::pow(internal.qy_at_gp[gp], 2) / internal.water_column_hgt_at_gp[gp] +
				Global::g * (0.5 * std::pow(internal.ze_at_gp[gp], 2) + internal.ze_at_gp[gp] * internal.bath_at_gp[gp]);
		}

		//skip dof = 0, which is a constant and thus trivially 0 NOT ALWAYS!
		for (uint dof = 1; dof < elt.data.get_ndof(); ++dof) {
			state.rhs_ze[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.ze_flux_at_gp[GlobalCoord::x]) +
				elt.IntegrationDPhi(GlobalCoord::y, dof, internal.ze_flux_at_gp[GlobalCoord::y]);

			state.rhs_qx[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.qx_flux_at_gp[GlobalCoord::x]) +
				elt.IntegrationDPhi(GlobalCoord::y, dof, internal.qx_flux_at_gp[GlobalCoord::y]);

			state.rhs_qy[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.qy_flux_at_gp[GlobalCoord::x]) +
				elt.IntegrationDPhi(GlobalCoord::y, dof, internal.qy_flux_at_gp[GlobalCoord::y]);
		}
	}

	template<typename ElementType>
	void Problem::source_kernel(const Stepper& stepper, ElementType& elt) {
		const uint stage = stepper.get_stage();

		auto& state = elt.data.state[stage];
		auto& internal = elt.data.internal;

		//note we assume that the values at gauss points have already been computed
		for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
			//compute contribution of hydrostatic pressure
			internal.qx_source_term_at_gp[gp] = Global::g * internal.bath_deriv_wrt_x_at_gp[gp] * internal.ze_at_gp[gp];
			internal.qy_source_term_at_gp[gp] = Global::g * internal.bath_deriv_wrt_y_at_gp[gp] * internal.ze_at_gp[gp];

			double u_at_gp = internal.qx_at_gp[gp] / internal.water_column_hgt_at_gp[gp];
			double v_at_gp = internal.qy_at_gp[gp] / internal.water_column_hgt_at_gp[gp];

			//compute bottom friction contribution
			double bottom_friction_stress = Global::Cf *
				std::hypot(u_at_gp, v_at_gp) / internal.water_column_hgt_at_gp[gp];

			internal.qx_source_term_at_gp[gp] -= bottom_friction_stress * internal.qx_at_gp[gp];
			internal.qy_source_term_at_gp[gp] -= bottom_friction_stress * internal.qy_at_gp[gp];
		}

		for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
			//current we don't have any source terms that affect ze
			state.rhs_qx[dof] += elt.IntegrationPhi(dof, internal.qx_source_term_at_gp);
			state.rhs_qy[dof] += elt.IntegrationPhi(dof, internal.qy_source_term_at_gp);
		}
	}

	template<typename InterfaceType>
	void Problem::interface_kernel(const Stepper& stepper, InterfaceType& intface) {
		const uint stage = stepper.get_stage();

		auto& state_in = intface.data_in.state[stage];
		auto& boundary_in = intface.data_in.boundary;

		auto& state_ex = intface.data_ex.state[stage];
		auto& boundary_ex = intface.data_ex.boundary;

		intface.ComputeUgpIN(state_in.ze, boundary_in.ze_at_gp);
		intface.ComputeUgpIN(state_in.qx, boundary_in.qx_at_gp);
		intface.ComputeUgpIN(state_in.qy, boundary_in.qy_at_gp);

		intface.ComputeUgpEX(state_ex.ze, boundary_ex.ze_at_gp);
		intface.ComputeUgpEX(state_ex.qx, boundary_ex.qx_at_gp);
		intface.ComputeUgpEX(state_ex.qy, boundary_ex.qy_at_gp);

		//assemble numerical fluxes
		for (uint gp = 0; gp < intface.data_in.get_ngp_boundary(); ++gp) {
			LLF_flux(boundary_in.ze_at_gp[gp], boundary_ex.ze_at_gp[gp],
				boundary_in.qx_at_gp[gp], boundary_ex.qx_at_gp[gp],
				boundary_in.qy_at_gp[gp], boundary_ex.qy_at_gp[gp],
				boundary_in.bath_at_gp[gp], intface.surface_normal[gp],
				boundary_in.ze_numerical_flux_at_gp[gp],
				boundary_in.qx_numerical_flux_at_gp[gp],
				boundary_in.qy_numerical_flux_at_gp[gp]);
		}

		//now compute contributions to the righthand side
		for (uint dof = 0; dof < intface.data_in.get_ndof(); ++dof) {
			state_in.rhs_ze[dof] -= intface.IntegrationPhiIN(dof, boundary_in.ze_numerical_flux_at_gp);
			state_in.rhs_qx[dof] -= intface.IntegrationPhiIN(dof, boundary_in.qx_numerical_flux_at_gp);
			state_in.rhs_qy[dof] -= intface.IntegrationPhiIN(dof, boundary_in.qy_numerical_flux_at_gp);
		}

		for (uint dof = 0; dof < intface.data_ex.get_ndof(); ++dof) {
			state_ex.rhs_ze[dof] += intface.IntegrationPhiEX(dof, boundary_in.ze_numerical_flux_at_gp);
			state_ex.rhs_qx[dof] += intface.IntegrationPhiEX(dof, boundary_in.qx_numerical_flux_at_gp);
			state_ex.rhs_qy[dof] += intface.IntegrationPhiEX(dof, boundary_in.qy_numerical_flux_at_gp);
		}
	}

	template<typename BoundaryType>
	void Problem::boundary_kernel(const Stepper& stepper, BoundaryType& bound) {
		const uint stage = stepper.get_stage();

		auto& state = bound.data.state[stage];
		auto& boundary = bound.data.boundary;

		bound.ComputeUgp(state.ze, boundary.ze_at_gp);
		bound.ComputeUgp(state.qx, boundary.qx_at_gp);
		bound.ComputeUgp(state.qy, boundary.qy_at_gp);

		double ze_ex, qx_ex, qy_ex;
		for (uint gp = 0; gp < bound.data.get_ngp_boundary(); ++gp) {
			bound.boundary_condition.set_ex(stepper, bound.surface_normal[gp],
				boundary.ze_at_gp[gp], boundary.qx_at_gp[gp], boundary.qy_at_gp[gp],
				ze_ex, qx_ex, qy_ex);

			LLF_flux(boundary.ze_at_gp[gp], ze_ex,
				boundary.qx_at_gp[gp], qx_ex,
				boundary.qy_at_gp[gp], qy_ex,
				boundary.bath_at_gp[gp], bound.surface_normal[gp],
				boundary.ze_numerical_flux_at_gp[gp],
				boundary.qx_numerical_flux_at_gp[gp],
				boundary.qy_numerical_flux_at_gp[gp]);
		}

		//now compute contributions to the righthand side
		for (uint dof = 0; dof < bound.data.get_ndof(); ++dof) {
			state.rhs_ze[dof] -= bound.IntegrationPhi(dof, boundary.ze_numerical_flux_at_gp);
			state.rhs_qx[dof] -= bound.IntegrationPhi(dof, boundary.qx_numerical_flux_at_gp);
			state.rhs_qy[dof] -= bound.IntegrationPhi(dof, boundary.qy_numerical_flux_at_gp);
		}
	}

	template<typename ElementType>
	void Problem::update_kernel(const Stepper& stepper, ElementType& elt) {
		const uint stage = stepper.get_stage();
		double dt = stepper.get_dt();

		auto& state = elt.data.state;
		auto& curr_state = elt.data.state[stage];
		auto& next_state = elt.data.state[stage + 1];

		curr_state.rhs_ze = elt.SolveLSE(curr_state.rhs_ze);
		curr_state.rhs_qx = elt.SolveLSE(curr_state.rhs_qx);
		curr_state.rhs_qy = elt.SolveLSE(curr_state.rhs_qy);

		std::fill(next_state.ze.begin(), next_state.ze.end(), 0);
		std::fill(next_state.qx.begin(), next_state.qx.end(), 0);
		std::fill(next_state.qy.begin(), next_state.qy.end(), 0);

		for (uint s = 0; s <= stage; ++s) {
			for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
				next_state.ze[dof] += stepper.ark[stage][s] * state[s].ze[dof]
					+ dt*stepper.brk[stage][s] * state[s].rhs_ze[dof];

				next_state.qx[dof] += stepper.ark[stage][s] * state[s].qx[dof]
					+ dt*stepper.brk[stage][s] * state[s].rhs_qx[dof];

				next_state.qy[dof] += stepper.ark[stage][s] * state[s].qy[dof]
					+ dt*stepper.brk[stage][s] * state[s].rhs_qy[dof];
			}
		}
	}

	template<typename ElementType>
	void Problem::swap_states_kernel(const Stepper& stepper, ElementType& elt) {
		uint n_stages = stepper.get_num_stages();
		auto& state = elt.data.state;

		std::swap(state[0].ze, state[n_stages].ze);
		std::swap(state[0].qx, state[n_stages].qx);
		std::swap(state[0].qy, state[n_stages].qy);
	}

	template<typename ElementType>
	void Problem::scrutinize_solution_kernel(const Stepper& stepper, ElementType& elt) {
		uint stage = stepper.get_stage();

		auto& state = elt.data.state[stage];

		for (auto& ze_mode : state.ze) {
			if (isnan(ze_mode)) {
				std::cerr << "Error: found isnan ze at Element " << elt.get_id();
				std::cerr << "       At stage: " << stage << "\n";
			}
		}

		for (auto& qx_mode : state.qx) {
			if (isnan(qx_mode)) {
				std::cerr << "Error: found isnan qx at Element " << elt.get_id();
				std::cerr << "       At stage: " << stage << "\n";
			}
		}

		for (auto& qy_mode : state.qy) {
			if (isnan(qy_mode)) {
				std::cerr << "Error: found isnan qy at Element " << elt.get_id();
				std::cerr << "       At stage: " << stage << "\n";
			}
		}

		for (auto& rhs_ze_mode : state.rhs_ze) {
			if (isnan(rhs_ze_mode)) {
				std::cerr << "Error: found isnan rhs_ze at Element " << elt.get_id();
				std::cerr << "       At stage: " << stage << "\n";
			}
		}

		for (auto& rhs_qx_mode : state.rhs_qx) {
			if (isnan(rhs_qx_mode)) {
				std::cerr << "Error: found isnan rhs_qx at Element " << elt.get_id();
				std::cerr << "       At stage: " << stage << "\n";
			}
		}

		for (auto& rhs_qy_mode : state.rhs_qy) {
			if (isnan(rhs_qy_mode)) {
				std::cerr << "Error: found isnan rhs_qy at Element " << elt.get_id();
				std::cerr << "       At stage: " << stage << "\n";
			}
		}
	}

	template<typename ElementType>
	void Problem::extract_VTK_data_kernel(ElementType& elt, Array2D<double>& cell_data, Array2D<double>& point_data) {
		elt.WriteCellDataVTK(elt.data.state[0].ze, cell_data[0]);
		elt.WriteCellDataVTK(elt.data.state[0].qx, cell_data[1]);
		elt.WriteCellDataVTK(elt.data.state[0].qy, cell_data[2]);
		elt.WriteCellDataVTK(elt.data.state[0].bath, cell_data[3]);

		elt.WritePointDataVTK(elt.data.state[0].ze, point_data[0]);
		elt.WritePointDataVTK(elt.data.state[0].qx, point_data[1]);
		elt.WritePointDataVTK(elt.data.state[0].qy, point_data[2]);
		elt.WritePointDataVTK(elt.data.state[0].bath, point_data[3]);
	}

	template<typename MeshType>
	void Problem::write_VTK_data_kernel(const Stepper& stepper, MeshType& mesh) {
		Array2D<double> cell_data;
		Array2D<double> point_data;

		cell_data.resize(4);
		point_data.resize(4);

		auto extract_VTK_data_kernel = [&cell_data, &point_data](auto& elt) {
			Problem::extract_VTK_data_kernel(elt, cell_data, point_data);
		};

		mesh.CallForEachElement(extract_VTK_data_kernel);

		std::string file_name = "output/data.vtk";
		std::ofstream file(file_name);

		file << "CELL_DATA " << (*cell_data.begin()).size() << '\n';
		file << "SCALARS ze_cell float 1\n";
		file << "LOOKUP_TABLE default\n";
		for (auto it = cell_data[0].begin(); it != cell_data[0].end(); it++) file << *it << '\n';

		file << "SCALARS qx_cell float 1\n";
		file << "LOOKUP_TABLE default\n";
		for (auto it = cell_data[1].begin(); it != cell_data[1].end(); it++) file << *it << '\n';

		file << "SCALARS qy_cell float 1\n";
		file << "LOOKUP_TABLE default\n";
		for (auto it = cell_data[2].begin(); it != cell_data[2].end(); it++) file << *it << '\n';

		file << "SCALARS bath_cell float 1\n";
		file << "LOOKUP_TABLE default\n";
		for (auto it = cell_data[3].begin(); it != cell_data[3].end(); it++) file << *it << '\n';

		file << "POINT_DATA " << (*point_data.begin()).size() << '\n';
		file << "SCALARS ze_point float 1\n";
		file << "LOOKUP_TABLE default\n";
		for (auto it = point_data[0].begin(); it != point_data[0].end(); it++) file << *it << '\n';

		file << "SCALARS qx_point float 1\n";
		file << "LOOKUP_TABLE default\n";
		for (auto it = point_data[1].begin(); it != point_data[1].end(); it++) file << *it << '\n';

		file << "SCALARS qy_point float 1\n";
		file << "LOOKUP_TABLE default\n";
		for (auto it = point_data[2].begin(); it != point_data[2].end(); it++) file << *it << '\n';

		file << "SCALARS bath_point float 1\n";
		file << "LOOKUP_TABLE default\n";
		for (auto it = point_data[3].begin(); it != point_data[3].end(); it++) file << *it << '\n';

		file.close();

		std::string file_name_geom = "output/geometry.vtk";
		std::string file_name_data = "output/data.vtk";

		std::ifstream file_geom(file_name_geom, std::ios_base::binary);
		std::ifstream file_data(file_name_data, std::ios_base::binary);

		uint n_step = (uint)(stepper.get_t_at_curr_stage() / stepper.get_dt());

		std::string file_name_merge = "output/mesh_data_" + std::to_string(n_step) + ".vtk";
		std::ofstream file_merge(file_name_merge, std::ios_base::binary);

		file_merge << file_geom.rdbuf() << file_data.rdbuf();
		file_merge.close();
	}

	template<typename ElementType>
	void Problem::extract_modal_data_kernel(ElementType& elt, std::vector<std::pair<uint, Array2D<double>>>& modal_data) {
		modal_data.push_back(
			std::make_pair(elt.GetID(), 
				Array2D<double>{elt.data.state[0].ze, elt.data.state[0].qx,
								elt.data.state[0].qy, elt.data.state[0].bath}
			)
		);
	}

	template<typename MeshType>
	void Problem::write_modal_data_kernel(const Stepper& stepper, MeshType& mesh) {
		std::vector<std::pair<uint, Array2D<double>>> modal_data;

		auto extract_modal_data_kernel = [&modal_data](auto& elt) {
			Problem::extract_modal_data_kernel(elt, modal_data);
		};

		mesh.CallForEachElement(extract_modal_data_kernel);

		std::string file_name = "output/modal_ze.txt";
		std::ofstream file(file_name, std::ios::app);

		file << std::to_string(stepper.get_t_at_curr_stage()) << '\n';
		for (auto it = modal_data.begin(); it != modal_data.end(); it++) {
			for(auto itt = (*it).second[0].begin(); itt != (*it).second[0].end(); itt++) {
				file << (*it).first << ' '<< std::scientific << (*itt) << '\n';
			}
		}

		file.close();

		file_name = "output/modal_qx.txt";
		file = std::ofstream(file_name, std::ios::app);

		file << std::to_string(stepper.get_t_at_curr_stage()) << '\n';
		for (auto it = modal_data.begin(); it != modal_data.end(); it++) {
			for(auto itt = (*it).second[1].begin(); itt != (*it).second[1].end(); itt++) {
				file << (*it).first << ' '<< std::scientific << (*itt) << '\n';
			}
		}

		file.close();

		file_name = "output/modal_qy.txt";
		file = std::ofstream(file_name, std::ios::app);

		file << std::to_string(stepper.get_t_at_curr_stage()) << '\n';
		for (auto it = modal_data.begin(); it != modal_data.end(); it++) {
			for(auto itt = (*it).second[2].begin(); itt != (*it).second[2].end(); itt++) {
				file << (*it).first << ' '<< std::scientific << (*itt) << '\n';
			}
		}

		file.close();

		file_name = "output/modal_bath.txt";
		file = std::ofstream(file_name, std::ios::app);

		file << std::to_string(stepper.get_t_at_curr_stage()) << '\n';
		for (auto it = modal_data.begin(); it != modal_data.end(); it++) {
			for(auto itt = (*it).second[3].begin(); itt != (*it).second[3].end(); itt++) {
				file << (*it).first << ' '<< std::scientific << (*itt) << '\n';
			}
		}

		file.close();
	}
}

#endif