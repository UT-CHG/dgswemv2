#ifndef SWE_DATA_HPP
#define SWE_DATA_HPP

#include "../../general_definitions.hpp"

namespace SWE {
	struct State {
		State(uint ndof)
		  : ze(ndof), qx(ndof), qy(ndof),
			rhs_ze(ndof), rhs_qx(ndof), rhs_qy(ndof)
		{}

		std::vector<double> ze;
		std::vector<double> qx;
		std::vector<double> qy;
		std::vector<double> rhs_ze;
		std::vector<double> rhs_qx;
		std::vector<double> rhs_qy;
	};

	struct Boundary
	{
		Boundary(uint n_gp)
		  : ze_at_gp(n_gp),
			qx_at_gp(n_gp),
			qy_at_gp(n_gp),
			bath_at_gp(n_gp),
			ze_numerical_flux_at_gp(n_gp),
			qx_numerical_flux_at_gp(n_gp),
			qy_numerical_flux_at_gp(n_gp)
		{}

		std::vector<double> ze_at_gp;
		std::vector<double> qx_at_gp;
		std::vector<double> qy_at_gp;

		std::vector<double> bath_at_gp;

		std::vector<double> ze_numerical_flux_at_gp;
		std::vector<double> qx_numerical_flux_at_gp;
		std::vector<double> qy_numerical_flux_at_gp;
	};

	struct Internal {
		Internal(uint n_gp)
		  : ze_flux_at_gp({ std::vector<double>(n_gp), std::vector<double>(n_gp) }),
			qx_flux_at_gp({ std::vector<double>(n_gp), std::vector<double>(n_gp) }),
			qy_flux_at_gp({ std::vector<double>(n_gp), std::vector<double>(n_gp) }),
			ze_source_term_at_gp(n_gp),
			qx_source_term_at_gp(n_gp),
			qy_source_term_at_gp(n_gp),
			ze_at_gp(n_gp),
			qx_at_gp(n_gp),
			qy_at_gp(n_gp),
			bath_at_gp(n_gp),
			bath_deriv_wrt_x_at_gp(n_gp),
			bath_deriv_wrt_y_at_gp(n_gp),
			water_column_hgt_at_gp(n_gp)
		{}

		std::array<std::vector<double>, 2> ze_flux_at_gp;
		std::array<std::vector<double>, 2> qx_flux_at_gp;
		std::array<std::vector<double>, 2> qy_flux_at_gp;

		std::vector<double> ze_source_term_at_gp;
		std::vector<double> qx_source_term_at_gp;
		std::vector<double> qy_source_term_at_gp;

		std::vector<double> ze_at_gp;
		std::vector<double> qx_at_gp;
		std::vector<double> qy_at_gp;

		std::vector<double> bath_at_gp;
		std::vector<double> bath_deriv_wrt_x_at_gp;
		std::vector<double> bath_deriv_wrt_y_at_gp;

		std::vector<double> water_column_hgt_at_gp;
	};

	struct Data {
		Data()
			: state(0, State(0)), internal(0), boundary(0)
		{}

		void resize(uint nstate) {
			this->state = std::vector<State>(nstate, State(ndof));
			this->internal = Internal(ngp_internal);
			this->boundary = Boundary(ngp_boundary);

			//SET INITIAL CONDITIONS FOR TESTING
			for (uint i = 0; i < this->ngp_internal; i++)
				this->internal.bath_at_gp[i] = 5.0;

			for (uint i = 0; i < this->ngp_boundary; i++)
				this->boundary.bath_at_gp[i] = 5.0;
		}

		std::vector<State> state;
		Internal internal;
		Boundary boundary;

		uint get_ngp_internal() { return ngp_internal; }
		uint get_ngp_boundary() { return ngp_boundary; }		
		uint get_ndof() { return ndof; }

		void set_ngp_internal(uint ngp) { this->ngp_internal = ngp; }
		void set_ngp_boundary(uint ngp) { this->ngp_boundary = ngp; }
		void set_ndof(uint ndof) { this->ndof = ndof; }

	private:
		uint ndof;
		uint ngp_boundary;
		uint ngp_internal;
	};
}

#endif