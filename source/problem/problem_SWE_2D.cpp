#include <cmath>
#include <algorithm> 

#include "../general_definitions.h"
#include "../class_problem.h"

PROBLEM::PROBLEM() {
	this->mesh = new MESH(2,0);
}

PROBLEM::~PROBLEM() {
	delete this->mesh;
}

void PROBLEM::Timestep() {
	std::vector <INTERFACE*> internal_interfaces(this->mesh->interfaces.find(INTERNAL_BOUNDARY)->second);

	for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
		it->second->ComputeInternalU(UA);
		it->second->ComputeBoundaryU(UA);

		it->second->ComputeInternalU(VA);
		it->second->ComputeBoundaryU(VA);

		it->second->ComputeInternalU(H);
		it->second->ComputeBoundaryU(H);
	}

	for (auto it = internal_interfaces.begin(); it != internal_interfaces.end(); it++) {
		this->InternalInterfaceDiffusion(*it);
	}

	int number_bf;
	double* RHS;

	for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
		number_bf = it->second->number_bf;
		RHS = it->second->RHS;

		//COMPUTING TXX
		for (int i = 0; i < number_bf; i++) {
			RHS[i] = (it->second->IntegrationInternalDPhiDX(UA, i) +
				it->second->IntegrationBoundaryNX(UA_AVG, i))*
				pow(it->second->u[SP][i], 2);
		}
		it->second->SolveLSE(TXX);

		it->second->ComputeInternalU(TXX);
		it->second->ComputeBoundaryU(TXX);

		//COMPUTING TXY
		for (int i = 0; i < number_bf; i++) {
			RHS[i] = it->second->IntegrationInternalDPhiDY(UA, i) +
				it->second->IntegrationBoundaryNY(UA_AVG, i);
		}
		it->second->SolveLSE(TXY);

		it->second->ComputeInternalU(TXY);
		it->second->ComputeBoundaryU(TXY);

		//COMPUTING TYX
		for (int i = 0; i < number_bf; i++) {
			RHS[i] = (it->second->IntegrationInternalDPhiDX(VA, i) +
				it->second->IntegrationBoundaryNX(VA_AVG, i))*
				pow(it->second->u[SP][i], 2);
		}
		it->second->SolveLSE(TYX);

		it->second->ComputeInternalU(TYX);
		it->second->ComputeBoundaryU(TYX);

		//COMPUTING TYY
		for (int i = 0; i < number_bf; i++) {
			RHS[i] = it->second->IntegrationInternalDPhiDY(VA, i) +
				it->second->IntegrationBoundaryNY(VA_AVG, i);
		}
		it->second->SolveLSE(TYY);

		it->second->ComputeInternalU(TYY);
		it->second->ComputeBoundaryU(TYY);
	}

	for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
		this->ComputeUVA(it->second);
		this->ComputeF(it->second);
	}

	for (auto it = internal_interfaces.begin(); it != internal_interfaces.end(); it++) {
		this->LLFNumericalFlux(*it);
	}

	for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
		number_bf = it->second->number_bf;
		RHS = it->second->RHS;
		
		//COMPUTING UA
		for (int i = 0; i < number_bf; i++) {
			RHS[i] = it->second->IntegrationInternalDPhiDX(F11, i) +
				it->second->IntegrationInternalDPhiDY(F12, i) -
				it->second->IntegrationInternalPhi(NUM_FLUX_UA, i);
		}
		it->second->SolveLSE(D_UA);

		//COMPUTING VA
		for (int i = 0; i < number_bf; i++) {
			RHS[i] = it->second->IntegrationInternalDPhiDX(F21, i) +
				it->second->IntegrationInternalDPhiDY(F22, i) -
				it->second->IntegrationInternalPhi(NUM_FLUX_VA, i);
		}
		it->second->SolveLSE(D_VA);

		//COMPUTING H
		for (int i = 0; i < number_bf; i++) {
			RHS[i] = it->second->IntegrationInternalDPhiDX(F31, i) +
				it->second->IntegrationInternalDPhiDY(F32, i) -
				it->second->IntegrationInternalPhi(NUM_FLUX_H, i);
		}
		it->second->SolveLSE(D_H);
	}

}

void PROBLEM::ComputeUVA(ELEMENT* element) {
	double** u = element->u_internal;
	
	for (int i = 0; i < element->number_gp_internal; i++) {
		u[U][i] = u[UA][i] / u[A][i];
		u[V][i] = u[VA][i] / u[A][i];
		u[A][i] = u[ZB][i] + u[H][i];
	}

	for (int j = 0; j < element->number_interfaces; j++) {
		u = element->u_boundary[j];

		for (int i = 0; i < element->number_gp_boundary; i++) {
			u[U][i] = u[UA][i] / u[A][i];
			u[V][i] = u[VA][i] / u[A][i];
			u[A][i] = u[ZB][i] + u[H][i];
		}
	}
}

void PROBLEM::ComputeF(ELEMENT* element) {
	double** u = element->u_internal;

	for (int i = 0; i < element->number_gp_internal; i++) {
		u[F11][i] = u[SP][i] * (u[UA][i] * u[U][i] + 0.5*GRAVITY*pow(u[A][i], 2)) - VISCOSITY*u[TXX][i];
		u[F21][i] = u[SP][i] * u[UA][i] * u[V][i] - VISCOSITY*u[TYX][i];
		u[F31][i] = u[SP][i] * u[UA][i];

		u[F12][i] = u[UA][i] * u[V][i] - VISCOSITY*u[TXY][i];
		u[F22][i] = u[VA][i] * u[V][i] + 0.5*GRAVITY*pow(u[A][i], 2) - VISCOSITY*u[TYY][i];
		u[F32][i] = u[VA][i];
	}

	for (int j = 0; j < element->number_interfaces; j++) {
		u = element->u_boundary[j];

		for (int i = 0; i < element->number_gp_boundary; i++) {
			u[F11][i] = u[SP][i] * (u[UA][i] * u[U][i] + 0.5*GRAVITY*pow(u[A][i], 2)) - VISCOSITY*u[TXX][i];
			u[F21][i] = u[SP][i] * u[UA][i] * u[V][i] - VISCOSITY*u[TYX][i];
			u[F31][i] = u[SP][i] * u[UA][i];

			u[F12][i] = u[UA][i] * u[V][i] - VISCOSITY*u[TXY][i];
			u[F22][i] = u[VA][i] * u[V][i] + 0.5*GRAVITY*pow(u[A][i], 2) - VISCOSITY*u[TYY][i];
			u[F32][i] = u[VA][i];
		}
	}
}

void PROBLEM::InternalInterfaceDiffusion(INTERFACE* interface) {
	double** u_in = interface->u_boundary_in;
	double** u_ex = interface->u_boundary_ex;

	int i_ex;

	for (int i = 0; i < interface->number_gp; i++) {
		i_ex = interface->number_gp - 1 - i;

		u_in[UA_AVG][i] = 0.5 * (u_in[UA][i] + u_ex[UA][i_ex]);
		u_ex[UA_AVG][i_ex] = u_in[UA_AVG][i];

		u_in[VA_AVG][i] = 0.5 * (u_in[VA][i] + u_ex[VA][i_ex]);
		u_ex[VA_AVG][i_ex] = u_in[VA_AVG][i];
	}
}

void PROBLEM::LLFNumericalFlux(INTERFACE* interface){
	double u_normal_in;
	double u_normal_ex;

	double sqrt_g_a_in;
	double sqrt_g_a_ex;

	double eigen_1;
	double eigen_2;
	double eigen_3;
	double eigen_4;
	double eigen_5;
	double eigen_6;

	double lambda;

	double** u_in = interface->u_boundary_in;
	double** u_ex = interface->u_boundary_ex;
	double* n_x = interface->normal_x;
	double* n_y = interface->normal_y;

	int i_ex;

	for (int i = 0; i < interface->number_gp; i++) {
		i_ex = interface->number_gp - 1 - i;

		u_normal_in = u_in[U][i] * n_x[i] + u_in[V][i] * n_y[i];
		sqrt_g_a_in = sqrt(GRAVITY*u_in[A][i] * (pow((n_x[i] * u_in[SP][i]), 2) + pow(n_y[i], 2)));

		u_normal_ex = u_ex[U][i_ex] * n_x[i_ex] + u_ex[V][i_ex] * n_y[i_ex];
		sqrt_g_a_ex = sqrt(GRAVITY*u_ex[A][i_ex] * (pow((n_x[i_ex] * u_ex[SP][i_ex]), 2) + pow(n_y[i_ex], 2)));

		eigen_1 = abs(u_normal_in + sqrt_g_a_in);
		eigen_2 = abs(u_normal_in);
		eigen_3 = abs(u_normal_in - sqrt_g_a_in);

		eigen_4 = abs(u_normal_ex + sqrt_g_a_ex);
		eigen_5 = abs(u_normal_ex);
		eigen_6 = abs(u_normal_ex - sqrt_g_a_ex);

		lambda = std::max(eigen_1, eigen_2);
		lambda = std::max(lambda, eigen_3);
		lambda = std::max(lambda, eigen_4);
		lambda = std::max(lambda, eigen_5);
		lambda = std::max(lambda, eigen_6);

		u_in[F11_AVG][i] = 0.5*(u_in[F11][i] + u_ex[F11][i_ex]);
		u_in[F12_AVG][i] = 0.5*(u_in[F12][i] + u_ex[F12][i_ex]);
		u_in[F21_AVG][i] = 0.5*(u_in[F21][i] + u_ex[F21][i_ex]);
		u_in[F22_AVG][i] = 0.5*(u_in[F22][i] + u_ex[F22][i_ex]);
		u_in[F31_AVG][i] = 0.5*(u_in[F31][i] + u_ex[F31][i_ex]);
		u_in[F32_AVG][i] = 0.5*(u_in[F32][i] + u_ex[F32][i_ex]);

		u_in[UA_JUMP][i] = u_in[UA][i] - u_ex[UA][i_ex];
		u_in[VA_JUMP][i] = u_in[VA][i] - u_ex[VA][i_ex];
		u_in[H_JUMP][i] = u_in[H][i] - u_ex[H][i_ex];

		u_in[NUM_FLUX_UA][i] = n_x[i]*u_in[F11_AVG][i] + n_y[i]*u_in[F12_AVG][i] + 
			0.5 * lambda * u_in[UA_JUMP][i];
		u_ex[NUM_FLUX_UA][i_ex] = -u_in[NUM_FLUX_UA][i];

		u_in[NUM_FLUX_VA][i] = n_x[i]*u_in[F21_AVG][i] + n_y[i]*u_in[F22_AVG][i] +
			0.5 * lambda * u_in[VA_JUMP][i];
		u_ex[NUM_FLUX_VA][i_ex] = -u_in[NUM_FLUX_VA][i];

		u_in[NUM_FLUX_H][i] = n_x[i]*u_in[F31_AVG][i] + n_y[i]*u_in[F32_AVG][i] +
			0.5 * lambda * u_in[H_JUMP][i];
		u_ex[NUM_FLUX_H][i_ex] = -u_in[NUM_FLUX_H][i];
	}
}