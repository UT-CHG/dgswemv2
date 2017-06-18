#include <cmath>
#include <algorithm> 

#include "class_problem_SWE_2D_LINEAR.h"

void PROBLEM::Timestep() {
	for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
		it->second->ComputeInternalU(U);
		it->second->ComputeBoundaryU(U);

		it->second->ComputeInternalU(V);
		it->second->ComputeBoundaryU(V);

		it->second->ComputeInternalU(H);
		it->second->ComputeBoundaryU(H);
	}

	for (auto it = this->mesh->interfaces.begin(); it != this->mesh->interfaces.end(); it++) {
		switch (it->first) {
		case INTERNAL:
			break;
		case OCEAN:
			for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
				this->OceanInterfaceSetBC(*(itt));
			}
			break;
		case LAND:	
			for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
				this->LandInterfaceSetBC(*(itt));
			}
			break;
		case FLOW:	
			for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
				this->FlowInterfaceSetBC(*(itt));
			}
			break;
		default:
			printf("\n");
			printf("PROBLEM Timestep - Fatal error!\n");
			printf("Undefined interface type = %d\n", it->first);
			exit(1);
		}
	}

	for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
		this->ComputeF(it->second);
		this->ComputeS(it->second);
	}

	for (auto it = this->mesh->interfaces.begin(); it != this->mesh->interfaces.end(); it++) {
		switch (it->first) {
		case INTERNAL:
			for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
				this->LLFNumericalFlux(*(itt));
			}
			break;
		case OCEAN:
		case LAND:
		case FLOW:
			for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
				this->ComputeBoundaryInterfaceF(*(itt));
				this->LLFNumericalFlux(*(itt));
			}
			break;
		default:
			printf("\n");
			printf("PROBLEM Timestep - Fatal error!\n");
			printf("Undefined interface type = %d\n", it->first);
			exit(1);
		}
	}

	for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
		int number_bf = it->second->number_bf;
        std::vector<double>& RHS = it->second->RHS;

        //COMPUTING U
        for (int i = 0; i < number_bf; i++) {
			RHS[i] = it->second->IntegrationInternalDPhi(X, F11, i) +
				it->second->IntegrationInternalDPhi(Y, F12, i) -
				it->second->IntegrationBoundaryPhi(NUM_FLUX_U, i) +
				it->second->IntegrationInternalPhi(S1, i);
        }
        it->second->SolveLSE(D_U);

        //COMPUTING V
        for (int i = 0; i < number_bf; i++) {
            RHS[i] = it->second->IntegrationInternalDPhi(X, F21, i) +
                it->second->IntegrationInternalDPhi(Y, F22, i) -
                it->second->IntegrationBoundaryPhi(NUM_FLUX_V, i) +
				it->second->IntegrationInternalPhi(S2, i);
        }

        it->second->SolveLSE(D_V);

        //COMPUTING H
        for (int i = 0; i < number_bf; i++) {
            RHS[i] = it->second->IntegrationInternalDPhi(X, F31, i) +
                it->second->IntegrationInternalDPhi(Y, F32, i) -
                it->second->IntegrationBoundaryPhi(NUM_FLUX_H, i) + 
				it->second->IntegrationInternalPhi(S3, i);
        }
        it->second->SolveLSE(D_H);
    }
}

void PROBLEM::ComputeF(Element<>* element) {
    Array2D<double>& u = element->u_gp;

    for (int i = 0; i < element->number_gp_internal; i++) {
		u[F11][i] = u[SP][i] * GRAVITY * u[H][i];
		u[F21][i] = 0;
		u[F31][i] = u[SP][i] * u[U][i] * u[ZB][i];

		u[F12][i] = 0;
		u[F22][i] = GRAVITY * u[H][i];
		u[F32][i] = u[V][i] * u[ZB][i];
	}

    for (int j = 0; j < element->number_boundaries; j++) {
        double** u = element->u_boundary[j];

        for (int i = 0; i < element->number_gp_boundary; i++) {
			u[F11][i] = u[SP][i] * GRAVITY * u[H][i];
			u[F21][i] = 0;
			u[F31][i] = u[SP][i] * u[U][i] * u[ZB][i];

			u[F12][i] = 0;
			u[F22][i] = GRAVITY * u[H][i];
			u[F32][i] = u[V][i] * u[ZB][i];
		}
    }
}

void PROBLEM::ComputeS(Element<>* element) {
	Array2D<double>& u = element->u_gp;

	double cf = 0; // FOR TESTING NO BOTTOM FRICTION 
	double tau_b;
	
	double f = 0; // FOR TESTING NO CORIOLIS

	for (int i = 0; i < element->number_gp_internal; i++) {
		tau_b = 0;//cf / sqrt(pow(u[U][i], 2) + pow(u[V][i], 2));

		u[S1][i] = - tau_b*u[U][i];

		u[S2][i] = - tau_b*u[V][i];
		
		u[S3][i] = 0;
	}
}

void PROBLEM::LLFNumericalFlux(INTERFACE* intface) {
    double sqrt_g_a_in;
    double sqrt_g_a_ex;

    double lambda;

    double** u_in = intface->u_boundary_in;
    double** u_ex = intface->u_boundary_ex;
    double* n_x = intface->normal[X];
    double* n_y = intface->normal[Y];

    int i_ex;

    for (int i = 0; i < intface->number_gp; i++) {
        i_ex = intface->number_gp - 1 - i;

        sqrt_g_a_in = sqrt(GRAVITY * u_in[ZB][i] * (pow((n_x[i] * u_in[SP][i]), 2) + pow(n_y[i], 2)));

        sqrt_g_a_ex = sqrt(GRAVITY * u_ex[ZB][i_ex] * (pow((n_x[i_ex] * u_ex[SP][i_ex]), 2) + pow(n_y[i_ex], 2)));

        lambda = std::max(sqrt_g_a_in, sqrt_g_a_ex);

        u_in[F11_AVG][i] = 0.5*(u_in[F11][i] + u_ex[F11][i_ex]);
        u_in[F12_AVG][i] = 0.5*(u_in[F12][i] + u_ex[F12][i_ex]);
        u_in[F21_AVG][i] = 0.5*(u_in[F21][i] + u_ex[F21][i_ex]);
        u_in[F22_AVG][i] = 0.5*(u_in[F22][i] + u_ex[F22][i_ex]);
        u_in[F31_AVG][i] = 0.5*(u_in[F31][i] + u_ex[F31][i_ex]);
        u_in[F32_AVG][i] = 0.5*(u_in[F32][i] + u_ex[F32][i_ex]);

        u_in[U_JUMP][i] = u_in[U][i] - u_ex[U][i_ex];
        u_in[V_JUMP][i] = u_in[V][i] - u_ex[V][i_ex];
        u_in[H_JUMP][i] = u_in[H][i] - u_ex[H][i_ex];

		u_in[NUM_FLUX_U][i] = n_x[i] * u_in[F11_AVG][i] + n_y[i] * u_in[F12_AVG][i] +
			0.5 * lambda * (u_in[U_JUMP][i] * pow(n_x[i], 2) + u_in[V_JUMP][i] * n_x[i] * n_y[i]);
        u_ex[NUM_FLUX_U][i_ex] = -u_in[NUM_FLUX_U][i];

		u_in[NUM_FLUX_V][i] = n_x[i] * u_in[F21_AVG][i] + n_y[i] * u_in[F22_AVG][i] +
			0.5 * lambda * (u_in[V_JUMP][i] * pow(n_y[i], 2) + u_in[U_JUMP][i] * n_x[i] * n_y[i]);
        u_ex[NUM_FLUX_V][i_ex] = -u_in[NUM_FLUX_V][i];

        u_in[NUM_FLUX_H][i] = n_x[i] * u_in[F31_AVG][i] + n_y[i] * u_in[F32_AVG][i] +
            0.5 * lambda * u_in[H_JUMP][i];
        u_ex[NUM_FLUX_H][i_ex] = -u_in[NUM_FLUX_H][i];
    }
}

void PROBLEM::OceanInterfaceSetBC(INTERFACE* intface) {
	double H_0 = 0.3;

	//if (this->t < 172800.0) H_0 = 0.3*this->t / 172800.0; //LINEAR RAMPING

	double H_ocean = H_0*cos(2*PI*this->t/43200); //FOR TESTING M2 TIDAL WAVE WITH PERIOD OF 12HOURS AND AMPLITUDE OF 0.3m

    double** u_in = intface->u_boundary_in;
    double** u_ex = intface->u_boundary_ex;

    int i_ex;
    for (int i = 0; i < intface->number_gp; i++) {
        i_ex = intface->number_gp - 1 - i;

        //H_ocean = u_in[H][i]; //NEED TO COMPUTE OCEAN ELEVATION
		u_ex[ZB][i_ex] = u_in[ZB][i];
        u_ex[SP][i_ex] = u_in[SP][i];
        u_ex[U][i_ex] = u_in[U][i];
        u_ex[V][i_ex] = u_in[V][i];
        u_ex[H][i_ex] = H_ocean; //BOUNDARY CONDITION SET
    }
}

void PROBLEM::LandInterfaceSetBC(INTERFACE* intface) {
    double n_x;
    double n_y;
    double t_x;
    double t_y;

    double un_in;
    double ut_in;
    double un_ex;
    double ut_ex;

    double** u_in = intface->u_boundary_in;
    double** u_ex = intface->u_boundary_ex;

    int i_ex;

    for (int i = 0; i < intface->number_gp; i++) {
        i_ex = intface->number_gp - 1 - i;

        n_x = intface->normal[X][i];
        n_y = intface->normal[Y][i];
        t_x = -n_y;
        t_y = n_x;

        un_in = u_in[U][i] * n_x + u_in[V][i] * n_y;
        ut_in = u_in[U][i] * t_x + u_in[V][i] * t_y;

        un_ex = -un_in;
        ut_ex = ut_in;

        u_ex[ZB][i_ex] = u_in[ZB][i];
        u_ex[SP][i_ex] = u_in[SP][i];
        u_ex[U][i_ex] = un_ex*n_x + ut_ex*t_x;
        u_ex[V][i_ex] = un_ex*n_y + ut_ex*t_y;
        u_ex[H][i_ex] = u_in[H][i];
    }
}

void PROBLEM::FlowInterfaceSetBC(INTERFACE* intface) {
    double n_x;
    double n_y;
    double t_x;
    double t_y;

    double qn_ex;

    double** u_in = intface->u_boundary_in;
    double** u_ex = intface->u_boundary_ex;

    int i_ex;

    for (int i = 0; i < intface->number_gp; i++) {
        i_ex = intface->number_gp - 1 - i;

        n_x = intface->normal[X][i];
        n_y = intface->normal[Y][i];
        t_x = -n_y;
        t_y = n_x;

        qn_ex = 0; //NEED TO COMPUTE QN_EX

        u_ex[ZB][i_ex] = u_in[ZB][i];
        u_ex[SP][i_ex] = u_in[SP][i];
        u_ex[U][i_ex] = qn_ex*n_x;
        u_ex[V][i_ex] = qn_ex*n_y;
        u_ex[H][i_ex] = u_in[H][i];
    }
}

void PROBLEM::ComputeBoundaryInterfaceF(INTERFACE* intface) {
    double** u_ex = intface->u_boundary_ex;

    for (int i = 0; i < intface->number_gp; i++) {
		u_ex[F11][i] = u_ex[SP][i] * GRAVITY * u_ex[H][i];
        u_ex[F21][i] = 0;
		u_ex[F31][i] = u_ex[SP][i] * u_ex[U][i] * u_ex[ZB][i];

        u_ex[F12][i] = 0;
		u_ex[F22][i] = GRAVITY * u_ex[H][i];
		u_ex[F32][i] = u_ex[V][i] * u_ex[ZB][i];
    }
}