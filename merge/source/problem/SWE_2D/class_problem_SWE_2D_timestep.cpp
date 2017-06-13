#include <cmath>
#include <algorithm> 

#include "class_problem_SWE_2D.h"

void PROBLEM::Timestep() {
	for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
		it->second->ComputeInternalU(UA);
		it->second->ComputeBoundaryU(UA);

		it->second->ComputeInternalU(VA);
		it->second->ComputeBoundaryU(VA);

		it->second->ComputeInternalU(H);
		it->second->ComputeBoundaryU(H);
	}

	for (auto it = this->mesh->interfaces.begin(); it != this->mesh->interfaces.end(); it++) {
		switch (it->first) {
		case INTERNAL:
			//for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
			//	this->InterfaceFlowAverage(*(itt));
			//}
			break;
		case OCEAN:
			for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
				this->OceanInterfaceSetBC(*(itt));
				//this->InterfaceFlowAverage(*(itt));
			}
			break;
		case LAND:
			for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
				this->LandInterfaceSetBC(*(itt));
				//this->InterfaceFlowAverage(*(itt));
			}
			break;
		case FLOW:
			for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
				this->FlowInterfaceSetBC(*(itt));
				//this->InterfaceFlowAverage(*(itt));
			}
			break;
		default:
			printf("\n");
			printf("PROBLEM Timestep - Fatal error!\n");
			printf("Undefined interface type = %d\n", it->first);
			exit(1);
		}
	}

	int number_bf;
	double* RHS;

	//for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
	//	number_bf = it->second->number_bf;
	//	RHS = it->second->RHS;

	//	//COMPUTING TXX
	//	for (int i = 0; i < number_bf; i++) {
	//		RHS[i] = (-it->second->IntegrationInternalDPhiDX(UA, i) +
	//			it->second->IntegrationBoundaryNX(UA_AVG, i))*
	//			pow(it->second->u[SP][i], 2);
	//	}
	//	it->second->SolveLSE(TXX);

	//	it->second->ComputeInternalU(TXX);
	//	it->second->ComputeBoundaryU(TXX);

	//	//COMPUTING TXY
	//	for (int i = 0; i < number_bf; i++) {
	//		RHS[i] = -it->second->IntegrationInternalDPhiDY(UA, i) +
	//			it->second->IntegrationBoundaryNY(UA_AVG, i);
	//	}
	//	it->second->SolveLSE(TXY);

	//	it->second->ComputeInternalU(TXY);
	//	it->second->ComputeBoundaryU(TXY);

	//	//COMPUTING TYX
	//	for (int i = 0; i < number_bf; i++) {
	//		RHS[i] = (-it->second->IntegrationInternalDPhiDX(VA, i) +
	//			it->second->IntegrationBoundaryNX(VA_AVG, i))*
	//			pow(it->second->u[SP][i], 2);
	//	}
	//	it->second->SolveLSE(TYX);

	//	it->second->ComputeInternalU(TYX);
	//	it->second->ComputeBoundaryU(TYX);

	//	//COMPUTING TYY
	//	for (int i = 0; i < number_bf; i++) {
	//		RHS[i] = -it->second->IntegrationInternalDPhiDY(VA, i) +
	//			it->second->IntegrationBoundaryNY(VA_AVG, i);
	//	}
	//	it->second->SolveLSE(TYY);

	//	it->second->ComputeInternalU(TYY);
	//	it->second->ComputeBoundaryU(TYY);
	//}

	for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
		this->ComputeUVA(it->second);
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
		number_bf = it->second->number_bf;
        RHS = it->second->RHS;

        //COMPUTING UA
        for (int i = 0; i < number_bf; i++) {
			RHS[i] = it->second->IntegrationInternalDPhi(X, F11, i) +
				it->second->IntegrationInternalDPhi(Y, F12, i) -
				it->second->IntegrationBoundaryPhi(NUM_FLUX_UA, i) -
				it->second->IntegrationInternalPhi(S1, i);
        }
        it->second->SolveLSE(D_UA);

        //COMPUTING VA
        for (int i = 0; i < number_bf; i++) {
            RHS[i] = it->second->IntegrationInternalDPhi(X, F21, i) +
                it->second->IntegrationInternalDPhi(Y, F22, i) -
                it->second->IntegrationBoundaryPhi(NUM_FLUX_VA, i) -
				it->second->IntegrationInternalPhi(S2, i);;
        }

        it->second->SolveLSE(D_VA);

        //COMPUTING H
        for (int i = 0; i < number_bf; i++) {
            RHS[i] = it->second->IntegrationInternalDPhi(X, F31, i) +
                it->second->IntegrationInternalDPhi(Y, F32, i) -
                it->second->IntegrationBoundaryPhi(NUM_FLUX_H, i) - 
				it->second->IntegrationInternalPhi(S3, i);
        }
        it->second->SolveLSE(D_H);
    }
}

void PROBLEM::InterfaceFlowAverage(INTERFACE* intface) {
    double** u_in = intface->u_boundary_in;
    double** u_ex = intface->u_boundary_ex;

    int i_ex;

    for (int i = 0; i < intface->number_gp; i++) {
        i_ex = intface->number_gp - 1 - i;

        u_in[UA_AVG][i] = 0.5 * (u_in[UA][i] + u_ex[UA][i_ex]);
        u_ex[UA_AVG][i_ex] = u_in[UA_AVG][i];

        u_in[VA_AVG][i] = 0.5 * (u_in[VA][i] + u_ex[VA][i_ex]);
        u_ex[VA_AVG][i_ex] = u_in[VA_AVG][i];
    }
}

void PROBLEM::ComputeUVA(ELEMENT* element) {
    double** u = element->u_internal;

    for (int i = 0; i < element->number_gp_internal; i++) {
        u[A][i] = u[ZB][i] + u[H][i];
        u[U][i] = u[UA][i] / u[A][i];
        u[V][i] = u[VA][i] / u[A][i];
    }

    for (int j = 0; j < element->number_boundaries; j++) {
        u = element->u_boundary[j];

        for (int i = 0; i < element->number_gp_boundary; i++) {
            u[A][i] = u[ZB][i] + u[H][i];
            u[U][i] = u[UA][i] / u[A][i];
            u[V][i] = u[VA][i] / u[A][i];
        }
    }
}

void PROBLEM::ComputeF(ELEMENT* element) {
    double** u = element->u_internal;

    for (int i = 0; i < element->number_gp_internal; i++) {
		u[F11][i] = u[SP][i] * (u[UA][i] * u[U][i] + 0.5*GRAVITY*pow(u[A][i], 2));// -VISCOSITY*u[TXX][i];
		u[F21][i] = u[SP][i] * u[UA][i] * u[V][i];// -VISCOSITY*u[TYX][i];
        u[F31][i] = u[SP][i] * u[UA][i];

		u[F12][i] = u[UA][i] * u[V][i];// -VISCOSITY*u[TXY][i];
		u[F22][i] = u[VA][i] * u[V][i] + 0.5*GRAVITY*pow(u[A][i], 2);// -VISCOSITY*u[TYY][i];
        u[F32][i] = u[VA][i];
    }

    for (int j = 0; j < element->number_boundaries; j++) {
        u = element->u_boundary[j];

        for (int i = 0; i < element->number_gp_boundary; i++) {
			u[F11][i] = u[SP][i] * (u[UA][i] * u[U][i] + 0.5*GRAVITY*pow(u[A][i], 2));// -VISCOSITY*u[TXX][i];
			u[F21][i] = u[SP][i] * u[UA][i] * u[V][i];// -VISCOSITY*u[TYX][i];
            u[F31][i] = u[SP][i] * u[UA][i];

			u[F12][i] = u[UA][i] * u[V][i];// -VISCOSITY*u[TXY][i];
			u[F22][i] = u[VA][i] * u[V][i] + 0.5*GRAVITY*pow(u[A][i], 2);// -VISCOSITY*u[TYY][i];
            u[F32][i] = u[VA][i];
		
			//printf("%f %f\n%f %f\n%f %f\n\n", u[F11][i], u[F12][i], u[F21][i], u[F22][i], u[F31][i], u[F32][i]);
		}
    }
}

void PROBLEM::ComputeS(ELEMENT* element) {
	double** u = element->u_internal;

	double cf = 0; // FOR TESTING NO BOTTOM FRICTION 
	double tau_b;
	
	double f = 0; // FOR TESTING NO CORIOLIS
	
	for (int i = 0; i < element->number_gp_internal; i++) {
		tau_b = 0;//cf / sqrt(pow(u[U][i], 2) + pow(u[V][i], 2));

		u[S1][i] = 0;// f*u[V][i] * u[A][i]
			//- tau_b*u[U][i];

		u[S2][i] = 0;// -f*u[U][i] * u[A][i]
			//- tau_b*u[V][i];

		u[S3][i] = 0;
	}	
	
	//element->ComputeInternalDUDX(ZB, DZBDX);
	//element->ComputeInternalDUDY(ZB, DZBDY);
	//
	//for (int i = 0; i < element->number_gp_internal; i++) {
	//	tau_b = 0;//cf / sqrt(pow(u[U][i], 2) + pow(u[V][i], 2));

	//	u[S1][i] = u[S1][i]
	//		- u[SP][i] * GRAVITY*u[A][i] * u[DZBDX][i];

	//	u[S2][i] = u[S2][i]
	//		- GRAVITY*u[A][i] * u[DZBDY][i];
	//}
}

void PROBLEM::LLFNumericalFlux(INTERFACE* intface) {
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

    double** u_in = intface->u_boundary_in;
    double** u_ex = intface->u_boundary_ex;
    double* n_x = intface->normal[X];
    double* n_y = intface->normal[Y];

    int i_ex;

    for (int i = 0; i < intface->number_gp; i++) {
        i_ex = intface->number_gp - 1 - i;

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

        u_in[NUM_FLUX_UA][i] = n_x[i] * u_in[F11_AVG][i] + n_y[i] * u_in[F12_AVG][i] +
            0.5 * lambda * u_in[UA_JUMP][i];
        u_ex[NUM_FLUX_UA][i_ex] = -u_in[NUM_FLUX_UA][i];

        u_in[NUM_FLUX_VA][i] = n_x[i] * u_in[F21_AVG][i] + n_y[i] * u_in[F22_AVG][i] +
            0.5 * lambda * u_in[VA_JUMP][i];
        u_ex[NUM_FLUX_VA][i_ex] = -u_in[NUM_FLUX_VA][i];

        u_in[NUM_FLUX_H][i] = n_x[i] * u_in[F31_AVG][i] + n_y[i] * u_in[F32_AVG][i] +
            0.5 * lambda * u_in[H_JUMP][i];
        u_ex[NUM_FLUX_H][i_ex] = -u_in[NUM_FLUX_H][i];
    }
}

void PROBLEM::OceanInterfaceSetBC(INTERFACE* intface) {
	double H_0 = 0.3;

	if (this->t < 172800.0) H_0 = 0.3*this->t / 172800.0; //LINEAR RAMPING

	double H_ocean = H_0*cos(2 * PI*this->t / 43200); //FOR TESTING M2 TIDAL WAVE WITH PERIOD OF 12HOURS AND AMPLITUDE OF 0.3m

    double** u_in = intface->u_boundary_in;
    double** u_ex = intface->u_boundary_ex;

    int i_ex;

    for (int i = 0; i < intface->number_gp; i++) {
        i_ex = intface->number_gp - 1 - i;

        //H_ocean = u_in[H][i]; //NEED TO COMPUTE OCEAN ELEVATION

        u_ex[ZB][i_ex] = u_in[ZB][i];
        u_ex[SP][i_ex] = u_in[SP][i];
        u_ex[UA][i_ex] = u_in[UA][i];
        u_ex[VA][i_ex] = u_in[VA][i];
        u_ex[H][i_ex] = H_ocean; //BOUNDARY CONDITION SET
    }
}

void PROBLEM::LandInterfaceSetBC(INTERFACE* intface) {
    double n_x;
    double n_y;
    double t_x;
    double t_y;

    double qn_in;
    double qt_in;
    double qn_ex;
    double qt_ex;

    double** u_in = intface->u_boundary_in;
    double** u_ex = intface->u_boundary_ex;

    int i_ex;

    for (int i = 0; i < intface->number_gp; i++) {
        i_ex = intface->number_gp - 1 - i;

        n_x = intface->normal[X][i];
        n_y = intface->normal[Y][i];
        t_x = -n_y;
        t_y = n_x;

        qn_in = u_in[UA][i] * n_x + u_in[VA][i] * n_y;
        qt_in = u_in[UA][i] * t_x + u_in[VA][i] * t_y;

        qn_ex = -qn_in;
        qt_ex = qt_in;

        u_ex[ZB][i_ex] = u_in[ZB][i];
        u_ex[SP][i_ex] = u_in[SP][i];
        u_ex[UA][i_ex] = qn_ex*n_x + qt_ex*t_x;
        u_ex[VA][i_ex] = qn_ex*n_y + qt_ex*t_y;
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
        u_ex[UA][i_ex] = qn_ex*n_x;
        u_ex[VA][i_ex] = qn_ex*n_y;
        u_ex[H][i_ex] = u_in[H][i];
    }
}

void PROBLEM::ComputeBoundaryInterfaceF(INTERFACE* intface) {
    double** u_in = intface->u_boundary_in;
    double** u_ex = intface->u_boundary_ex;

    int i_in;

    for (int i = 0; i < intface->number_gp; i++) {
        i_in = intface->number_gp - 1 - i;

        u_ex[A][i] = u_ex[ZB][i] + u_ex[H][i];
        u_ex[U][i] = u_ex[UA][i] / u_ex[A][i];
        u_ex[V][i] = u_ex[VA][i] / u_ex[A][i];

        u_ex[F11][i] = u_ex[SP][i] * (u_ex[UA][i] * u_ex[U][i] + 0.5*GRAVITY*pow(u_ex[A][i], 2)) - VISCOSITY*u_in[TXX][i_in];
        u_ex[F21][i] = u_ex[SP][i] * u_ex[UA][i] * u_ex[V][i] - VISCOSITY*u_in[TYX][i_in];
        u_ex[F31][i] = u_ex[SP][i] * u_ex[UA][i];

        u_ex[F12][i] = u_ex[UA][i] * u_ex[V][i] - VISCOSITY*u_in[TXY][i_in];
        u_ex[F22][i] = u_ex[VA][i] * u_ex[V][i] + 0.5*GRAVITY*pow(u_ex[A][i], 2) - VISCOSITY*u_in[TYY][i_in];
        u_ex[F32][i] = u_ex[VA][i];
    }
}