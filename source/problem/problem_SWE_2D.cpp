#include "../class_mesh_v2.h"
#include "problem_SWE_2D.h"

void problem_timestep(MESH* mesh) {
	std::vector <INTERFACE*> internal_interfaces(mesh->interfaces.find(INTERNAL_BOUNDARY)->second);

	for (auto it = mesh->elements.begin(); it != mesh->elements.end(); it++) {
		it->second->ComputeInternalU(H);
		it->second->ComputeBoundaryU(H);

		it->second->ComputeInternalU(UA);
		it->second->ComputeBoundaryU(UA);

		it->second->ComputeInternalU(VA);
		it->second->ComputeBoundaryU(VA);

		it->second->ComputeF();
	}

	for (auto it = internal_interfaces.begin(); it != internal_interfaces.end(); it++) {
		(*it)->ComputeAverageU(H, H_AVG);
		(*it)->ComputeAverageU(UA, UA_AVG);
		(*it)->ComputeAverageU(VA, VA_AVG);
	}


}

void problem_compute_f(int number_gp, double** U) {
	for (int i = 0; i < number_gp; i++) {
		U[F11][i] = U[SP][i] * (U[UA][i] * U[UA][i] / (U[ZB][i] + U[H][i]) 
			+ 0.5*GRAVITY*(U[ZB][i] + U[H][i]) * (U[ZB][i] + U[H][i]));
		U[F21][i] = U[SP][i] * U[UA][i] * U[VA][i] / (U[ZB][i] + U[H][i]);
		U[F31][i] = U[SP][i] * U[UA][i];

		U[F12][i] = U[UA][i] * U[VA][i] / (U[ZB][i] + U[H][i]);
		U[F22][i] = (U[VA][i] * U[VA][i] / (U[ZB][i] + U[H][i]) 
			+ 0.5*GRAVITY*(U[ZB][i] + U[H][i]) * (U[ZB][i] + U[H][i]));
		U[F32][i] = U[VA][i];
	}
}