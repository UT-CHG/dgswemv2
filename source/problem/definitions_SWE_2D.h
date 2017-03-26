#ifndef SWE_2D
#define SWE_2D

#define GRAVITY 9.81
#define VISCOSITY 0

#define SIZE_U 9U
//VARIABLES TO STORE IN U
#define SP 0U
#define ZB 1U
#define UA 2U
#define VA 3U
#define H 4U
#define TXX 5U
#define TXY 6U
#define TYX 7U
#define TYY 8U

#define D_UA 5U
#define D_VA 6U
#define D_H 7U

#define SIZE_U_INTERNAL 14
//INTERNAL VARIABLES TO BE CALCULATED AND STORED AT TIMESTEP:
//VARIABLES STORED IN U +
#define F11 5U
#define F12 6U
#define F21 7U
#define F22 8U
#define F31 9U
#define F32 10U
#define U 11U
#define V 12U
#define A 13U

#define SIZE_U_BOUNDARY 14
//BOUNDARY VARIABLES TO BE CALCULATED AND STORED AT TIMESTEP
//VARIABLES STORED IN U_INTERNAL +
//BOUNDARY FLAGS DIFFUSION STEP
#define UA_AVG 11U
#define VA_AVG 12U
#define H_AVG 13U

//BOUNDARY FLAGS HYPER STEP
#define F11_AVG 5U
#define F12_AVG 6U
#define F21_AVG 7U
#define F22_AVG 8U
#define F31_AVG 9U
#define F32_AVG 10U
#define UA_JUMP 11U
#define VA_JUMP 12U
#define H_JUMP 13U

#define NUM_FLUX_UA 11U
#define NUM_FLUX_VA 12U
#define NUM_FLUX_H 13U

#define INTERNAL 0U
#define OCEAN 1U
#define LAND 2U
#define FLOW 3U

//TIMESTEPPING
#define U0 0U
#define K1 1U
#define K2 2U
#define K3 3U
#define K4 4U

#endif
