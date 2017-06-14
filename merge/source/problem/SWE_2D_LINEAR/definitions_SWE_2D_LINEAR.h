#ifndef SWE_2D_LINEAR
#define SWE_2D_LINEAR

#define GRAVITY 9.81
#define VISCOSITY 0

#define SIZE_U 8U
//VARIABLES TO STORE IN U
#define SP 0U
#define ZB 1U
#define U 2U
#define V 3U
#define H 4U

#define D_U 5U
#define D_V 6U
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

#define DZBDX 2U
#define DZBDY 3U

#define S1 11U
#define S2 12U
#define S3 13U

#define SIZE_U_BOUNDARY 14
//BOUNDARY VARIABLES TO BE CALCULATED AND STORED AT TIMESTEP
//VARIABLES STORED IN U_INTERNAL +
//BOUNDARY FLAGS DIFFUSION STEP
#define U_AVG 11U
#define V_AVG 12U
#define H_AVG 13U

//BOUNDARY FLAGS HYPER STEP
#define F11_AVG 5U
#define F12_AVG 6U
#define F21_AVG 7U
#define F22_AVG 8U
#define F31_AVG 9U
#define F32_AVG 10U
#define U_JUMP 11U
#define V_JUMP 12U
#define H_JUMP 13U

#define NUM_FLUX_U 5U
#define NUM_FLUX_V 6U
#define NUM_FLUX_H 7U

//INTERFACE TYPES
#define INTERNAL 0U
#define OCEAN 1U
#define LAND 2U
#define FLOW 3U

#endif