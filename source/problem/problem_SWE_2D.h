#ifndef SWE_2D
#define SWE_2D

#define GRAVITY 9.81

#define SIZE_U 15U

#define UA 0U
#define VA 1U
#define H 2U
#define ZB 3U
#define SP 4U

#define UA_AVG 5U
#define VA_AVG 6U
#define H_AVG 7U

#define TXX 5U
#define TXY 6U
#define TYX 7U
#define TYY 8U

#define F11 9U
#define F12 10U
#define F21 11U
#define F22 12U
#define F31 13U
#define F32 14U

#define INTERNAL_BOUNDARY 0U

void problem_compute_f(int number_gp, double** U);

#endif
