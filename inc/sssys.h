#ifndef SSSYS_H
#define SSSYS_H

typedef struct {
	int Nx; // number of states
	int Nu; // number of inputs
	int Ny; // number of outputs
	float* x; // states
	float* y; // outputs
	float* A, * B, * C, * D; // system matrices
} Statespace;

// Initialization of a state-space system. 
// The matrices are stored in row-major order.
void statespace_init(Statespace* ss, int Nx, int Nu, int Ny,
	float* x, float* y, float* A, float* B, float* C, float* D);

// Conversion of discrete or continuous-time SISO transfer function
// into state-space form. num and den are the arrays of transfer-function
// numerator and denominator coefficients, respectively. Their length is 
// Nx + 1. The coefficients are specified in order of descending power.
// The function returns 1 if the system is not SISO, otherwise, in case
// of success, it returns 0.
int statespace_tf2ss(Statespace* ss, const float* num, const float* den);

// Discretization of a continuous-time state-space system. 
// This function converts the matrices from continuous to discrete form
// with sample time Ts. The available methods are:
// 1) "e" (forward Euler);
// 2) "b" (backward Euler);
// 3) "t" (Tustin);
// 4) "z" (zero-order hold);
// 5) "f" (first-order hold).
// The length of the workspace (work) depends on the selected method:
// 1) "e" -> 0;
// 2) "b" -> Nx*(Nx+Nu);
// 3) "t" -> Nx*(2*Nx+Nu);
// 4) "z" -> (Nx+Nu)*(4*(Nx+Nu)+1);
// 5) "f" -> (Nx+2*Nu)*(4*(Nx+2*Nu)+1).
// Use the dedicated macros to set the appropriate workspace memory.
// The function returns 1 if the system is singular, otherwise, in case
// of success, it returns 0.
int statespace_c2d(Statespace* ss, const float Ts, const char* method, float* work);

// Update of the discrete system response with inputs u. 
// The length of the workspace (work) is Nx.
void statespace_iter(Statespace* ss, float* u, float* work);

// Macro for setting the appropriate workspace memory for backward Euler discretization.
#define STATESPACE_WORKLEN_BWE(Nx, Nu) ((Nx)*((Nx)+(Nu)))
// Macro for setting the appropriate workspace memory for Tustin discretization.
#define STATESPACE_WORKLEN_TUS(Nx, Nu) ((Nx)*(2*(Nx)+(Nu)))
// Macro for setting the appropriate workspace memory for zero-order hold discretization.
#define STATESPACE_WORKLEN_ZOH(Nx, Nu) (((Nx)+(Nu))*(4*((Nx)+(Nu))+1))
// Macro for setting the appropriate workspace memory for first-order hold discretization.
#define STATESPACE_WORKLEN_FOH(Nx, Nu) (((Nx)+2*(Nu))*(4*((Nx)+2*(Nu))+1))

#endif