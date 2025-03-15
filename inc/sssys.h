#ifndef SSSYS_H
#define SSSYS_H

typedef struct {
	int Nx;						// number of states
	int Nu;						// number of inputs
	int Ny;						// number of outputs
	float* x;					// states
	float* y;					// outputs
	float* A, * B, * C, * D;	// system matrices
} State_space;

// Initialization of a state-space system. 
// The matrices are stored in row-major order.
void state_space_init(State_space* ss, int Nx, int Nu, int Ny,
	float* x, float* y, float* A, float* B, float* C, float* D);

// Conversion of SISO transfer function coefficients into state-space form.
// num and den are the arrays of numerator and denominator coefficients, 
// respectively. Their length is Nx + 1. The coefficients are specified in
// order of descending power.
int state_space_convert_tf_coeffs(State_space* ss, float* num, float* den);

// Discretization of a state-space system. 
// This function converts the matrices from continuous to discrete form
// with sample time Ts. The available methods are:
// 1) "euler" (forward Euler);
// 2) "backward" (backward Euler);
// 3) "tustin" (Tustin);
// 4) "zoh" (zero-order hold);
// 5) "foh" (first-order hold).
// The length of the workspace (work) depends on the selected method:
// 1) "euler" -> 0;
// 2) "backward" -> Nx*(Nx+Nu);
// 3) "tustin" -> Nx*(2*Nx+Nu);
// 4) "zoh" -> (Nx+Nu)*(4*(Nx+Nu)+1);
// 5) "foh" -> (Nx+2*Nu)*(4*(Nx+2*Nu)+1).
int state_space_discretize(State_space* ss, float Ts, char* method, float* work);

// Update of the discrete system response with inputs u. 
// The length of the workspace (work) is Nx.
void state_space_iter(State_space* ss, float* u, float* work);

#endif