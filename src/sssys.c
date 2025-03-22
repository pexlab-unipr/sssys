#include "sssys.h"
#include "detego.h"

void statespace_init(Statespace* ss, int Nx, int Nu, int Ny,
	float* x, float* y, float* A, float* B, float* C, float* D)
{
	ss->Nx = Nx;
	ss->Nu = Nu;
	ss->Ny = Ny;
	ss->x = x;
	ss->y = y;
	ss->A = A;
	ss->B = B;
	ss->C = C;
	ss->D = D;
}

int statespace_tf2ss(Statespace* ss, float* num, float* den)
{
	int i;
	const int Nx = ss->Nx;
	float d0 = den[0], n0 = num[0] / d0;

	if (ss->Nu != 1 || ss->Ny != 1) {
		return 1;
	}

	for (i = 0; i < Nx; i++) {
		ss->A[i] = -den[i + 1] / d0;
		ss->B[i] = !i;
		ss->C[i] = num[i + 1] / d0 - n0 * den[i + 1] / d0;
	}

	for (i = Nx; i < Nx * Nx; i++) {
		ss->A[i] = (float)((i % (Nx + 1)) == Nx);
	}

	ss->D[0] = n0;

	return 0;
}

int statespace_c2d(Statespace* ss, float Ts, char* method, float* work)
{
	int i, j, k, p;
	const int Nx = ss->Nx;
	const int Nu = ss->Nu;
	const int Ny = ss->Ny;
	const int Nxx = Nx * Nx;
	const int Nxu = Nx * Nu;
	float* A = ss->A;
	float* B = ss->B;
	float* C = ss->C;
	float* D = ss->D;
	Matrixf L = { 0 };
	Matrixf R = { 0 };
	Matrixf Ad = { 0 };
	Matrixf Bd = { 0 };
	Matrixf Ad1 = { 0 };
	Matrixf Ad2 = { 0 };

	switch (*method)
	{
	case 'e': // Forward Euler

		// Ad = I+A*Ts
		// Bd = B*Ts
		// Cd = C
		// Dd = D

		// Ad = I+A*Ts
		for (i = 0; i < Nxx; i++) {
			A[i] = A[i] * Ts + !(i % (Nx + 1));
		}

		// Bd = B*Ts
		for (i = 0; i < Nxu; i++) {
			B[i] *= Ts;
		}

		break;

	case 'b': // Backward Euler

		// Ad = (I-A*Ts)^(-1) 
		// Bd = (I-A*Ts)^(-1)*B*Ts
		// Cd = C*(I-A*Ts)^(-1) = C*Ad
		// Dd = C*(I-A*Ts)^(-1)*B*Ts+D = C*Bd+D

		// R = [I, B*Ts]
		for (k = j = 0; j < Nx; j++) {
			for (i = 0; i < Nx; i++) {
				work[k++] = (float)(i == j);
			}
		}
		for (j = 0; j < Nu; j++) {
			for (i = 0; i < Nx; i++) {
				work[k++] = B[i * Nu + j] * Ts;
			}
		}
		matrixf_init(&R, Nx, Nx + Nu, work, 0);

		// L = I-A*Ts
		for (i = 0; i < Nxx; i++) {
			A[i] = A[i] * (-1.0f) * Ts + !(i % (Nx + 1));
		}
		matrixf_init(&L, Nx, Nx, A, 1);

		// Solve L*X = R for X, where X = [Ad, Bd]
		if (matrixf_solve_lu(&L, &R)) {
			return 1;
		}
		matrixf_init(&Ad, Nx, Nx, work, 0);
		matrixf_init(&Bd, Nx, Nu, work + Nxx, 0);
		matrixf_transpose(&Ad);
		matrixf_transpose(&Bd);
		for (i = 0; i < Nxx; i++) {
			A[i] = Ad.data[i];
		}
		for (i = 0; i < Nxu; i++) {
			B[i] = Bd.data[i];
		}

		// Dd = C*Bd+D
		for (i = 0; i < Ny; i++) {
			for (j = 0; j < Nu; j++) {
				for (k = 0; k < Nx; k++) {
					D[i * Nu + j] += C[i * Nx + k] * B[k * Nu + j];
				}
			}
		}

		// Cd = C*Ad
		for (i = 0; i < Ny; i++) {
			for (j = 0; j < Nx; j++) {
				work[j + Nxx] = 0;
				for (k = 0; k < Nx; k++) {
					work[j + Nxx] += C[i * Nx + k] * Ad.data[k * Nx + j];
				}
			}
			for (j = 0; j < Nx; j++) {
				C[i * Nx + j] = work[j + Nxx];
			}
		}

		break;

	case 't': // Tustin

		// Ad = (I-A*Ts/2)^(-1)*(I+A*Ts/2) 
		// Bd = (I-A*Ts/2)^(-1)*B*Ts
		// Cd = C*(I-A*Ts/2)^(-1)
		// Dd = C*(I-A*Ts/2)^(-1)*B*Ts/2+D = C*Bd/2+D

		// R = [I, A*Ts/2, B*Ts]
		for (k = j = 0; j < Nx; j++) {
			for (i = 0; i < Nx; i++) {
				work[k++] = (float)(i == j);
			}
		}
		for (j = 0; j < Nx; j++) {
			for (i = 0; i < Nx; i++) {
				work[k++] = A[i * Nx + j] * 0.5f * Ts;
			}
		}
		for (j = 0; j < Nu; j++) {
			for (i = 0; i < Nx; i++) {
				work[k++] = B[i * Nu + j] * Ts;
			}
		}
		matrixf_init(&R, Nx, Nx + Nx + Nu, work, 0);

		// L = I-A*Ts/2
		for (i = 0; i < Nxx; i++) {
			A[i] = A[i] * (-0.5f) * Ts + !(i % (Nx + 1));
		}
		matrixf_init(&L, Nx, Nx, A, 1);

		// Solve L*X = R for X, where X = [Ad1, Ad2, Bd] and Ad = Ad1 + Ad2
		if (matrixf_solve_lu(&L, &R)) {
			return 1;
		}
		matrixf_init(&Ad1, Nx, Nx, work, 0);
		matrixf_init(&Ad2, Nx, Nx, work + Nxx, 0);
		matrixf_init(&Bd, Nx, Nu, work + Nxx + Nxx, 0);
		matrixf_transpose(&Ad1);
		matrixf_transpose(&Ad2);
		matrixf_transpose(&Bd);
		for (i = 0; i < Nxx; i++) {
			A[i] = Ad1.data[i] + Ad2.data[i];
		}
		for (i = 0; i < Nxu; i++) {
			B[i] = Bd.data[i];
		}

		// Dd = C*Bd/2+D
		for (i = 0; i < Ny; i++) {
			for (j = 0; j < Nu; j++) {
				for (k = 0; k < Nx; k++) {
					D[i * Nu + j] += C[i * Nx + k] * B[k * Nu + j] * 0.5f;
				}
			}
		}

		// Cd = C*Ad1
		for (i = 0; i < Ny; i++) {
			for (j = 0; j < Nx; j++) {
				work[j + Nxx] = 0;
				for (k = 0; k < Nx; k++) {
					work[j + Nxx] += C[i * Nx + k] * Ad1.data[k * Nx + j];
				}
			}
			for (j = 0; j < Nx; j++) {
				C[i * Nx + j] = work[j + Nxx];
			}
		}

		break;

	case 'z': // Zero-order hold

		// Ad = exp(A*Ts)
		//      / Ts
		// Bd = | (exp(A*(Ts-t))dt)*B
		//      / 0
		// or [Ad, Bd; 0, I] = exp([A*Ts, B*Ts; 0, 0])
		// Cd = C
		// Dd = D

		// Assemble L = [A*Ts, B*Ts; 0, 0]
		matrixf_init(&L, Nx + Nu, Nx + Nu, work, 0);
		for (i = 0; i < L.size[0]; i++) {
			for (j = 0; j < L.size[1]; j++) {
				if (j < Nx) {
					at(&L, i, j) = (i < Nx) ? A[i * Nx + j] * Ts : 0;
				}
				else {
					at(&L, i, j) = (i < Nx) ? B[i * Nu + j - Nx] * Ts : 0;
				}
			}
		}

		// Compute exp(L) = [Ad, Bd; 0, I]
		if (matrixf_exp(&L, work + L.size[0] * L.size[1])) {
			return 1;
		}

		// Get the blocks Ad and Bd
		for (i = 0; i < Nx; i++) {
			for (j = 0; j < L.size[1]; j++) {
				if (j < Nx) {
					A[i * Nx + j] = at(&L, i, j);
				}
				else {
					B[i * Nu + j - Nx] = at(&L, i, j);
				}
			}
		}

		break;

	case 'f': // First-order hold

		// [Ad, G1, G2; 0, I, I; 0, 0, I] = exp([A*Ts, B*Ts, 0; 0, 0, I; 0, 0, 0])
		// Bd = G1+Ad*G2-G2;
		// Cd = C
		// Dd = D+C*G2

		// Assemble L = [A*Ts, B*Ts, 0; 0, 0, I; 0, 0, 0]
		matrixf_init(&L, Nx + 2 * Nu, Nx + 2 * Nu, work, 0);
		for (i = 0; i < L.size[0] * L.size[1]; i++) {
			work[i] = 0;
		}
		for (i = 0; i < Nu; i++) {
			at(&L, Nx + i, Nx + Nu + i) = 1;
		}
		for (i = 0; i < Nx; i++) {
			for (j = 0; j < Nx + Nu; j++) {
				if (j < Nx) {
					at(&L, i, j) = A[i * Nx + j] * Ts;
				}
				else {
					at(&L, i, j) = B[i * Nu + j - Nx] * Ts;
				}
			}
		}

		// Compute exp(L) = [Ad, G1, G2; 0, I, I; 0, 0, I]
		p = L.size[0] * L.size[1];
		if (matrixf_exp(&L, work + p)) {
			return 1;
		}

		// Get the blocks Ad, G1 and G2
		for (i = 0; i < Nx; i++) {
			for (j = 0; j < L.size[1]; j++) {
				if (j < Nx) {
					A[i * Nx + j] = at(&L, i, j);
				}
				else if (j < Nx + Nu) {
					B[i * Nu + j - Nx] = at(&L, i, j);
				}
				else {
					work[i * Nu + j - Nx - Nu + p] = at(&L, i, j);
				}
			}
		}

		// Bd = B+Ad*G2-G2
		for (i = 0; i < Nx; i++) {
			for (j = 0; j < Nu; j++) {
				B[i * Nu + j] -= work[i * Nu + j + p];
				for (k = 0; k < Nx; k++) {
					B[i * Nu + j] += A[i * Nx + k] * work[k * Nu + j + p];
				}
			}
		}

		// Dd = D+C*G2
		for (i = 0; i < Ny; i++) {
			for (j = 0; j < Nu; j++) {
				for (k = 0; k < Nx; k++) {
					D[i * Nu + j] += C[i * Nx + k] * work[k * Nu + j + p];
				}
			}
		}

		break;
	}

	return 0;
}

void statespace_iter(Statespace* ss, float* u, float* work)
{
	int i, j;
	const int Nx = ss->Nx;
	const int Nu = ss->Nu;
	const int Ny = ss->Ny;
	float* A = ss->A;
	float* B = ss->B;
	float* C = ss->C;
	float* D = ss->D;
	float tmp;

	for (i = 0; i < Nx; i++) {
		work[i] = ss->x[i];
	}

	// y[k] = Cd*x[k] + Dd*u[k]
	for (i = 0; i < Ny; i++) {
		tmp = 0;
		for (j = 0; j < Nx; j++) {
			tmp += C[i * Nx + j] * work[j];
		}
		for (j = 0; j < Nu; j++) {
			tmp += D[i * Nu + j] * u[j];
		}
		ss->y[i] = tmp;
	}

	// x[k+1] = Ad*x[k] + Bd*u[k]
	for (i = 0; i < Nx; i++) {
		tmp = 0;
		for (j = 0; j < Nx; j++) {
			tmp += A[i * Nx + j] * work[j];
		}
		for (j = 0; j < Nu; j++) {
			tmp += B[i * Nu + j] * u[j];
		}
		ss->x[i] = tmp;
	}
}