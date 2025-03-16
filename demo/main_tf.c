/*This example shows how to run the continuous-time transfer function

	6e-06 s^4 + 0.004987 s^3 + 0.987 s^2 + 394.8 s + 3.948e-07
H = ----------------------------------------------------------
			0.0001 s^4 + 0.04 s^3 + 9.87 s^2 + 3948 s

Use the following script to check the results.

clear; close all; clc;

s = tf('s');

sys = ...
	(6e-06*s^4 + 0.004987*s^3 + 0.987*s^2 + 394.8*s + 3.948e-07)/...
	(0.0001*s^4 + 0.04*s^3 + 9.87*s^2 + 3948*s);

t = 0:1e-5:2;
u = sin(2*pi*40*t);
y = lsim(sys, u, t);

p = "./";
fid = fopen(p + "t.bin");
t_ = fread(fid,'single');
fclose(fid);

fid = fopen(p + "y.bin");
y_ = fread(fid,'single');
fclose(fid);

figure; plot(t,y,t_,y_)
*/
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include "sssys.h"

#define Nx 4

int main()
{
	float num[] = { 6e-06f, 0.004987f, 0.987f, 394.8f, 3.948e-07f };
	float den[] = { 0.0001f, 0.04f, 9.87f, 3948.0f, 0.0f };
	float A[Nx * Nx], B[Nx], C[Nx], D[1], x[Nx] = { 0 };
	float work[STATESPACE_WORKLEN_FOH(Nx, 1)];
	float t = 0, Ts = 1e-4f, tend = 2, u = 0, y = 0;
	int i, n = (int)(tend / Ts) + 1;
	Statespace ss = { 0 };
	FILE* t_file = fopen("../t.bin", "wb");
	FILE* y_file = fopen("../y.bin", "wb");

	statespace_init(&ss, 4, 1, 1, x, &y, A, B, C, D);
	statespace_tf2ss(&ss, num, den);
	statespace_c2d(&ss, Ts, "tustin", work);

	for (i = 0; i < n; i++)
	{
		t = i * Ts;
		u = sinf(2 * (float)M_PI * 40 * t);
		statespace_iter(&ss, &u, work);
		fwrite(&t, sizeof(float), 1, t_file);
		fwrite(&y, sizeof(float), 1, y_file);
	}

	fclose(t_file);
	fclose(y_file);

	return 0;
}