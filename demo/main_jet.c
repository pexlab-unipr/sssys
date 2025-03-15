// This code is based on the Matlab example 
// https://www.mathworks.com/help/control/ug/yaw-damper-design-for-a-747-jet-aircraft.html
#include <stdio.h>
#include "sssys.h"
#include "detego.h"

static void disp(float* A, int rows, int cols)
{
    int i, j;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++)
            printf("%11.7f ", A[i * cols + j]);
        printf("\n");
    }
    printf("\n");
}

#define Nx 4
#define Nu 2
#define Ny 3
#define WORK_LENGTH ((Nx + 2 * Nu) * (4 * (Nx + 2 * Nu) + 1))

int main()
{
    float A[Nx * Nx] = {
        -0.0558f, -0.9968f, 0.0802f, 0.0415f,
        0.5980f, -0.1150f, -0.0318f, 0,
        -3.0500f, 0.3880f ,-0.4650f, 0,
        0, 0.0805f, 1.0000f, 0 };
    float B[Nx * Nu] = {
        0.0073f, 0,
        -0.4750f, 0.0077f,
        0.1530f, 0.1430f,
        0, 0 };
    float C[Ny * Nx] = {
        0, 1, 0, 0,
        0, 0, 0, 1,
        1, 0, 0, 0};
    float D[Ny * Nu] = {
        0, 0,
        0, 0,
        0, 0};
    float x[Nx] = { 0 };
    float work[WORK_LENGTH] = { 0 };
    State_space ss = { 0 };
    int i = 0;
    float u[] = { 10, 100 }, y[] = { 0, 0, 0 };
    float Ts = 0.1f;
    FILE* y1_file = fopen("../yaw_rate.bin", "wb");
    FILE* y2_file = fopen("../bank_angle.bin", "wb");
    FILE* y3_file = fopen("../sideslip_angle.bin", "wb");

    state_space_init(&ss, Nx, Nu, Ny, x, y, A, B, C, D);
    state_space_discretize(&ss, Ts, "foh", work);

    printf("Ad = \n"); disp(A, Nx, Nx);
    printf("Bd = \n"); disp(B, Nx, Nu);
    printf("Cd = \n"); disp(C, Ny, Nx);
    printf("Dd = \n"); disp(D, Ny, Nu);

    while (i++ * Ts <= 10)
    {
        state_space_iter(&ss, u, work);
        fwrite(y + 0, sizeof(float), 1, y1_file);
        fwrite(y + 1, sizeof(float), 1, y2_file);
        fwrite(y + 2, sizeof(float), 1, y3_file);
    }

    fclose(y1_file);
    fclose(y2_file);
    fclose(y3_file);

    return 0;
}