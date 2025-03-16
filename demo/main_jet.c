/* This code is based on the Matlab example
https://www.mathworks.com/help/control/ug/yaw-damper-design-for-a-747-jet-aircraft.html
The following Matlab script can be use to validate the resaults.
clear; close all; clc;

%% Continuous time system
A = [-0.0558   -0.9968    0.0802    0.0415
    0.5980   -0.1150   -0.0318         0
    -3.0500    0.3880   -0.4650         0
    0    0.0805    1.0000         0];
B = [ 0.0073         0
    -0.4750    0.0077
    0.1530    0.1430
    0         0];
C = [0     1     0     0
    0     0     0     1
    1   0   0   0];
D = [0     0
    0     0
    0   0];

t_ref = 0:1e-5:10;
u = [10; 100];
sys_ref = ss(A,B,C,D);
y_ref = lsim(sys_ref, repmat(u, size(t_ref)), t_ref);

%% Discrete time system
Ts = 0.1;
t = 0:Ts:10;

sysc = ss(A,B,C,D);
sysd = c2d(sysc,Ts,'zoh');

Ad = (sysd.A);
Bd = (sysd.B);
Cd = (sysd.C);
Dd = (sysd.D);

x = zeros(size(A,1),1);
y = zeros(size(C,1), length(t));
for k = 1:length(t)
    y(:,k) = Cd*x + Dd*u;
    x = Ad*x + Bd*u;
end

%% Data import
p = "./";
fid = fopen(p + "yaw_rate.bin");
y1 = fread(fid,'single');
fclose(fid);

fid = fopen(p + "bank_angle.bin");
y2 = fread(fid,'single');
fclose(fid);

fid = fopen(p + "sideslip_angle.bin");
y3 = fread(fid,'single');
fclose(fid);

%% Figures
figure
plot(t_ref, y_ref','k', t, y', 'r',...
    t, y1(:), 'b-.', t, y2(:), 'b-.', t, y3(:), 'b-.', 'LineWidth', 1);
grid on

figure
subplot(3,1,1); plot(y1(:) - y(1,:)'); grid on;
subplot(3,1,2); plot(y2(:) - y(2,:)'); grid on;
subplot(3,1,3); plot(y3(:) - y(3,:)'); grid on;
*/

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
    float work[STATESPACE_WORKLEN_FOH(Nx, Nu)] = {0};
    Statespace ss = { 0 };
    int i = 0;
    float u[] = { 10, 100 }, y[] = { 0, 0, 0 };
    float Ts = 0.1f;
    FILE* y1_file = fopen("../yaw_rate.bin", "wb");
    FILE* y2_file = fopen("../bank_angle.bin", "wb");
    FILE* y3_file = fopen("../sideslip_angle.bin", "wb");

    statespace_init(&ss, Nx, Nu, Ny, x, y, A, B, C, D);
    statespace_c2d(&ss, Ts, "zoh", work);

    printf("Ad = \n"); disp(A, Nx, Nx);
    printf("Bd = \n"); disp(B, Nx, Nu);
    printf("Cd = \n"); disp(C, Ny, Nx);
    printf("Dd = \n"); disp(D, Ny, Nu);

    while (i++ * Ts <= 10)
    {
        statespace_iter(&ss, u, work);
        fwrite(y + 0, sizeof(float), 1, y1_file);
        fwrite(y + 1, sizeof(float), 1, y2_file);
        fwrite(y + 2, sizeof(float), 1, y3_file);
    }

    fclose(y1_file);
    fclose(y2_file);
    fclose(y3_file);

    return 0;
}