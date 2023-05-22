/*
 * One-dimensional space
 * 1st order GKS
 * K = n + 2 = 4 (1-D case)
 * gamma = 1.4; sigma = dt / dx
 * */

#include <iostream>
#include <cmath>
using namespace std;
#define dx 0.0025
#define dt 0.000038
#define gamma 1.4
#define MaxStep 1000
#define PI acos(-1)

void InitialFlowVar(int, int, double [], double [], double []);
void CalF(int, double [], double [], double [], double[], double[], double[]);
void CalFlowVar(int, double, double [],
                double [], double [], double[], double[], double[]);
double CalculateLambda(double, double, double);
void UpdateBoundary(int, int, double [], double [], double []);

int main() {
    // Sod shock tube: flag = 1, dt = 0.00014, dx = 0.01
    // Sjogreen test case: flag = 2, dt = 0.0001, dx = 0.01
    // Woodward-Colella test case: flag = 3, dt = 0.000038, dx = 0.0025
    int N, flag;
    N = 1 / dx + 1;
    flag = 3;

    double x[N], Rho[N], RhoU[N], RhoE[N];
    double FR[N-1], FU[N-1], FE[N-1], p[N], V[N];
    double sigma = dt / dx;

    // Initial the Flow Var
    InitialFlowVar(N, flag, Rho, RhoU, RhoE);

    // Time Increment
    for (int n = 1; n <= MaxStep; n++)
    {
        for (int i = 1; i <= N - 1; i++)
            CalF(i, Rho, RhoU, RhoE, FR, FU, FE);
        for (int i = 1; i <= N - 2; i++)
            CalFlowVar(i, sigma, Rho, RhoU, RhoE, FR, FU, FE);

        // Update the Boundary Condition
        UpdateBoundary(N, flag, Rho, RhoU, RhoE);
    }

    // Calculate the Velocity & Pressure
    // V = RhoU / Rho
    // RhoE = p / (gamma - 1) + 0.5 * Rho * V ^ 2
    for (int i = 0; i <= N - 1; i++)
    {
        V[i] = RhoU[i] / Rho[i];
        p[i] = (gamma - 1) * (RhoE[i] - 0.5 * Rho[i] * pow(RhoU[i] / Rho[i], 2));
        cout << p[i] << endl;
    }
    return 0;
}

void InitialFlowVar(int N, int flag, double Rho[], double RhoU[], double RhoE[])
{
    // Initial Condition
    if (flag == 1)
    {
        for (int i = 0; i <= N - 1; i++)
        {
            if ((double)i * dx <= 0.5)
            {
                Rho[i] = 1; RhoU[i] = 0; RhoE[i] = 2.5;
            }
            else
            {
                Rho[i] = 0.125; RhoU[i] = 0; RhoE[i] = 0.25;
            }
        }
    }
    else if (flag == 2)
    {
        for (int i = 0; i <= N - 1; i++)
        {
            if ((double)i * dx <= 0.5)
            {
                Rho[i] = 1; RhoU[i] = -2; RhoE[i] = 3;
            }
            else
            {
                Rho[i] = 1; RhoU[i] = 2; RhoE[i] = 3;
            }
        }
    }
    else
    {
        for (int i = 0; i <= N - 1; i++)
        {
            if ((double)i * dx <= 0.1)
            {
                Rho[i] = 1; RhoU[i] = 0; RhoE[i] = 2500;
            }
            else if ((double)i * dx > 0.9)
            {
                Rho[i] = 1; RhoU[i] = 0; RhoE[i] = 250;
            }
            else
            {
                Rho[i] = 1; RhoU[i] = 0; RhoE[i] = 0.025;
            }
        }
    }
}

void CalF(int i, double Rho[], double RhoU[], double RhoE[],
                double FR[], double FU[], double FE[])
{
    double lambda0, lambda1, U0, U1, K;
    K = 4;
    U0 = RhoU[i-1] / Rho[i-1];
    U1 = RhoU[i] / Rho[i];
    lambda0 = CalculateLambda(Rho[i-1], RhoU[i-1], RhoE[i-1]);
    lambda1 = CalculateLambda(Rho[i], RhoU[i], RhoE[i]);

    FR[i-1] = Rho[i-1] * (0.5 * U0 * erfc(-sqrt(lambda0) * U0)
            + 0.5 * exp(-lambda0 * U0 * U0) / sqrt(PI * lambda0))
            + Rho[i] * (0.5 * U1 * erfc(sqrt(lambda1) * U1)
            - 0.5 * exp(-lambda1 * U1 * U1) / sqrt(PI * lambda1));

    FU[i-1] = Rho[i-1] * ((U0 * U0 / 2 + 1 / (4 * lambda0)) * erfc(-sqrt(lambda0) * U0)
            + 0.5 * U0 * exp(-lambda0 * U0 * U0) / sqrt(PI * lambda0))
            + Rho[i] * ((U1 * U1 / 2 + 1 / (4 * lambda1)) * erfc(sqrt(lambda1) * U1)
            - 0.5 * U1 * exp(-lambda1 * U1 * U1) / sqrt(PI * lambda1));

    FE[i-1] = Rho[i-1] * ((U0 * U0 * U0 / 4 + (K + 3) / (8 * lambda0) * U0)
            * erfc(-sqrt(lambda0) * U0) + (U0 * U0 / 4 + (K + 2) / (8 * lambda0))
            * exp(-lambda0 * U0 * U0) / sqrt(PI * lambda0))
            + Rho[i] * ((U1 * U1 * U1 / 4 + (K + 3) / (8 * lambda1) * U1)
            * erfc(sqrt(lambda1) * U1) - (U1 * U1 / 4 + (K + 2) / (8 * lambda1))
            * exp(-lambda1 * U1 * U1) / sqrt(PI * lambda1));

}

void CalFlowVar(int i, double sigma, double Rho[], double RhoU[]
                      , double RhoE[], double FR[], double FU[], double FE[])
{
    Rho[i] += sigma * (FR[i-1] - FR[i]);
    RhoU[i] += sigma * (FU[i-1] - FU[i]);
    RhoE[i] += sigma * (FE[i-1] - FE[i]);
}

double CalculateLambda(double Rho, double RhoU, double RhoE)
{
    double lambda, K;
    K = 4;
    lambda = ((K + 1) / 4) * (Rho / (RhoE - 0.5 * Rho * pow(RhoU/Rho, 2)));
    return lambda;
}

void UpdateBoundary(int N, int flag, double Rho[], double RhoU[], double RhoE[])
{
    Rho[0] = Rho[1]; Rho[N-1] = Rho[N-2];
    RhoE[0] = RhoE[1]; RhoE[N-1] = RhoE[N-2];
    if (flag == 1 || flag == 2)
    {
        RhoU[0] = RhoU[1]; RhoU[N-1] = RhoU[N-2];
    }
    else
    {
        RhoU[0] = -RhoU[1]; RhoU[N-1] = -RhoU[N-2];
    }
}
