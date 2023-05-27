/*
 * One dimensional space
 * 2nd order GKS
 * K = n + 2 = 4 (1-D case)
 * gamma = 1.4
 * */

#include <iostream>
#include <cmath>
using namespace std;

#define dx 0.01
#define dt 0.0001
#define gamma 1.4
#define MaxStep 1000
#define PI acos(-1)

void InitialFlowVar(int, int, double [], double [], double []);
void Vanleer(int, double [], double [], double [], double [], double []);
void CalFlux(int, const double [], double [], double [], double [],
             double [], double [], double [], double []);
void UpdateBoundary(int, int, double [], double [], double []);
double Sign(double);

int main() {
    // Sod shock tube: flag = 1, dt = 0.00014, dx = 0.01
    // Sjogreen test case: flag = 2, dt = 0.0001, dx = 0.01
    // Woodward-Colella test case: flag = 3, dt = 0.000038, dx = 0.0025
    int N, flag;
    N = 1 / dx + 1;
    flag = 2;

    double Rho[N+4], RhoU[N+4], RhoE[N+4], p[N+4], U[N+4];
    double FR[N+1], FU[N+1], FE[N+1], a[6], Par[6];

    // Initial the Flow Var
    InitialFlowVar(N, flag, Rho, RhoU, RhoE);

    for (int n = 1; n <= MaxStep; n++)
    {
        for (int i = 0; i <= N; i++)
        {
            Vanleer(i+1, Rho, RhoE, RhoU, a, Par);
            CalFlux(i+1, a, Rho, RhoU, RhoE, FR, FU, FE, Par);
        }

        for (int i = 2; i <= N+1; i++)
        {
            Rho[i] += (FR[i-2] - FR[i-1]) / dx;
            RhoU[i] += (FU[i-2] - FU[i-1]) / dx;
            RhoE[i] += (FE[i-2] - FE[i-1]) / dx;
        }
        UpdateBoundary(N, flag, Rho, RhoU, RhoE);
    }
    for (int i = 0; i <= N + 3; i++)
    {
        U[i] = RhoU[i] / Rho[i];
        p[i] = (gamma - 1) * (RhoE[i] - 0.5 * Rho[i] * pow(U[i], 2));
        cout << p[i] << endl;
    }
    return 0;
}

void CalFlux(int i, const double a[], double Rho[], double RhoU[], double RhoE[],
             double FR[], double FU[], double FE[], double Par[])
{
    double Fl, Fr, F1, F2, lambda, U;
    double Rhol, Rhor, RhoUl, RhoUr, RhoEl, RhoEr;
    // Calculate Flux RhoE
    Rhor = Rho[i+1] - 0.5*dx*Par[3];
    RhoUr = RhoU[i+1] - 0.5*dx*Par[4];
    RhoEr = RhoE[i+1] - 0.5*dx*Par[5];
    U = RhoUr / Rhor;
    if (RhoEr / Rhor - 0.5 * U * U < 0)
    {
        Rhor = Rho[i+1];
        RhoUr = RhoU[i+1];
        RhoEr = RhoE[i+1];
        U = RhoUr / Rhor;
    }
    lambda = 1.25 * (Rhor / (RhoEr - 0.5 * Rhor * U * U));
    F1 = exp(-U*U*lambda)*dt*sqrt(lambda)*(dt*(32*a[4]+97*a[5]*U)+
         4*(-12+dt*U*(9*a[3]+U*(13*a[4]+11*a[5]*U)))*lambda+4*U*U*
         (-4+dt*U*(2*a[3]+U*(2*a[4]+a[5]*U)))*lambda*lambda)/sqrt(PI);
    F2 = -0.5*dt*(4*lambda*(dt*(7*a[3]+27*a[4]*U)+4*U*(-7+dt*U*(5*a[3]+7*a[4]*U))*lambda
         +4*U*U*U*(-2+dt*U*(a[3]+a[4]*U))*lambda*lambda)+a[5]*dt*(63+234*U*U*lambda+
         92*pow(U,4)*lambda*lambda+8*pow(U,6)*pow(lambda,3)))*erfc(U*sqrt(lambda));
    Fr = Rhor * (F1 + F2) / (64 * pow(lambda,3));

    Rhol = Rho[i] + 0.5*dx*Par[0];
    RhoUl = RhoU[i] + 0.5*dx*Par[1];
    RhoEl = RhoE[i] + 0.5*dx*Par[2];
    U = RhoUl / Rhol;
    if (RhoEl / Rhol - 0.5 * U * U < 0)
    {
        Rhol = Rho[i];
        RhoUl = RhoU[i];
        RhoEl = RhoE[i];
        U = RhoUl / Rhol;
    }
    lambda = 1.25 * (Rhol / (RhoEl - 0.5 * Rhol * U * U));
    F1 = -exp(-U*U*lambda)*dt*sqrt(lambda)*(dt*(32*a[1]+97*a[2]*U)+
         4*(-12+dt*U*(9*a[0]+U*(13*a[1]+11*a[2]*U)))*lambda+4*U*U*
         (-4+dt*U*(2*a[0]+U*(2*a[1]+a[2]*U)))*lambda*lambda)/sqrt(PI);
    F2 = -0.5*dt*(4*lambda*(dt*(7*a[0]+27*a[1]*U)+4*U*(-7+dt*U*(5*a[0]+7*a[1]*U))*lambda
         +4*U*U*U*(-2+dt*U*(a[0]+a[1]*U))*lambda*lambda)+a[2]*dt*(63+234*U*U*lambda+
         92*pow(U,4)*lambda*lambda+8*pow(U,6)*pow(lambda,3)))*erfc(-U*sqrt(lambda));
    Fl = Rhol * (F1 + F2) / (64 * pow(lambda,3));

    FE[i-1] = Fr + Fl;

    // Calculate Flux RhoU
    Rhor = Rho[i+1] - 0.5*dx*Par[3];
    RhoUr = RhoU[i+1] - 0.5*dx*Par[4];
    RhoEr = RhoE[i+1] - 0.5*dx*Par[5];
    U = RhoUr / Rhor;
    if (RhoEr / Rhor - 0.5 * U * U < 0)
    {
        Rhor = Rho[i+1];
        RhoUr = RhoU[i+1];
        RhoEr = RhoE[i+1];
        U = RhoUr / Rhor;
    }
    lambda = 1.25 * (Rhor / (RhoEr - 0.5 * Rhor * U * U));
    F1 = exp(-U*U*lambda)*dt*(2*lambda*(dt*(2*a[3]+5*a[4]*U)+2*U*
    (-2+dt*U*(a[3]+a[4]*U))*lambda)+a[5]*dt*(8+U*U*lambda*
    (13+2*U*U*lambda)))/sqrt(PI);
    F2 = -0.5*dt*sqrt(lambda)*(3*dt*(2*a[4]+9*a[5]*U)+4*(-2+dt*U*(3*a[3]+
    U*(6*a[4]+7*a[5]*U)))*lambda+4*U*U*(-4+dt*U*(2*a[3]+U*(2*a[4]+a[5]*U)))
    *lambda*lambda)*erfc(U*sqrt(lambda));
    Fr = Rhor * (F1 + F2) / (16 * pow(lambda,2.5));

    Rhol = Rho[i] + 0.5*dx*Par[0];
    RhoUl = RhoU[i] + 0.5*dx*Par[1];
    RhoEl = RhoE[i] + 0.5*dx*Par[2];
    U = RhoUl / Rhol;
    if (RhoEl / Rhol - 0.5 * U * U < 0)
    {
        Rhol = Rho[i];
        RhoUl = RhoU[i];
        RhoEl = RhoE[i];
        U = RhoUl / Rhol;
    }
    lambda = 1.25 * (Rhol / (RhoEl - 0.5 * Rhol * U * U));
    F1 = -exp(-U*U*lambda)*dt*(2*lambda*(dt*(2*a[0]+5*a[1]*U)+2*U*
         (-2+dt*U*(a[0]+a[1]*U))*lambda)+a[2]*dt*(8+U*U*lambda*
         (13+2*U*U*lambda)))/sqrt(PI);
    F2 = -0.5*dt*sqrt(lambda)*(3*dt*(2*a[1]+9*a[2]*U)+4*(-2+dt*U*(3*a[0]+
         U*(6*a[1]+7*a[2]*U)))*lambda+4*U*U*(-4+dt*U*(2*a[0]+U*(2*a[1]+a[2]*U)))
         *lambda*lambda)*erfc(-U*sqrt(lambda));
    Fl = Rhol * (F1 + F2) / (16 * pow(lambda,2.5));

    FU[i-1] = Fr + Fl;

    // Calculate Flux Rho
    Rhor = Rho[i+1] - 0.5*dx*Par[3];
    RhoUr = RhoU[i+1] - 0.5*dx*Par[4];
    RhoEr = RhoE[i+1] - 0.5*dx*Par[5];
    U = RhoUr / Rhor;
    if (RhoEr / Rhor - 0.5 * U * U < 0)
    {
        Rhor = Rho[i+1];
        RhoUr = RhoU[i+1];
        RhoEr = RhoE[i+1];
        U = RhoUr / Rhor;
    }
    lambda = 1.25 * (Rhor / (RhoEr - 0.5 * Rhor * U * U));
    F1 = exp(-U*U*lambda)*dt*sqrt(lambda)*(4*(-2+a[3]*dt*U)*lambda+a[5]*dt*U*(9+2*U*U*lambda)
                                           +4*a[4]*(dt+dt*U*U*lambda))/sqrt(PI);
    F2 = -0.5*dt*(4*lambda*(dt*(a[3]+3*a[4]*U)+2*U*(-2+dt*U*(a[3]+a[4]*U))*lambda)+
                  a[5]*dt*(7+4*U*U*lambda*(5+U*U*lambda)))*erfc(U*sqrt(lambda));
    Fr = Rhor * (F1 + F2) / (16 * pow(lambda, 2));

    Rhol = Rho[i] + 0.5*dx*Par[0];
    RhoUl = RhoU[i] + 0.5*dx*Par[1];
    RhoEl = RhoE[i] + 0.5*dx*Par[2];
    U = RhoUl / Rhol;
    if (RhoEl / Rhol - 0.5 * U * U < 0)
    {
        Rhol = Rho[i];
        RhoUl = RhoU[i];
        RhoEl = RhoE[i];
        U = RhoUl / Rhol;
    }
    lambda = 1.25 * (Rhor / (RhoEl - 0.5 * Rhol * U * U));
    F1 = -exp(-U*U*lambda)*dt*sqrt(lambda)*(4*(-2+a[0]*dt*U)*lambda+a[2]*dt*U*(9+2*U*U*lambda)
                                            +4*a[1]*(dt+dt*U*U*lambda))/sqrt(PI);
    F2 = -0.5*dt*(4*lambda*(dt*(a[0]+3*a[1]*U)+2*U*(-2+dt*U*(a[0]+a[1]*U))*lambda)+
                  a[2]*dt*(7+4*U*U*lambda*(5+U*U*lambda)))*erfc(-U*sqrt(lambda));
    Fl = Rhol * (F1 + F2)/ (16 * pow(lambda, 2));

    FR[i-1] = Fr + Fl;
}

void InitialFlowVar(int N, int flag, double Rho[], double RhoU[], double RhoE[])
{
    // Initial Condition
    if (flag == 1)
    {
        for (int i = 0; i <= N + 3; i++)
        {
            if ((double)i * dx <= 0.5 + 2 * dx)
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
        for (int i = 0; i <= N + 3; i++)
        {
            if ((double)i * dx <= 0.5 + 2 * dx)
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
        for (int i = 0; i <= N + 3; i++)
        {
            if ((double)i * dx <= 0.1 + 2 * dx)
            {
                Rho[i] = 1; RhoU[i] = 0; RhoE[i] = 2500;
            }
            else if ((double)i * dx > 0.9 + 2 * dx)
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

void UpdateBoundary(int N, int flag, double Rho[], double RhoU[], double RhoE[])
{
    Rho[0] = Rho[3]; Rho[1] = Rho[2]; Rho[N+3] = Rho[N]; Rho[N+2] = Rho[N+1];
    RhoE[0] = RhoE[3]; RhoE[1] = RhoE[2]; RhoE[N+3] = RhoE[N]; RhoE[N+2] = RhoE[N+1];
    RhoU[0] = RhoU[3]; RhoU[1] = RhoU[2]; RhoU[N+3] = RhoU[N]; RhoU[N+2] = RhoU[N+1];
    if (flag == 3) {
        RhoU[0] = -RhoU[3];
        RhoU[1] = -RhoU[2];
        RhoU[N+3] = -RhoU[N];
        RhoU[N+2] = -RhoU[N+1];
    }
}

void Vanleer(int i, double Rho[], double RhoE[], double RhoU[], double a[], double Par[])
{
    // a[]: Include 6 parameters(a1l,a2l,a3l,a1r,a2r,a3r)
    // K = 4;
    // Left side
    double s, r, Rhox, RhoEx, RhoUx, lambda, lambdax, Ex, Ux, U;
    s = (Rho[i+1]-Rho[i]) / dx;
    r = (Rho[i]-Rho[i-1]) / dx;
    Rhox = (Sign(s)+Sign(r)) * abs(s * r) / (abs(s) + abs(r));
    if (s == 0 && r == 0)
        Rhox = 0;

    s = (RhoU[i+1]-RhoU[i]) / dx;
    r = (RhoU[i]-RhoU[i-1]) / dx;
    RhoUx = (Sign(s)+Sign(r)) * abs(s * r) / (abs(s) + abs(r));
    if (s == 0 && r == 0)
        RhoUx = 0;


    s = (RhoE[i+1]-RhoE[i]) / dx;
    r = (RhoE[i]-RhoE[i-1]) / dx;
    RhoEx = (Sign(s)+Sign(r)) * abs(s * r) / (abs(s) + abs(r));
    if (s == 0 && r == 0)
        RhoEx = 0;

    U = RhoU[i] / Rho[i];
    lambda = 1.25 * (Rho[i] / (RhoE[i] - 0.5 * Rho[i] * U * U));
    Ex = RhoEx / Rho[i] - Rhox*RhoE[i]/(Rho[i]*Rho[i]);
    Ux = -Rhox*U/Rho[i]+RhoUx/Rho[i];
    lambdax = 1.25 / pow(RhoE[i]/Rho[i]-0.5*U*U,2)*(-Ex+U*Ux);
    a[0] = Rhox / Rho[i] - 2 * lambda * Ux * U +
           (5/(2*lambda)-U*U) * lambdax;
    a[1] = 2 * lambda * Ux + 2 * U * lambdax;
    a[2] = -2 * lambdax;

    Par[0] = Rhox;
    Par[1] = RhoUx;
    Par[2] = RhoEx;

    // Right side
    s = (Rho[i+2]-Rho[i+1]) / dx;
    r = (Rho[i+1]-Rho[i]) / dx;
    Rhox = (Sign(s)+Sign(r)) * abs(s * r) / (abs(s) + abs(r));
    if (s == 0 && r == 0)
        Rhox = 0;

    s = (RhoU[i+2]-RhoU[i+1]) / dx;
    r = (RhoU[i+1]-RhoU[i]) / dx;
    RhoUx = (Sign(s)+Sign(r)) * abs(s * r) / (abs(s) + abs(r));
    if (s == 0 && r == 0)
        RhoUx = 0;

    s = (RhoE[i+2]-RhoE[i+1]) / dx;
    r = (RhoE[i+1]-RhoE[i]) / dx;
    RhoEx = (Sign(s)+Sign(r)) * abs(s * r) / (abs(s) + abs(r));
    if (s == 0 && r == 0)
        RhoEx = 0;


    U = RhoU[i+1] / Rho[i+1];
    lambda = 1.25 * (Rho[i+1] / (RhoE[i+1] - 0.5 * Rho[i+1] * U * U));
    Ex = RhoEx / Rho[i+1] - Rhox*RhoE[i+1]/(Rho[i+1]*Rho[i+1]);
    Ux = -Rhox*RhoU[i+1]/(Rho[i+1]*Rho[i+1])+RhoUx/Rho[i+1];
    lambdax = 1.25 / pow(RhoE[i+1]/Rho[i+1]-0.5*pow(RhoU[i+1]/Rho[i+1],2),2)
              *(-Ex+RhoU[i+1]*Ux/Rho[i+1]);
    a[3] = Rhox / Rho[i+1] - 2 * lambda * Ux * RhoU[i+1] / Rho[i+1] +
           (5/(2*lambda)-pow(RhoU[i+1]/Rho[i+1],2)) * lambdax;
    a[4] = 2 * lambda * Ux + 2 * RhoU[i+1]/Rho[i+1] * lambdax;
    a[5] = -2 * lambdax;
    Par[3] = Rhox;
    Par[4] = RhoUx;
    Par[5] = RhoEx;
}


double Sign(double a)
{
    if (a > 0)
        return 1;
    else if(a < 0)
        return -1;
    else
        return 0;
}


