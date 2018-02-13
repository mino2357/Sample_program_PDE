#pragma once

constexpr int N = 512;
constexpr double alpha = 1.0;
constexpr double beta  = 0.022;
constexpr double Re    = 1000.0;
constexpr double Lx =  2.0;
constexpr double Ly =  1.0;
constexpr double xstart = 0.0;
constexpr double ystart = 0.0;
constexpr double dx = Lx / N;
constexpr double dy = Ly / N;
constexpr double dt = dx*dx*dx*dx;
constexpr double pi = 3.14159265358979323846264338327950288;
constexpr double tLimit = 2000000.0 * pi;
double Tol = 10e-8;
constexpr int INTV = 1000000;
constexpr int plus = 1;
constexpr double vecLen = static_cast<double>(N);
constexpr int istart = 0;
constexpr int jstart = 0;
constexpr double nu = 0.0;
int Lim = 100;
