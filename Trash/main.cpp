#include <iostream>
#include <limits>
#include <iomanip>
#include <cmath>
#include <Eigen/Sparse>
//#include <Eigen/Dense>

constexpr int N = 1024;
constexpr double L = 1.;
constexpr double dx = L / (N + 1);
constexpr double xStart = 0 + dx;
constexpr double xEnd   = L - dx;
constexpr double boundaryA = 0;
constexpr double boundaryB = 1;

using VectorNd = Eigen::Matrix<double, N, 1>;

template <typename T = double>
constexpr T func(const T& x) noexcept{
    return -x;
}

template <typename T = double>
constexpr T exactSol(const T& x) noexcept{
    return - 1. / 6. * x * (x * x - 7.);
}

int main(){
	std::cout << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	//std::cout << std::fixed << std::setprecision(6);
    
    Eigen::SparseMatrix<double> A(N, N);
    VectorNd b;
    VectorNd x;

    for(int i=0; i<N; ++i){
        A.insert(i, i) = - 2.;
        b(i) = dx * dx * func<>(xStart + i * dx);
    }

    b(0)   = dx * dx * func<>(0) - boundaryA;
    b(N-1) = dx * dx * func<>(L) - boundaryB;

    for(int i=0; i<N-1; ++i){
        A.insert(i, i+1) = 1.;
        A.insert(i+1, i) = 1.;
    }

    A.makeCompressed();

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

    solver.analyzePattern(A);
    
    solver.factorize(A);

    x = solver.solve(b);
    
    
    std::cout << xStart - dx << " " << boundaryA - exactSol<>(xStart - dx) << std::endl;
    for(int i{}; i<N; ++i){
        std::cout << xStart + i * dx << " " << x(i) - exactSol<>(xStart + i * dx) << std::endl;
    }
    std::cout << xEnd + dx << " " << boundaryB - exactSol<>(xEnd + dx) << std::endl;
    
/*
    std::cout << xStart - dx << " " << boundaryA << std::endl;
    for(int i{}; i<N; ++i){
        std::cout << xStart + i * dx << " " << x(i) << std::endl;
    }
    std::cout << xEnd + dx << " " << boundaryB << std::endl;
*/
}
