/*
 * 拡散方程式を解く．
 *
 * Takaaki MINOMO
 */

#include <iostream>
#include <cmath>


//template<typename T = double>
//constexpr T pi = 3.14159265358979323846264338327;

//区間をN等分する．
const int N = 100;
//円周率
const double pi = 3.14159265358979323846264338327;
//区間
const double L = pi;

int main(){

    auto *u = new double[N+1];

    double dx = L / N;
    double x  = 0.;
    
    for(int i=0; i<=N; ++i){
        u[i] = std::sin(x + i*dx);
    }

    for(int i=0; i<=N; ++i){
        std::cout << x + i*dx << " " << u[i] << std::endl;
    }

    delete[] u;
}
