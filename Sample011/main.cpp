/*
 * cavity flow test.
 *
 * kazakami method. SMAC method.
 * 
 * Takaaki MINOMO.
 */

#include <iostream>
#include <cmath>
#include <array>
#include <limits>
#include <iomanip>
#include "parameter.hpp"
#include "gnuplot.hpp"
#include "extendedArray.hpp"

using Array = sksat::array_wrapper<N + 1, double>;

namespace mino2357{
    
    template <typename T = double>
    void setInitFuncF(Array& f){
        T x;
        for(int i=0; i<=N; ++i){
           	x = xstart + i * dx;
			f[i] = cos(pi * x);
        }
    }

	void succU(Array& u1, Array& u2){
		for(int i=0; i<=N; ++i){
			if(u1[i] > 0){
				u2[i] = u1[i] - alpha * dt * u1[i] * (u1[i] - u1[i-1]) / dx 
						- beta * beta * dt * ( - u1[i-2] + 2.0 * u1[i-1] - 2.0 * u1[i+1] + u1[i+2]) / (2.0 * dx * dx * dx);
			}else{
				u2[i] = u1[i] - alpha * dt * u1[i] * (u1[i+1] - u1[i]) / dx 
						- beta * beta * dt * ( - u1[i-2] + 2.0 * u1[i-1] - 2.0 * u1[i+1] + u1[i+2]) / (2.0 * dx * dx * dx);
			}
		}
	}
}


int main(){
    std::cout << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 - 8);

    auto u1 = Array{};
    auto u2 = Array{};
    
	auto gp = mino2357::gnuplot();    
	mino2357::setInitFuncF(u1);

	gp.print(u1);	

    std::cout << "Enter!" << std::endl;
    getchar();

	double t = 0.0;

    for(int it=0; t<=tLimit; ++it){

		mino2357::succU(u1, u2);
		t += dt;

        if(it%INTV == 0){
			std::cout << t << std::endl;
			gp.print(u2);
            //gp.printP(p2);
            //gp.speed(f3, g3);
            //gp.vector(f3, g3);
            //gp.vectorWithP(f3, g3, p2);
            //gp.vectorWithSpeed(f3, g3);
            //gp.vectorWithSpeedLog10(f3, g3);
            //gp.vectorWithLog10P(f3, g3, p2);
            //gp.multiplot(f3, g3, p2, t);
            //gp.multiplotMakePNG(f3, g3, p2, t);
			//gp.printRho(rho2);	
        }
		for(int i=0; i<=N; ++i){
			u1[i] = u2[i];
		}
    }
}
