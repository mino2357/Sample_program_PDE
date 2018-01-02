/*
 * 移流拡散方程式を解く．(周期境界条件)
 *
 * 移流部分：CIP method ( F(X) = a3 X^3 + a2 X^2 + a1 X + a0, X = x - x_i
 *
 * uOldのデータからuNewを構成していく．
 *
 * Takaaki MINOMO
 * 私のコメントには先頭に mino) がつきます
 *
 * 改善（プログラムの監修）:Rittai_3D（仮名）
 * わたしのコメントには先頭に Rittai) がつきます
 */

#include <iostream>
#include <array>
#include <cmath>
#include <sprout/cmath.hpp>

constexpr int N = 512;
constexpr double pi = sprout::math::acos(-1.0);
constexpr double L = pi;
constexpr double c = 1.0;
constexpr double D = 0.01;
constexpr double tLimit = 1000;
constexpr double dx = L / N;
constexpr double dt = 0.001;
constexpr double x  = 0.;
constexpr int INTV = 50;

namespace mino2357{
    template <typename T = double>
    constexpr T func(T x) {
        return std::exp(- 40* (x - L/2.) * (x - L/2.));
        //if(x > L/2.) return 0.;
        //return 1.0;    
    }
}

int main()
{
    std::array<double, N+1> uNew, uOld;
	
    //初期条件   
    for(int i=0; i<=N; ++i) {
        uOld[i] = mino2357::func(x + i * dx);
    }

    //初期の時刻．
    double t = 0;

    //std::cout << (dt * c) / dx << std::endl;
    //安定性条件のチェック CFL条件
    if(1.0 < std::abs((dt * c) / dx) || D * dt / (dx * dx) > 0.5){
        std::cout << std::abs((dt * c) / dx) << " " << D * dt / (dx * dx) << std::endl;
        std::cout << "安定性条件を満たしていません．" << std::endl;
        return 0;
    }

    /**********************************************************************/
    /*                 可視化の設定(gnuplot)                              */
    /**********************************************************************/
    std::FILE *gp = popen( "gnuplot -persist", "w" );
    fprintf(gp, "set xr [0:%f]\n", L);
    fprintf(gp, "set yr [-0.5:1.5]\n");
    //fprintf(gp, "set size square\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "unset key\n");
    //fprintf(gp, "set term dumb\n");
    
    //初期条件描画
    fprintf(gp, "plot '-' w l\n");
    for(int i=0; i<=N; ++i){
        fprintf(gp, "%f %f\n", x + i * dx , uOld[i]);
    }
    fprintf(gp, "e\n");
    fflush(gp);

    getchar();

    //タイムループ
    for(int it = 0; t<tLimit; ++it) {   

        //uをdt時間だけ進める
        if(c >= 0){
            for(int i=1; i<N; ++i){
                uNew[i] = uOld[i] - c * dt / dx * (uOld[i] - uOld[i-1]) + D * dt / (dx * dx) * (uOld[i-1] - 2.0 * uOld[i] + uOld[i+1]);
            }
            uNew[0] = uOld[0] - c * dt / dx * (uOld[0] - uOld[N]) + D * dt / (dx * dx) * (uOld[N] - 2.0 * uOld[0] + uOld[1]);
            uNew[N] = uOld[N] - c * dt / dx * (uOld[N] - uOld[N-1]) + D * dt / (dx * dx) * (uOld[N-1] - 2.0 * uOld[N] + uOld[0]);
        }else if(c < 0){
            for(int i=1; i<N; ++i){
                uNew[i] = uOld[i] - c * dt / dx * (uOld[i+1] - uOld[i]) + D * dt / (dx * dx) * (uOld[i-1] - 2.0 * uOld[i] + uOld[i+1]);
            }
            uNew[0] = uOld[0] - c * dt / dx * (uOld[1] - uOld[0]) + D * dt / (dx * dx) * (uOld[N] - 2.0 * uOld[0] + uOld[1]);
            uNew[N] = uOld[N] - c * dt / dx * (uOld[0] - uOld[N]) + D * dt / (dx * dx) * (uOld[N-1] - 2.0 * uOld[N] + uOld[0]);
        }


        //uNewの描画
        if(it%INTV == 0){
            fprintf(gp, "plot '-' w l\n");
            for( int i=0 ; i<=N; ++i ){
                fprintf(gp, "%f %f\n", x + i * dx , uNew[i]);
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }

        //更新
        std::copy(begin(uNew), end(uNew), begin(uOld));
        
        //時刻の更新
        t += dt;
    }

    //FILEポインタの解放
    pclose(gp);
}
