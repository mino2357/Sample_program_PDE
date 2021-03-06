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

// Rittai)個人的に、定数には static もつけます.
// mino)そうします．
// mino)constexpr解禁しました．

constexpr int N = 512;
constexpr double pi = sprout::math::acos(-1.0);
constexpr double L = pi;
constexpr double c = 1.0;
constexpr double D = 0.01;
constexpr double tLimit = 1000;
constexpr double dx = L / N;
constexpr double dt = 0.001;
constexpr double x  = 0.;
//描画のインターバル
constexpr int INTV = 50;

namespace mino2357{
    namespace coeff{

        template <typename T = double>
        constexpr T a(int i, const std::array<T, N+1>& u, const std::array<T, N+1>& g){
            if(c >= 0){
                if(i == 0){
                    return - 2.0 * (u[N] - u[0]) / (dx * dx * dx) + (g[0] + g[N] ) / (dx * dx);
                }
                return - 2.0 * (u[i-1] - u[i]) / (dx * dx * dx) + (g[i] + g[i-1] ) / (dx * dx);
            }else{
                if(i == N){
                    return - 2.0 * (u[N] - u[0]) / (dx * dx * dx) + (g[N] + g[0] ) / (dx * dx);
                }
                return - 2.0 * (u[i+1] - u[i]) / (dx * dx * dx) + (g[i] + g[i+1] ) / (dx * dx);
            }
        }

        template < typename T = double >
        constexpr T b(int i, const std::array< T, N+1>& u, const std::array<T, N+1>& g){
            if(c >= 0){
                if(i == 0){
                    return - 3.0 * (u[0] - u[N]) / (dx * dx) - (2.0 * g[0] + g[N]) / dx;
                }
                return - 3.0 * (u[i] - u[i-1]) / (dx * dx) - (2.0 * g[i] + g[i-1]) / dx;
            }else{
                if(i == N){
                    return - 3.0 * (u[N] - u[0]) / (dx * dx) - (2.0 * g[N] + g[0]) / dx;
                }
                return - 3.0 * (u[i] - u[i+1]) / (dx * dx) - (2.0 * g[i] + g[i+1]) / dx;
            }
        }
    }

    template <typename T = double>
    constexpr T makeNextU(int i, const std::array<T, N+1>& u, const std::array<T, N+1>& g){
        T z = c * dt;
        if(c < 0){z = - 1.0 * z;}
        return coeff::a<>(i, u, g) * z * z * z + coeff::b<>(i, u, g) * z * z + g[i] * z + u[i];
    }
    
    template <typename T = double>
    constexpr T makeNextG(int i, const std::array<T, N+1>& u, const std::array<T, N+1>& g){
        T z = c * dt;
        if(c < 0) z = - 1.0 * z;
        return 3.0 * coeff::a<>(i, u, g) * z * z + 2.0 * coeff::b<>(i, u, g) * z + g[i];
    }
    
    template <typename T = double>
    constexpr T func(T x) {
        //return std::exp(- 40* (x - L/2.) * (x - L/2.));
        if(x > L/2.) return 0.;
        return 1.;
        
    }

    template <typename T = double>
    constexpr void makeGrad(const std::array<T, N+1>& u, std::array<T, N+1>& g){
        for(int i=1; i<N; ++i){
            g[i] = (u[i+1] - u[i-1]) / (2.0 * dx);
        }
        g[0] = (u[1] - u[N])   / (2.0 * dx);
        g[N] = (u[0] - u[N-1]) / (2.0 * dx);
    }
}

int main()
{
	
    //値とその勾配のarray
    std::array<double, N+1> uNew, uOld, uInt;
    std::array<double, N+1> gNew, gOld, gInt;
	
    //初期条件   
    for(int i=0; i<=N; ++i) {
        uOld[i] = mino2357::func(x + i * dx);
    }

    //はじめの勾配を中心差分で求めておく
    mino2357::makeGrad(uOld, gOld);

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
    fprintf(gp, "set yr [-1.0:2.0]\n");
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
        for(int i=0; i<=N; ++i) {
            uInt[i] = mino2357::makeNextU(i, uOld, gOld);
            gInt[i] = mino2357::makeNextG(i, uOld, gOld);
        }
        
        for(int i=1; i<N; ++i) {
            uNew[i] = uInt[i] + dt * D * (uInt[i+1] - 2.0 * uInt[i] + uInt[i-1]) / (dx * dx);
            gNew[i] = gInt[i] + dt * D * (gInt[i+1] - 2.0 * gInt[i] + gInt[i-1]) / (dx * dx);
        }
        uNew[0] = uInt[0] + dt * D * (uInt[1] - 2.0 * uInt[0] + uInt[N]) / (dx * dx);
        uNew[N] = uInt[N] + dt * D * (uInt[0] - 2.0 * uInt[N] + uInt[N-1]) / (dx * dx);
        gNew[0] = gInt[0] + dt * D * (gInt[1] - 2.0 * gInt[0] + gInt[N]) / (dx * dx);
        gNew[N] = gInt[N] + dt * D * (gInt[0] - 2.0 * gInt[N] + gInt[N-1]) / (dx * dx);
        
        //uNewの描画
        if(it%INTV == 0){
            fprintf(gp, "plot '-' w l\n");
            for( int i=0 ; i<=N; ++i ){
                fprintf(gp, "%f %f\n", x + i * dx , uNew[i]);
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }

        // mino) Rittai_3Dさんはこう書くようだ．
        std::copy(begin(uNew), end(uNew), begin(uOld));
        std::copy(begin(gNew), end(gNew), begin(gOld));
        
        //時刻の更新
        t += dt;
    }

    //FILEポインタの解放
    pclose(gp);
}
