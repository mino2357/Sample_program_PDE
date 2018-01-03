/*
 * 非粘性バーガース方程式を解く．(周期境界条件)
 *
 * 移流部分は風上差分法．
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

//パラメータ
constexpr int N = 512;
constexpr double pi = sprout::math::acos(-1.0);
constexpr double L = pi;
constexpr double c = 1.0;
constexpr double tLimit = 1000;
constexpr double dx = L / N;
constexpr double dt = 0.0001;
constexpr double x  = 0.;
constexpr int INTV = 100;

namespace mino2357{
    //初期条件
    template <typename T = double>
    constexpr T func(T x)noexcept {
        return 0.5 * std::sin(2.0 * x) + 1.0;
        //return std::exp(- 40* (x - L/2.) * (x - L/2.));
        //if(x > L/2.) return 0.;
        //return 1.;
        
    }
}

int main()
{
	
    //値とその勾配のarray
    std::array<double, N+1> uNew, uOld;
	
    //初期条件   
    for(int i=0; i<=N; ++i) {
        uOld[i] = mino2357::func(x + i * dx);
    }

    //初期の時刻．
    double t = 0;

    //安定性条件のチェック CFL条件
    if(1.0 < std::abs((dt * c) / dx)){
        std::cout << "安定性条件を満たしていません．" << std::endl;
        return 0;
    }

    /**********************************************************************/
    /*                 可視化の設定(gnuplot)                              */
    /**********************************************************************/
    std::FILE *gp = popen( "gnuplot -persist", "w" );
    fprintf(gp, "set xr [0:%f]\n", L);
    fprintf(gp, "set yr [-1:2]\n");
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

    std::cout << "Enterキーを押してください．" << std::endl;//実はなんでも良い（何でもとは言っていない）
    getchar();

    //タイムループ
    for(int it = 0; t<tLimit; ++it) {   

        //uをdtだけ進める
        for(int i=1; i<N; ++i){
            if(c*uOld[i] >= 0){
                uNew[i] = uOld[i] - c * uOld[i] * dt / dx * (uOld[i] - uOld[i-1]);
            }else{
                uNew[i] = uOld[i] - c * uOld[i] * dt / dx * (uOld[i+1] - uOld[i]);
            }    
        }

        if(c * uOld[0] >= 0){
            uNew[0] = uOld[0] - c * uOld[0] * dt / dx * (uOld[0] - uOld[N]);
        }else if(c * uOld[0] < 0){
            uNew[0] = uOld[0] - c * uOld[0] * dt / dx * (uOld[1] - uOld[0]);
        }
            
        if(c * uOld[N] >= 0){
            uNew[N] = uOld[N] - c * uOld[N] * dt / dx * (uOld[N] - uOld[N-1]);
        }else if(c * uOld[0] < 0){
            uNew[N] = uOld[N] - c * uOld[N] * dt / dx * (uOld[N+1] - uOld[N]);
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
