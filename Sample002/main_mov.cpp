/*
 * 移流方程式を解く．(周期境界条件)
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

// Rittai)個人的に、定数には static もつけます.
// mino)そうします．

//区間をN等分する．
static const int N = 512;
//円周率
// Rittai)ちなみに <cmath> に M_PI ってのがあります
static const double pi = acos(-1.);
//区間
static const double L = pi;
//境界条件
static const double a = 0.;
static const double b = 0.;
//速度
static const double c = 1.0;
//計算する時間
static const double tLimit = 100.0;
//描画のインターバル
static const int INTV = 10;

//初期条件
double func(double x) {
    //return std::exp(- 40* (x - L/2.) * (x - L/2.));
    if(x > L/2.) return 0.;
    return 1.;
}

int main()
{
    //分割した微小区間幅
    double dx = L / N;
    //微小時間
    double dt = 0.0001;
    //スタート地点（空間方向）
    double x  = 0.;
	
    // Rittai)配列より std::array. range-based for に index が欲しい！
	//        元コードの疑問点：普通の配列で十分なのに new double[ N + 1 ] するのは何故？
	// mino)
    //        配列と言えばそれしか思い浮かばなかった．std::array使います．
    std::array< double, N + 1 > uNew,uOld;
	
    //初期条件   
    for( int i=0 ; i<=N ; ++i ) {
        uOld[i] = func( x + i * dx );
    }

    //初期の時刻．
    double t = 0;

    //std::cout << (dt * c) / dx << std::endl;
    //安定性条件のチェック
    if( 1. < (dt * c) / dx ){
        std::cout << "安定性条件を満たしていません．" << std::endl;
        return 0;
    }

    /**********************************************************************/
    /*                 可視化の設定(gnuplot)                              */
    /**********************************************************************/
    std::FILE *gp = popen( "gnuplot -persist", "w" );
    fprintf(gp, "set xr [0:%f]\n", L);
    fprintf(gp, "set yr [-1:2]\n");
    fprintf(gp, "set size square\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "unset key\n");
    
    //初期条件描画
    fprintf(gp, "plot '-' w l\n");
    for( int i=0 ; i<=N; ++i ){
        fprintf(gp, "%f %f\n", x + i * dx , uOld[i]);
    }
    fprintf(gp, "e\n");
    fflush(gp);

    getchar();

    //タイムループ
    for( int it = 0 ; t<tLimit ; ++it ) {
        
        //uをdt時間だけ進める
        for( int i=1 ; i<N ; ++i ) {
            uNew[i] = uOld[i] - c * dt / ( 2. * dx ) * ( uOld[i+1] - uOld[i-1] );
        }
        //境界条件の計算（PBC）
        uNew[0] = uOld[0] - c * dt / ( 2. * dx ) * ( uOld[1] - uOld[N] );
        uNew[N] = uOld[N] - c * dt / ( 2. * dx ) * ( uOld[0] - uOld[N-1] );

        //uNewの描画
        if( it%INTV == 0 ){
            fprintf(gp, "plot '-' w l\n");
            for( int i=0 ; i<=N; ++i ){
                fprintf(gp, "%f %f\n", x + i * dx , uNew[i]);
            }
            fprintf(gp, "e\n");
            fflush(gp);
        }

        // mino) Rittai_3Dさんはこう書くようだ．
        std::copy( begin( uNew ), end( uNew ), begin( uOld ));
        
        //時刻の更新
        t += dt;
    }

    //FILEポインタの解放
    pclose(gp);
}
