/*
 * 拡散方程式を解く．
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
static const int N = 100;
//円周率
// Rittai)ちなみに <cmath> に M_PI ってのがあります
static const double pi = acos(-1.);
//区間
static const double L = pi;
//境界条件
static const double a = 0.;
static const double b = 0.;
//拡散係数
static const double D = 1.;
//計算する時間
static const double tLimit = 1000.0;

//初期条件
double func(double x) {
    return std::sin(x);
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

    //境界条件
    uOld[0] = a;
    uOld[N] = b;

    //初期の時刻．
    double t = 0;

    //安定性条件のチェック
    if( D * dt / ( dx * dx ) > 0.5 ){
        std::cout << "安定性条件を満たしていません．" << std::endl;
        return 0;
    }

    //タイムループ
    for( int it = 0 ; t<tLimit ; ++it ) {
        
        //uをdt時間だけ進める
        for( int i=1 ; i<N ; ++i ) {
            uNew[i] = uOld[i] + D * dt / ( dx * dx ) * (uOld[i-1] - 2. * uOld[i] + uOld[i+1] );
        }

        // uNew -> uOld へ更新
        /*for( int i=0 ; i<=N ; ++i ) {
            uOld[i] = uNew[i];
        }*/
        // mino) Rittai_3Dさんはこう書くようだ．
        std::copy( begin( uNew ), end( uNew ), begin( uOld ));
        
        //時刻の更新
        t += dt;
    }
}
