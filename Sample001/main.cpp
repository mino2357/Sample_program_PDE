/*
 * 拡散方程式を解く．
 *
 * Takaaki MINOMO
 * 私のコメントには先頭に mino) がつきます
 *
 * 改悪:Rittai_3D（仮名）
 *  わたしのコメントには先頭に Rittai) がつきます
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

int main()
{
    //分割した微小区間幅
    double dx = L / N;
    //スタート地点（空間方向）
    double x  = 0.;
	
    // Rittai)配列より std::array. range-based for に index が欲しい！
	//        元コードの疑問点：普通の配列で十分なのに new double[ N + 1 ] するのは何故？
	// mino)
    //        配列と言えばそれしか思い浮かばなかった．std::array使います．
    std::array< double, N + 1 > u;
	
    //初期条件   
    for( int i=0 ; i<=N ; ++i ) {
        u[i] = std::sin( x + i * dx );
    }

    //境界条件
    u[0] = a;
    u[N] = b;


    //デバッグ用出力
    for( int i=0 ; i<=N ; ++i ) {
        std::cout << x + i * dx << " " << u[i] << std::endl;
    }
}

