#pragma once
#include <iostream>
#include <cmath>
#include <array>

namespace rittai3d{
	namespace utility{
		// [ minimum, maximum ) の範囲でラップアラウンド
		template <typename T>
		constexpr T wrap_around(T value, T minimum, T maximum){
			const T n = (value - minimum) % (maximum - minimum);
			return n >= 0 ? (n + minimum) : (n + maximum); 
		}
	}
}

namespace sksat {
	template<std::size_t Num, typename T = double>
	class array_wrapper{
	private:
		std::array<T, Num> arr;
	public:
		constexpr array_wrapper() : arr() {}
		~array_wrapper() = default;

		T& operator[](int i){
			i = rittai3d::utility::wrap_around(i, 0, static_cast<int>(Num));
			return arr[i];
		}
	};
}
