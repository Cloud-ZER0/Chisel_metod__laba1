#include "Matrix.h"
#include <iostream>

int main() {
	//Matrix den{ {2.16, 1.96, 1.56}, {3.55, 3.23, 2.78}, {4.85, 4.47, 3.97} };
	//Matrix chel{ {3.78, 3.44, 3.02}, {4.33, 3.88, 3.39}, {4.76, 4.24, 3.71} };
	Matrix me{ {3.72, 3.47, 3.06}, {4.47, 4.10, 3.63}, {4.96, 4.53, 4.01} };
	std::vector<double> me_f{ 30.74, 36.8, 40.79 };
	//std::vector<double> den_f{ 13.16, 21.73, 29.75 };
	//std::vector<double> chel_f{ 46.81, 53.43, 58.73 };
	me.add_free(std::move(me_f));
	auto res = me.Gaus();
	int x = 1;
	for (auto i : res) {
		std::cout <<"x" << x << '=' << i << '\n';
		++x;
	}
	return 0;
}