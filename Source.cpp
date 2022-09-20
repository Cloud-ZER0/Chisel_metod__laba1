#include "Matrix.h"
#include <iostream>

int main() {
	Matrix m{ {3.72, 3.47, 3.06}, {4.47, 4.10, 3.63}, {4.96, 4.53, 4.01} };
	std::vector<double> f{ 30.74, 36.80, 40.79 };
	m.add_free(std::move(f));
	auto res = m.Gaus();
	for (auto i : res) {
		std::cout << i << ' ';
	}
	return 0;
}