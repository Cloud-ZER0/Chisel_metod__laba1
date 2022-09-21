#include "Matrix.h"
#include <iostream>
Matrix::Matrix(Matrix&& other) noexcept : mtr(std::move(other.mtr)), size(std::move(other.size)) {}

Matrix::Matrix(const Matrix& other): mtr(other.mtr), size(other.size) {}

Matrix::Matrix(const std::initializer_list<std::initializer_list<double>>& li) {
	std::size_t size = 0;
	std::size_t sz = 0;
	for (auto it = li.begin(); it != li.end(); ++it) {
		mtr.emplace_back(it->size(), 0);
		for (auto jt = it->begin(); jt != it->end(); ++jt) {
			mtr[size].at(sz) = *jt;
			++sz;
		}
		++size;
		sz = 0;
	}
	for (std::size_t i = 0; i < mtr.size(); ++i) {
		this->size += mtr[i].size();
	}
}	
Matrix& Matrix::operator = (Matrix&& other) noexcept {
	mtr = std::move(other.mtr);
	size = std::move(other.size);
	return *this;
}
void Matrix::add_free(std::vector<double>&& vec) noexcept {
	free = std::move(vec);
}

void Matrix::helper(const std::size_t pos1, const std::size_t pos2, const double value, bool sum = false) {

	if (pos1 > mtr[0].size() || pos2 > mtr[0].size()) {
		std::cerr << "Invalid position\n", exit(-5);
	}
	if (size == 0) {
		std::cerr << "Matrix is emty\n", exit(-1);
	}

	for (std::size_t i = 0; i < mtr[0].size(); ++i) {
		auto val = mtr[pos2].at(i);
		if (!sum) {
			mtr[pos1].at(i) -= (val * abs(value));
		}
		else {
			mtr[pos1].at(i) += (val * abs(value));
		}
	}
}

std::vector<double> Matrix::Gaus() {
	
	for (std::size_t i = 0; i < mtr.size(); ++i) {
		mtr[i].push_back(free[i]);
	}

	auto decl = [this](std::size_t pos) -> void {
		double val = mtr[pos].at(pos);
		for (std::size_t i = 0; i < mtr[pos].size(); ++i) {
			if (mtr[pos].at(i) != 0) {
				mtr[pos].at(i) /= val;
			}
		}
	};

	decl(0);
	for (std::size_t i = 1; i < mtr.size(); ++i) {
		helper(i, 0, mtr[i].at(0));
	}

	decl(1);
	helper(0, 1, mtr[0].at(1), ((mtr[0].at(1)) > 0 ? false : true));
	helper(2, 1, mtr[2].at(1), ((mtr[2].at(1)) > 0 ? false : true)); //true

	decl(2);
	helper(0, 2, mtr[0].at(2), ((mtr[0].at(2)) > 0 ? false : true)); // true
	helper(1, 2, mtr[1].at(2), ((mtr[1].at(2)) > 0 ? false : true));

	return std::vector<double>{mtr[0].at(3), mtr[1].at(3), mtr[2].at(3)};
}

std::vector<double> Matrix::Kholetskiy() {	

	std::size_t sz = mtr[0].size();
	auto decl = [this, sz](std::size_t pos1, std::size_t pos2, double val) {
		for (pos1; pos1 < sz; ++pos1) {
			mtr[pos2].at(pos1) /= val;
		}
	};
	decl(1, 0, mtr[0].at(0)); 

	for (std::size_t i = 1; i < sz; ++i) {
		for (std::size_t j = 1; j < sz; ++j) {
			if (i == 2 and j == 2) {
				auto v1 = mtr[i].at(0);
				auto v2 = mtr[0].at(j);
				auto v3 = mtr[i - 1].at(j);
				auto v4 = mtr[i].at(j - 1);
				mtr[i].at(j) -= (mtr[i].at(0) * mtr[0].at(j) + mtr[i - 1].at(j) * mtr[i].at(j - 1));
				break;
			}
			mtr[i].at(j) -= mtr[i].at(0) * mtr[0].at(j);
		}
		if (i + 1 < sz) {
			decl(i + 1, i, mtr[i].at(i)); 
		}
	}


	Matrix C(*this);

	for (std::size_t i = 1; i < sz; ++i) {
		mtr[0].at(i) = 0;
	}
	mtr[1].at(2) = 0;

	for (std::size_t i = 0; i < sz; ++i) {
		C.mtr[i].at(i) /= C.mtr[i].at(i);
	}

	for (std::size_t i = 1; i < sz; ++i) {
		for (std::size_t j = 0; j < sz; ++j) {
			if (C.mtr[i].at(j) != 1 and C.mtr[i].at(j) != C.mtr[1].at(2)) {
				C.mtr[i].at(j) = 0;
			}
		}
	}
	
	std::vector<double> y{};
	y.reserve(3);

	y.push_back(free[0] / mtr[0].at(0));
	y.push_back((free[1] - mtr[1].at(0) * y[0]) / mtr[1].at(1));
	y.push_back((free[2] - mtr[2].at(0) * y[0] - mtr[2].at(1) * y[1]) / mtr[2].at(2));

	std::vector<double> x{};
	x.reserve(3);

	x.push_back(y[2]);
	x.push_back(y[1] - (C.mtr[1].at(2) * y[2]));
	x.push_back(y[0] -(C.mtr[0].at(1) * x[1] + C.mtr[0].at(2) * y[2]));

	std::reverse(x.begin(), x.end());

	return x;
}