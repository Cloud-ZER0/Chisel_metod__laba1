#pragma once

#include <vector>

class Matrix
{
private: 
	std::vector<std::vector<double>> mtr{};
	std::vector<double> free{};
	std::size_t size{};

	void helper(const std::size_t pos, const std::size_t pos2, const double value, bool sum);
public:
	Matrix() = default;

	explicit  Matrix(Matrix&& other) noexcept;

	explicit Matrix(const Matrix& other);

	explicit  Matrix(const std::initializer_list<std::initializer_list<double>>& li);

	Matrix& operator = (Matrix&& other) noexcept;

	void add_free(std::vector<double>&& vec) noexcept;

	std::vector<double> Gaus();

	std::vector <double> Kholetskiy();

	~Matrix() = default;
};

