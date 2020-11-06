#include <utility>
#include <cmath>
#include <vector>
#include <iostream>
#include <array>

#include "math_func.h"
#include "print_func.h"
#include "random_generator.h"
#include "motion.h"

void test_cholesky() {
	constexpr size_t n = 3;
	std::vector<std::vector<double>> cov_matrix(n, std::vector<double>(n,0));
	cov_matrix[0][0] = 4;
	cov_matrix[0][1] = 12;
	cov_matrix[0][2] = -16;
	cov_matrix[1][0] = 12;
	cov_matrix[1][1] = 37;
	cov_matrix[1][2] = -43;
	cov_matrix[2][0] = -16;
	cov_matrix[2][1] = -43;
	cov_matrix[2][2] = 98;

	print_vector(cov_matrix);

	std::vector<std::vector<double>> low_triangular = \
		cholesky_factorisation(cov_matrix);

	std::cout << std::endl;
	print_vector(low_triangular);
}

// void test_covariance_matrix() {
// 	RandomGeneratorLewis<double> randgen(0);

// 	std::array<double, 3> mu = { 0,0,0 };
// 	Tensor2<double, 3> cov_matrix = { {0} };
// 	cov_matrix[0][0] = 4;
// 	cov_matrix[0][1] = 12;
// 	cov_matrix[0][2] = -16;
// 	cov_matrix[1][0] = 12;
// 	cov_matrix[1][1] = 37;
// 	cov_matrix[1][2] = -43;
// 	cov_matrix[2][0] = -16;
// 	cov_matrix[2][1] = -43;
// 	cov_matrix[2][2] = 98;

// 	print_array(cov_matrix);

// 	Tensor2<double, 3> low_triangular = \
// 		cholesky_factorisation(cov_matrix);

// 	std::cout << std::endl;
// 	print_array(low_triangular);

// 	std::array<double, 3> z = randgen.normal_array(mu, cov_matrix);

// 	print_array(z);
// }


void test_brownian_motion_1d() {
	const size_t n_steps(10);
	const double final_time(10.0);
	BrownianMotion<double> bm(final_time, n_steps);
	std::vector<double> v1 = bm.motion(1,4);
	std::vector<double> mu = { 1,0,0,1,0 };
	std::vector<double> cov = { 1,2,2,1,2 };
	std::vector<double> v2 = bm.motion(mu, cov);
	print_vector(v1);
	print_vector(v2);
}

int main(int, char**) {

	test_brownian_motion_1d();


    return 0;
}
