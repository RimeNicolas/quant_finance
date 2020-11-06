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

void test_brownian_motion_d() {
	const size_t n(11);
	const size_t d(2);
	BrownianMotionMultiDim<double> bm(10.0, n);
	std::vector<double> mu = { 0.0,0 };
	std::vector<std::vector<double>> cov(d, std::vector<double>(d,0));
	cov[0][0] = 4;
	cov[0][1] = 0.0;
	cov[1][0] = 0.0;
	cov[1][1] = 1;
	std::vector<std::vector<double>> v2 = bm.motion(mu, cov);
	for (size_t i(0); i < n; i++) {
		print_vector(v2[i]);
	}
}

int main(int, char**) {

	test_brownian_motion_d();


    return 0;
}
