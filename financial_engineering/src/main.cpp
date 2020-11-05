#include "math_func.h"
#include "print_func.h"
#include "random_generator.h"

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

int main(int, char**) {

    test_cholesky();

    return 0;
}
