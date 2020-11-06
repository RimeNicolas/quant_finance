#include <vector>
#include <utility>

template<typename T>
class BrownianMotion {
public:
	BrownianMotion(const T final_time, const std::size_t n_steps);
	std::vector<T> motion(const T mu, const T sigma);
	std::vector<T> motion(const std::vector<T>& mu, const std::vector<T>& sigma);

protected:
	T final_time;
	T delta_t;
	T sqrt_delta_t;
	const size_t n_steps;
	RandomGeneratorLewis<T> rgl;
};

template<typename T>
inline BrownianMotion<T>::BrownianMotion(const T final_time, const std::size_t n_steps)
	: final_time(final_time), n_steps(n_steps), rgl(0) {
	this->delta_t = this->final_time / (this->n_steps - 1);
	std::cout << "delta t = " << delta_t << std::endl;
	this->sqrt_delta_t = std::sqrt(this->delta_t);
}

template<typename T>
std::vector<T> inline BrownianMotion<T>::motion(const T mu, const T sigma) {
	std::vector<T> v(this->n_steps, 0);
	std::pair<T, T> z_pair = { 0,0 };
	T z(0);
	for (size_t i(1); i < v.size(); i++) {
		if (i % 2 == 1) {
			z_pair = this->rgl.box_muller_opt();
			z = z_pair.first;
		}
		else {
			z = z_pair.second;
		}
		v[i] = v[i - 1] + mu * delta_t + sigma * sqrt_delta_t * z;
	}
	return v;
}

template<typename T>
inline std::vector<T> BrownianMotion<T>::motion(const std::vector<T>& mu, const std::vector<T>& sigma) {
	std::vector<T> v(this->n_steps, 0);
	std::pair<T, T> z_pair = { 0,0 };
	T z(0);
	for (size_t i(1); i < v.size(); i++) {
		if (i % 2 == 1) {
			z_pair = this->rgl.box_muller_opt();
			z = z_pair.first;
		}
		else {
			z = z_pair.second;
		}
		v[i] = v[i - 1] + mu[i] * delta_t + sigma[i] * sqrt_delta_t * z;
	}
	return v;
}

template<typename T>
class BrownianMotionMultiDim : public BrownianMotion<T> {
public:
	BrownianMotionMultiDim(const T final_time);
	// std::vector<std::vector<T>> motion(const std::vector<T> mu, const std::vector<T> sigma);
};

template<typename T>
inline BrownianMotionMultiDim<T>::BrownianMotionMultiDim(const T final_time) 
	: BrownianMotion<T>(final_time) {}

// template<typename T>
// inline Tensor2d<T> BrownianMotionMultiDim<T, n, d>::motion(const Tensor1<T, d> mu, const Tensor2<T, d> sigma) {
// 	Tensor1<T, d> z = { 0 };
// 	Tensor2d<T, n, d> v = { 0 };
// 	for (size_t i(1); i < n; i++) {
// 		z = rgl.gauss_array(mu, sigma);
// 		for (size_t j(0); j < d; j++) {
// 			v[i][j] = v[i - 1][j] + mu[j] * delta_t + sqrt_delta_t * z[j];
// 		}
// 	}
// 	return v;
// }