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
	BrownianMotionMultiDim(const T final_time, const std::size_t n_steps);
	std::vector<std::vector<T>> motion(const std::vector<T>& mu, const std::vector<std::vector<T>>& sigma);
};

template<typename T>
inline BrownianMotionMultiDim<T>::BrownianMotionMultiDim(const T final_time, const std::size_t n_steps) 
	: BrownianMotion<T>(final_time, n_steps) {}

template<typename T>
inline std::vector<std::vector<T>> BrownianMotionMultiDim<T>::motion(const std::vector<T>& mu, const std::vector<std::vector<T>>& sigma) {
    T dim = mu.size();
	std::vector<T> z(dim, 0);
    std::vector<std::vector<T>> v(this->n_steps, std::vector<T>(dim, 0));
	for (size_t i(1); i < this->n_steps; i++) {
		z = this->rgl.gaussian_vector(std::vector<T>(dim,0), sigma);
		for (size_t j(0); j < dim; j++) {
			v[i][j] = v[i - 1][j] + mu[j] * this->delta_t + this->sqrt_delta_t * z[j];
		}
	}
	return v;
}