constexpr double m_pi = 3.14159265358979323846;

template<typename T>
std::vector<std::vector<T>> cholesky_factorisation(const std::vector<std::vector<T>>& matrix);

template<typename T>
class RandomGenerator {
public:
	RandomGenerator(unsigned int seed, const unsigned int mod,
		const unsigned int a, const unsigned int c);
	unsigned int randbasics();
	T randuniform01();

protected:
	unsigned int seed;
	const unsigned int mod;
	T mod_inv;
	const unsigned int a;
	const unsigned int c;
};

template<typename T>
inline RandomGenerator<T>::RandomGenerator(unsigned int seed, const unsigned int mod,
	const unsigned int a, const unsigned int c) 
	: mod(mod), a(a), c(c) {
	this->seed = seed & this->mod;
	this->mod_inv = 1.0 / static_cast<T>(this->mod);
}

template<typename T>
inline unsigned int RandomGenerator<T>::randbasics() {
	seed = (seed * a + c) & mod;
	return seed;
}

template<typename T>
T RandomGenerator<T>::randuniform01() {
	return randbasics() * this->mod_inv;
}

//==================================================

template<typename T>
class RandomGeneratorLewis : public RandomGenerator<T> {
public:
	RandomGeneratorLewis(unsigned int seed);
	T random_exponential(const T theta);
	std::pair<T, T> box_muller();
	std::pair<T, T> box_muller_opt();

	std::vector<T> gaussian_vector(const std::vector<T>& mu, const std::vector<std::vector<T>>& cov);
};

template<typename T>
inline RandomGeneratorLewis<T>::RandomGeneratorLewis(unsigned int seed)
	: RandomGenerator<T>(seed, 0x7fffffff, 16807, 12345U) {}

template<typename T>
T inline RandomGeneratorLewis<T>::random_exponential(const T theta) {
	return -theta * std::log(this->randuniform01());
}

template<typename T>
std::pair<T, T> inline RandomGeneratorLewis<T>::box_muller() {
	const T r = -2 * std::log(this->randuniform01());
	const T v = 2 * m_pi * this->randuniform01();
	return std::pair<T, T>(std::sqrt(r)*std::cos(v), std::sqrt(r)*std::sin(v));
}

template<typename T>
std::pair<T, T> inline RandomGeneratorLewis<T>::box_muller_opt() {
	T x(2), u1(0), u2(0);
	while (x > 1) {
		u1 = 2 * this->randuniform01() - 1;
		u2 = 2 * this->randuniform01() - 1;
		x = u1 * u1 + u2 * u2;
	}
	const double y(std::sqrt(-2 * std::log(x) / x));
	return std::pair<T, T>(u1 * y, u2 * y);
}

template<typename T>
inline std::vector<T> RandomGeneratorLewis<T>::gaussian_vector(const std::vector<T>& mu, const std::vector<std::vector<T>>& cov){
	std::vector<std::vector<T>> A = cholesky_factorisation<T>(cov);
	std::vector<T> random_array = mu;
	std::size_t dim = mu.size();
	T z1(0);
	std::pair<T, T> z = { 0,0 };
	for (size_t i(0); i < dim; i++) {
		if (i % 2 == 0) {
			z = this->box_muller();
			z1 = z.first;
		}
		else {
			z1 = z.second;
		}
		for (size_t j(0); j < dim; j++) {
			random_array[j] += A[j][i] * z1;
		}
	}
	return random_array;
}