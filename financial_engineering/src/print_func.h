#include <iostream>
#include <vector>
#include <array>

template<typename T, typename U>
void print_pair(std::pair<T, U> p) {
	std::cout << p.first << ", " << p.second << std::endl;
}

template<typename T>
void print_array(const T a[], const size_t size) {
	if (size == 0)
		return;
	std::cout << '(';
	for (size_t i(0); i < size-1; i++) {
		std::cout << a[i] << ", ";
	}
	std::cout << a[size-1] << ')' << std::endl;
}

template<typename T, const size_t size>
void print_array(const std::array<T, size>& a) {
	print_array(a.data(), a.size());
}

template<typename T, const size_t size>
void print_array(const std::array<std::array<T, size>, size>& a) {
	for (const auto& el : a) {
		print_array(el.data(), el.size());
	}
}

template<typename T> 
void print_vector(const std::vector<T>& v) {
	std::cout << '(';
	for (std::size_t i(0) ; i < v.size() - 1; i++) {
		std::cout << v[i] << ", ";
	}
	std::cout << v[v.size() - 1] << ")\n";
}

template<typename T>
void print_vector(const std::vector<std::vector<T>>& v) {
	for (const auto& el : v) {
		print_vector(el);
	}
}