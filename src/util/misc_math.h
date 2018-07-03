#pragma once

#include <tuple>
#include <cmath>

inline float dot_product(float *v_1, float *v_2)
{
	return *(v_1+0) * *(v_2+0) + *(v_1+1) * *(v_2+1) + *(v_1+2) * *(v_2+2);
}

inline float dot_product(float *v_1, std::tuple<float, float, float> v_2)
{
	return *(v_1+0) * std::get<0>(v_2) + *(v_1+1)
		* std::get<1>(v_2) + *(v_1+2) * std::get<2>(v_2);
}

inline float dot_product(std::tuple<float, float, float> v_1, float *v_2)
{
	return *(v_2+0) * std::get<0>(v_1) + *(v_2+1)
		* std::get<1>(v_1) + *(v_2+2) * std::get<2>(v_1);
}

inline float dot_product(std::tuple<float, float, float> v_1, std::tuple<float, float, float> v_2)
{
	return std::get<0>(v_1) * std::get<0>(v_2) + std::get<1>(v_1)
		* std::get<1>(v_2) + std::get<2>(v_1) * std::get<2>(v_2);
}

inline std::tuple<float, float, float> difference_of_vectors(float *v_1, float *v_2)
{
	return std::make_tuple(*v_1-*v_2,*(v_1+1)-*(v_2+1),*(v_1+2)-*(v_2+2));
}

// ==========================================================
// ========= Comute power of float at compiletime ==========
// ==========================================================
template <unsigned int n>
constexpr float get_power(float x)
{
    return x * get_power<n-1>(x);
}

// We ignore the unused parameter warning
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
template <>
constexpr float get_power<0>(float x)
{
    return 1.0;
}
#pragma GCC diagnostic pop

/// Fast Sqrt with floating point precission not really precise
float inline fastSqrt1(const float& x)
{
	unsigned int i = *(unsigned int*) &x;
	// adjust bias
	i  += 127 << 23;
	// approximation of square root
	i >>= 1;
	return *(float*) &i;
}

/// Slower than fastSqrt1 but almost that precise than std::sqrt and 30% faster
float inline fastSqrt2(const float& x)
{
	union
	{
		int i;
		float x;
	} u;
	u.x = x;
	u.i = (1 << 29) + (u.i >> 1) - (1 << 22);

	// Two Babylonian Steps (simplified from:)
	// u.x = 0.5f * (u.x + x/u.x);
	// u.x = 0.5f * (u.x + x/u.x);
	u.x = u.x + x/u.x;
	u.x = 0.25f * u.x + x/u.x;

	return u.x;
}
