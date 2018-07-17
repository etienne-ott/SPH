#include "kernel/kernel.h"
#include "kernel/wendland.h"
#include <cmath>

Wendland::Wendland(float h, int N, float mass) : Kernel(h, N, mass) {
    _fac2 = 8.f / (_h * _h * _h) * 3.f / (2.f * M_PI);
}

float Wendland::ValueOf(float r) {
    const float q = r / (0.5f * _h);

	if (q >= 2.f) return 0.f;

	const float h1 = 2.f / _h;
	const float tmp = 1.f - 0.5 * q;

	return 21.f / (16.f * M_PI)* h1 * h1 * h1
        * tmp * tmp * tmp * tmp * (2.f * q + 1.f);
}

void Wendland::FOD(float rx, float ry, float rz, float r, float* ret) {
    const float q = r / (0.5f * _h);

	if (q >= 2.f || q < 0.0001f) {
        ret[0] = 0.f;
        ret[1] = 0.f;
        ret[2] = 0.f;
        return;
    }

    float magn = 1.f - 0.5f * q;
	const float h1 = 2.f / _h;
	const float fac = 21.f / (16.f * M_PI) * h1 * h1 * h1;
	const float val = -5.f * q * magn * magn * magn * h1 / r;
	magn = val * fac;

    ret[0] = magn * rx;
    ret[1] = magn * ry;
    ret[2] = magn * rz;
}
