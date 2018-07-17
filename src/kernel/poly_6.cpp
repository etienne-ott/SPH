#include "kernel/kernel.h"
#include "kernel/poly_6.h"
#include <cmath>

Poly6::Poly6(float h, int N, float mass) : Kernel(h, N, mass) {
    _fac2 = 8.f / (_h * _h * _h) * 3.f / (2.f * M_PI);
}

float Poly6::ValueOf(float r) {
    if (r > _h) return 0.f;

	return 315.0f / (64.0f * M_PI * std::pow(_h, 9)) * std::pow(_h * _h - r * r, 3);
}

void Poly6::FOD(float rx, float ry, float rz, float r, float* ret) {
    if (r >= _h) {
        ret[0] = 0.f;
        ret[1] = 0.f;
        ret[2] = 0.f;
        return;
    }

    float magn = -945.f / (32.f * M_PI * std::pow(_h, 9)) * std::pow(_h * _h - r * r, 3);
    ret[0] = magn * rx;
    ret[1] = magn * ry;
    ret[2] = magn * rz;
}
