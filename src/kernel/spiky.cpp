#include "kernel/kernel.h"
#include "kernel/spiky.h"
#include <cmath>

Spiky::Spiky(float h, int N, float mass) : Kernel(h, N, mass) {
    _fac2 = 8.f / (_h * _h * _h) * 3.f / (2.f * M_PI);
}

float Spiky::ValueOf(float r) {
    if (r > _h) return 0.f;

	return 15.0f / (M_PI * std::pow(_h, 6)) * std::pow(_h - r, 3);
}

void Spiky::FOD(float rx, float ry, float rz, float r, float* ret) {
    if (r >= _h) {
        ret[0] = 0.f;
        ret[1] = 0.f;
        ret[2] = 0.f;
        return;
    }

    float magn = -45.0f / (M_PI * std::pow(_h, 6)) * std::pow(_h - r, 2) / r;
    ret[0] = magn * rx;
    ret[1] = magn * ry;
    ret[2] = magn * rz;
}
