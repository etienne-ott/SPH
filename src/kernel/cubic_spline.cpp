#include "kernel.h"
#include "cubic_spline.h"
#include <cmath>

CubicSpline::CubicSpline(float h, int N, float mass) : Kernel(h, N, mass) {
    _fac2 = 8.f / (_h * _h * _h) * 3.f / (2.f * M_PI);
}

float CubicSpline::ValueOf(float r) const {
    float q = 2.f * r / _h;

    if (q >= 2.f) {
        return 0.f;
    } else if (q >= 1.f) {
        return _fac2 * (1.0f / 6.0f * (2.0f - q) * (2.0f - q) * (2.0f - q));
    } else {
        return _fac2 * (2.0f / 3.0f - q * q + 0.5f * q * q * q);
    }
}

void CubicSpline::FOD(float rx, float ry, float rz, float r, float* ret) {
    if (r >= _h || r < 0.0001f) {
        ret[0] = 0.0f;
        ret[1] = 0.0f;
        ret[2] = 0.0f;
        return;
    }

    float q = 2.f * r / _h;
    if (q >= 1.f) {
        ret[0] = _fac2 * 0.5f * (q - 2.f) * (q - 2.f) * rx / q;
        ret[1] = _fac2 * 0.5f * (q - 2.f) * (q - 2.f) * ry / q;
        ret[2] = _fac2 * 0.5f * (q - 2.f) * (q - 2.f) * rz / q;
        return;
    }

    ret[0] = _fac2 * -0.5f * q * (3.f * q - 4.f) * rx / q;
    ret[1] = _fac2 * -0.5f * q * (3.f * q - 4.f) * ry / q;
    ret[2] = _fac2 * -0.5f * q * (3.f * q - 4.f) * rz / q;
}

void CubicSpline::SOD(float rx, float ry, float rz, float r, float* ret) {
    return;
}