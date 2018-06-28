#include "kernel.hpp"
#include "gaussian.hpp"
#include <cmath>

Gaussian::Gaussian(float h, int N, float mass) : Kernel(h, N, mass) {

}

float Gaussian::ValueOf(float r) const {
    float q = r / _h;

    if (q > 1.0) {
        return 0.0;
    } else {
        return _fac1 * exp(-4.0 * q * q);
    }
}

void Gaussian::FOD(float rx, float ry, float rz, float r, float* ret) {
    float q = r / _h;

    if (q > 1.0) {
        ret[0] = 0;
        ret[1] = 0;
        ret[2] = 0;
    } else {
        float tmp = -8.0 * _fac1 * exp(-4.0 * q * q) / (_h * _h);

        ret[0] = tmp * rx;
        ret[1] = tmp * ry;
        ret[2] = tmp * rz;
    }
}

void Gaussian::SOD(float rx, float ry, float rz, float r, float* ret) {
    float q = r / _h;
    if (q > 1.0) {
        for (int i = 0; i < 9; i++) {
            ret[i] = 0.0;
        }
    } else {
        float tmp = 64 * _fac1 * exp(-4.0 * q * q) / (_h * _h * _h * _h);

        ret[0] = -8.0 * _fac1 * (exp(-4.0 * q * q) * (1.0 - rx) + rx / _fac1) / (_h * _h);
        ret[1] = tmp * rx * ry;
        ret[2] = tmp * rx * rz;
        ret[3] = ret[1];
        ret[4] = -8.0 * _fac1 * (exp(-4.0 * q * q) * (1.0 - ry) + ry / _fac1) / (_h * _h);
        ret[5] = tmp * ry * rz;
        ret[6] = ret[2];
        ret[7] = ret[5];
        ret[8] = -8.0 * _fac1 * (exp(-4.0 * q * q) * (1.0 - rz) + rz / _fac1) / (_h * _h);
    }
}