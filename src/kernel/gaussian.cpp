#include "kernel/kernel.h"
#include "kernel/gaussian.h"
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
