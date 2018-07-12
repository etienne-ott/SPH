#include "kernel/kernel.h"
#include <cmath>

Kernel::Kernel(float h, int N, float mass) {
    _h = h;
    _N = N;
    _mass = mass;
    _fac1 = 1.0 / (h * h * h);
}

float Kernel::InterpolateDensity(float rx, float ry, float rz, float* density, float* position) {
    float sum = 0.0;
    float distance = 0.0;

    for (int j = 0; j < _N; j++) {
        distance = pow(
            (position[j * 3] - rx) * (position[j * 3] - rx)
                + (position[j * 3 + 1] - ry) * (position[j * 3 + 1] - ry)
                + (position[j * 3 + 2] - rz) * (position[j * 3 + 2] - rz),
            0.5
        );
        sum += _mass * density[j] * this->ValueOf(distance);
    }

    return sum;
}

void Kernel::SetH(float h) {
    _h = h;
    _fac1 = 1.0 / (h * h * h);
}