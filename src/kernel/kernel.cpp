#include "kernel.hpp"
#include <cmath>

Kernel::Kernel(double h, int N, double mass) {
    _h = h;
    _N = N;
    _mass = mass;
    _fac1 = 1.0 / (4.0 * h * h * h * N);
}

double Kernel::InterpolateDensity(double rx, double ry, double rz, double* density, double* position) const {
    double sum = 0.0;
    double distance = 0.0;

    for (int j = 0; j < _N; j++) {
        distance = pow(
            (position[j * 3] - rx) * (position[j * 3] - rx)
                + (position[j * 3 + 1] - ry) * (position[j * 3 + 1] - ry)
                + (position[j * 3 + 2] - rz) * (position[j * 3 + 2] - rz),
            0.5
        );
        sum += _mass * density[j] * this->Function(distance);
    }

    return sum;
}

void Kernel::SetH(double h) {
    _h = h;
    _fac1 = 1.0 / (4.0 * h * h * h * _N);
}