#include "kernel.hpp"
#include <cmath>

Kernel::Kernel(double h, int N, double mass) {
    _h = h;
    _N = N;
    _mass = mass;
    _fac1 = 1.0 / (4.0 * h * h * h * N);
}

double Kernel::Function(double r) const {
    double v = r / _h;
    if (v >= 0.0 && v < 1.0) {
        v = (4.0 - 6.0 * r * r + 3.0 * r * r * r);
    } else if (v >= 1.0 && v < 2.0) {
        v = pow((2.0 - r), 3.0);
    } else {
        return 0.0;
    }
    return v * _fac1;
}

void Kernel::FOD(double rx, double ry, double rz, double r, double* ret) {
    double v = r / _h;

    if (v >= 0.0 && v < 1.0) {
        ret[0] = -12.0 * rx  + 9.0 * rx * r;
        ret[1] = -12.0 * ry  + 9.0 * ry * r;
        ret[2] = -12.0 * rz  + 9.0 * rz * r;
    } else if (v >= 1.0 && v < 2.0) {
        ret[0] = -3.0 * rx * (4 / r + r - 4);
        ret[1] = -3.0 * ry * (4 / r + r - 4);
        ret[2] = -3.0 * rz * (4 / r + r - 4);
    } else {
        ret[0] = 0.0;
        ret[1] = 0.0;
        ret[2] = 0.0;
        return;
    }

    ret[0] *= _fac1;
    ret[1] *= _fac1;
    ret[2] *= _fac1;
}

void Kernel::SOD(double rx, double ry, double rz, double r, double* ret) {
    double v = r / _h, ir = 1.0 / r;

    if (v >= 0.0 && v < 1.0) {
        ret[0] = -12 + 9 * r + 9 * rx * rx * ir;
        ret[1] = 9 * ry * ir;
        ret[2] = 9 * rz * ir;
        ret[3] = 9 * rx * ir;
        ret[4] = -12 + 9 * r + 9 * ry * ry * ir;
        ret[5] = 9 * rz * ir;
        ret[6] = 9 * rx * ir;
        ret[7] = 9 * ry * ir;
        ret[8] = -12 + 9 * r + 9 * rz * rz * ir;
    } else if (v >= 1.0 && v < 2.0) {
        ret[0] = -12 * ir + 12 - 15 * r - 3 * rx * rx * ir;
        ret[1] = -3 * rx * (4 * r / ry + ry * ir);
        ret[2] = -3 * rx * (4 * r / rz + rz * ir);
        ret[3] = -3 * ry * (4 * r / rx + rx * ir);
        ret[4] = -12 * ir + 12 - 15 * r - 3 * ry * ry * ir;
        ret[5] = -3 * ry * (4 * r / rz + rz * ir);
        ret[6] = -3 * rz * (4 * r / rx + rx * ir);
        ret[7] = -3 * rz * (4 * r / ry + ry * ir);
        ret[8] = -12 * ir + 12 - 15 * r - 3 * ry * ry * ir;
    } else {
        for (int i = 0; i < 9; i++) {
            ret[i] = 0.0;
        }
        return;
    }

    for (int i = 0; i < 9; i++) {
        ret[i] *= _fac1;
    }
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