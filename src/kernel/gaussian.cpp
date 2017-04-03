#include "kernel.hpp"
#include "gaussian.hpp"
#include <cmath>

Gaussian::Gaussian(double h, int N, double mass) : Kernel(h, N, mass) {

}

double Gaussian::ValueOf(double r) const {
    double q = r / _h;

    if (q > 1.0) {
        return 0.0;
    } else {
        return _fac1 * exp(-4.0 * q * q);
    }
}

void Gaussian::FOD(double rx, double ry, double rz, double r, double* ret) {
    double q = r / _h;

    if (q > 1.0) {
        ret[0] = 0;
        ret[1] = 0;
        ret[2] = 0;
    } else {
        double tmp = -8.0 * _fac1 * exp(-4.0 * q * q) / (_h * _h);

        ret[0] = tmp * rx;
        ret[1] = tmp * ry;
        ret[2] = tmp * rz;
    }
}

void Gaussian::SOD(double rx, double ry, double rz, double r, double* ret) {
    double q = r / _h;
    if (q > 1.0) {
        for (int i = 0; i < 9; i++) {
            ret[i] = 0.0;
        }
    } else {
        double tmp = 64 * _fac1 * exp(-4.0 * q * q) / (_h * _h * _h * _h);

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