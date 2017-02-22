#include "kernel.hpp"
#include <cmath>

Kernel::Kernel(double h, int N) {
	_h = h;
	_N = N;
	_fac1 = 4.0 * h * h * h * N;
}

double Kernel::Function(double r) {
	double v = r / _h;
    if (v >= 0.0 && v < 1.0) {
        v = (4.0 - 6.0 * r * r + 3.0 * r * r * r);
    } else if (v >= 1.0 && v < 2.0) {
        v = pow((2.0 - r), 3.0);
    } else {
        return 0.0;
    }
    return v / _fac1;
}

void Kernel::FOD(double rx, double ry, double rz, double r, double* ret) {
	double v = r / _h, ir = 1.0 / r;

    if (v >= 0.0 && v < 1.0) {
        ret[0] = -12.0 * rx * ir + 9.0 * rx * pow(ir, 2.0);
        ret[1] = -12.0 * ry * ir + 9.0 * ry * pow(ir, 2.0);
        ret[2] = -12.0 * rz * ir + 9.0 * rz * pow(ir, 2.0);
    } else if (v >= 1.0 && v < 2.0) {
        ret[0] = -3.0 * pow(2.0 - r, 2.0) * rx * ir;
        ret[1] = -3.0 * pow(2.0 - r, 2.0) * ry * ir;
        ret[2] = -3.0 * pow(2.0 - r, 2.0) * rz * ir;
    } else {
        ret[0] = 0.0;
        ret[1] = 0.0;
        ret[2] = 0.0;
    }

    ret[0] /= _fac1;
    ret[1] /= _fac1;
    ret[2] /= _fac1;
}