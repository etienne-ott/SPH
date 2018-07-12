#include "kernel/kernel.h"
#include "kernel/bidomain.h"
#include <cmath>

Bidomain::Bidomain(float h, int N, float mass) : Kernel(h, N, mass) {

}

float Bidomain::ValueOf(float r) {
    float v = r / _h;
    if (v >= 0.0 && v < 1.0) {
        v = (4.0 - 6.0 * r * r + 3.0 * r * r * r);
    } else if (v >= 1.0 && v < 2.0) {
        v = pow((2.0 - r), 3.0);
    } else {
        return 0.0;
    }
    return v * _fac1;
}

void Bidomain::FOD(float rx, float ry, float rz, float r, float* ret) {
    float v = r / _h;

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
