#pragma once

#include "kernel/kernel.h"

class Wendland: public Kernel {
private:
    float _fac2;
public:
    /// @see Kernel::Kernel
    Wendland(float h, int N, float mass);

    /// @see Kernel::ValueOf
    float ValueOf(float r);

    /// @see Kernel::FOD
    void FOD(float rx, float ry, float rz, float r, float* ret);
};
