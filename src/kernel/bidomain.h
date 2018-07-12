#pragma once

#include "kernel/kernel.h"

class Bidomain: public Kernel {
public:
    /// @see Kernel::Kernel
    Bidomain(float h, int N, float mass);

    /// @see Kernel::ValueOf
    float ValueOf(float r) const;

    /// @see Kernel::FOD
    void FOD(float rx, float ry, float rz, float r, float* ret);
};
