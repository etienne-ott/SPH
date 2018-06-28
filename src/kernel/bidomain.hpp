#include "kernel.hpp"

#ifndef __KERNEL_BIDOMAIN_HPP
#define __KERNEL_BIDOMAIN_HPP

class Bidomain: public Kernel {
public:
    /// @see Kernel::Kernel
    Bidomain(float h, int N, float mass);

    /// @see Kernel::ValueOf
    float ValueOf(float r) const;

    /// @see Kernel::FOD
    void FOD(float rx, float ry, float rz, float r, float* ret);

    /// @see Kernel::SOD
    void SOD(float rx, float ry, float rz, float r, float* ret);
};
#endif // __KERNEL_BIDOMAIN_HPP
