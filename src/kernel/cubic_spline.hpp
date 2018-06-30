#include "kernel.hpp"

#ifndef __KERNEL_CUBIC_SPLINE_HPP
#define __KERNEL_CUBIC_SPLINE_HPP

class CubicSpline: public Kernel {
private:
    float _fac2;
public:
    /// @see Kernel::Kernel
    CubicSpline(float h, int N, float mass);

    /// @see Kernel::ValueOf
    float ValueOf(float r) const;

    /// @see Kernel::FOD
    void FOD(float rx, float ry, float rz, float r, float* ret);

    /// @see Kernel::SOD
    void SOD(float rx, float ry, float rz, float r, float* ret);
};
#endif // __KERNEL_CUBIC_SPLINE_HPP
