#include "kernel.hpp"

#ifndef __KERNEL_GAUSSIAN_HPP
#define __KERNEL_GAUSSIAN_HPP

class Gaussian: public Kernel {
public:
    /// @see Kernel::Kernel
    Gaussian(double h, int N, double mass);

    /// @see Kernel::ValueOf
    double ValueOf(double r) const;

    /// @see Kernel::FOD
    void FOD(double rx, double ry, double rz, double r, double* ret);

    /// @see Kernel::SOD
    void SOD(double rx, double ry, double rz, double r, double* ret);
};
#endif // __KERNEL_GAUSSIAN_HPP
