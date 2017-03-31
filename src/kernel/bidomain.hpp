#include "kernel.hpp"

#ifndef __KERNEL_BIDOMAIN_HPP
#define __KERNEL_BIDOMAIN_HPP

class Bidomain: public Kernel {
public:
    /// @see Kernel::Kernel
    Bidomain(double h, int N, double mass);

    double Function(double r) const;

    void FOD(double rx, double ry, double rz, double r, double* ret);

    void SOD(double rx, double ry, double rz, double r, double* ret);
};
#endif // __KERNEL_BIDOMAIN_HPP
