#ifndef __KERNEL_HPP
#define __KERNEL_HPP

class Kernel {
public:
    /// Constructor.
    ///
    /// @param h double The range factor for the kernel
    /// @param N int The number of particles
    Kernel(double h, int N);

    /// Returns the value of the kernel evaluated with the given distance.
    ///
    /// @param r double The distance to evaluate
    /// @return double The value of the kernel function
    double Function(double r) const;

    /// Returns the first order derivate of the kernel function evaluated with
    /// the given distance r and the individual components rx, ry and rz. The
    /// distance r could be calculated within the kernel, but since it is
    /// often necessarily calculated beforehand anyway, it is expected to be
    /// given as a parameter.
    ///
    /// @param rx double The x component of the distance vector
    /// @param ry double The y component of the distance vector
    /// @param rz double The z component of the distance vector
    /// @param r double The scalar distance value
    /// @param ret double* Output vector (3 dimensional)
    void FOD(double rx, double ry, double rz, double r, double* ret);

private:
    /// @var _h double The range factor used by the kernel
    double _h;

    /// @var _N int The number of particles
    int _N;

    /// @var _fac1 double A precalculated factor used in the kernel
    double _fac1;
};
#endif // __KERNEL_HPP
