#ifndef __KERNEL_HPP
#define __KERNEL_HPP

class Kernel {
public:
    /// Constructor.
    ///
    /// @param h double The range factor for the kernel
    /// @param N int The number of particles
    /// @param mass double The mass of a particle
    Kernel(double h, int N, double mass);

    /// Returns the value of the kernel evaluated with the given distance.
    ///
    /// @param r double The distance to evaluate
    /// @return double The value of the kernel function
    double Function(double r) const;

    /// Returns the first order derivative of the kernel function evaluated with
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

    /// Returns the second order derivative of the kernel function evaluated with
    /// the given distance r and the individual components rx, ry and rz. The
    /// distance r could be calculated within the kernel, but since it is
    /// often necessarily calculated beforehand anyway, it is expected to be
    /// given as a parameter.
    ///
    /// @param rx double The x component of the distance vector
    /// @param ry double The y component of the distance vector
    /// @param rz double The z component of the distance vector
    /// @param r double The scalar distance value
    /// @param ret double* Output vector (9 dimensional)
    void SOD(double rx, double ry, double rz, double r, double* ret);

    /// Interpolates the density at position (rx,ry,rz) using the kernel
    /// the kernel function.
    ///
    /// @param rx double The x component of the interpolation position
    /// @param ry double The y component of the interpolation position
    /// @param rz double The z component of the interpolation position
    /// @param density double* The density values of the particles
    /// @param position double* The position of the particles
    /// @return double The interpolated density value at position (rx,ry,rz)
    double InterpolateDensity(double rx, double ry, double rz, double* density, double* position) const;

private:
    /// @var _h double The range factor used by the kernel
    double _h;

    /// @var _N int The number of particles
    int _N;

    /// @var _mass double The mass of a particle
    double _mass;

    /// @var _fac1 double A precalculated factor used in the kernel
    double _fac1;
};
#endif // __KERNEL_HPP
