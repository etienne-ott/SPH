#pragma once

class Kernel {
public:
    /// Constructor.
    ///
    /// @param h float The range factor for the kernel
    /// @param N int The number of particles
    /// @param mass float The mass of a particle
    Kernel(float h, int N, float mass);

    /// Returns the value of the kernel evaluated with the given distance.
    ///
    /// @param r float The distance to evaluate
    /// @return float The value of the kernel function
    virtual float ValueOf(float r) const = 0;

    /// Returns the first order derivative of the kernel function evaluated with
    /// the given distance r and the individual components rx, ry and rz. The
    /// distance r could be calculated within the kernel, but since it is
    /// often necessarily calculated beforehand anyway, it is expected to be
    /// given as a parameter.
    ///
    /// @param rx float The x component of the distance vector
    /// @param ry float The y component of the distance vector
    /// @param rz float The z component of the distance vector
    /// @param r float The scalar distance value
    /// @param ret float* Output vector (3 dimensional)
    virtual void FOD(float rx, float ry, float rz, float r, float* ret) = 0;

    /// Interpolates the density at position (rx,ry,rz) using the kernel
    /// the kernel function.
    ///
    /// @param rx float The x component of the interpolation position
    /// @param ry float The y component of the interpolation position
    /// @param rz float The z component of the interpolation position
    /// @param density float* The density values of the particles
    /// @param position float* The position of the particles
    /// @return float The interpolated density value at position (rx,ry,rz)
    float InterpolateDensity(float rx, float ry, float rz, float* density, float* position) const;

    /// Sets the smoothing length uses by the kernel to the given value.
    ///
    /// @param h float The new smoothing length.
    void SetH(float h);

protected:
    /// @var _h float The range factor used by the kernel
    float _h;

    /// @var _N int The number of particles
    int _N;

    /// @var _mass float The mass of a particle
    float _mass;

    /// @var _fac1 float A precalculated factor used in the kernel
    float _fac1;
};
