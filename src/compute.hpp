#include "kernel.hpp"

#ifndef __COMPUTE_HPP
#define __COMPUTE_HPP

class Compute {
public:
    /// Constructor.
    ///
    /// @param N int The number of parameters
    /// @param kernel Kernel* The kernel used for calculations
    Compute(int N, Kernel* kernel);

    /// Destructor. Destroys the data fields properly, that were created
    /// during initialization.
    ~Compute();

    /// Calculates one timestep of the fluid simulations, which includes
    /// calculating the forces on each particles, enforcing boundary
    /// conditions, then integrating the particle velocities and positions
    /// to the next state.
    void Timestep();

    /// Returns the particle positions as consecutive x, y and z components,
    /// for an overall number of 3*N doubles.
    /// 
    /// @return double* The particle positions
    double* GetPosition();

    /// Returns the particle densities as consecutive values for each particle.
    /// 
    /// @return double* The particle densities
    double* GetDensity();

private:
    /// @var _N int The number of particles.
    int _N;

    /// @var _kernel Kernel* The kernel to use for calculations.
    Kernel* _kernel;

    /// @var _vec1 double* A temporary vector used in calculations.
    double* _vec1;

    /// @var _vec1 double* A temporary vector used in calculations.
    double* _vec2;

    /// @var _position double* The particle positions in x, y and z coordinates.
    double* _position;

    /// @var _velocity double* The particle velocities in x, y and z components.
    double* _velocity;

    /// @var _force double* The sum of all forces acting on the particles in x,
    ///     y and z components.
    double* _force;

    /// @var _density double* The density of the fluid particles.
    double* _density;

    /// @var _pressure double* The pressure of the fluid particles themselves,
    ///     as opposed to the pressure acting on them by other particles.
    double* _pressure;
};
#endif // __COMPUTE_HPP
