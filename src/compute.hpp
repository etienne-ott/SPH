#include "kernel.hpp"
#include "parameter.hpp"

#ifndef __COMPUTE_HPP
#define __COMPUTE_HPP

class Compute {
public:
    /// Constructor.
    ///
    /// @param param Parameter* The parameter object
    /// @param kernel Kernel* The kernel used for calculations
    Compute(Parameter* param, Kernel* kernel);

    /// Destructor. Destroys the data fields properly, that were created
    /// during initialization.
    ~Compute();

    /// Calculates the current densities of the particles from their
    /// position and the kernel. Also returns the average density of
    /// the fluid.
    ///
    /// @return double The average density of the fluid
    double CalculateDensity();

    /// Calculates the pressure at the particle positions.
    void CalculatePressure();

    /// Calculates the forces that act on the particles. This consists of
    /// multiple factors like gravity, pressure and viscosity.
    void CalculateForces();

    /// Does the velocity integration for the leap-frog integration scheme.
    ///
    /// @param firstStep bool Flag if we are doing the first time step, which
    ///   requires some special handling
    void VelocityIntegration(bool firstStep);

    /// Does the position integration for the leap-frog integration scheme. If
    /// a particle would leave the domain boundaries, it will be reflected
    /// instead and a dampening force applied. In the default case the boundaries
    /// are the unit box that is infinitely extended towards the positive z axis.
    void PositionIntegration();

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
    /// @var _param Parameter* The parameter object containing the values
    /// of all necessary parameters.
    Parameter* _param;

    /// @var _kernel Kernel* The kernel to use for calculations.
    Kernel* _kernel;

    /// @var _isFirstStep bool Flag to indicate if we are doing the first time
    ///     step, which requires special handling.
    bool _isFirstStep;

    /// @var _initPotNrg double Holds the initial potential energy of the system,
    /// or at least an approximation
    double _initPotNrg;

    /// @var _vec1 double* A temporary 3D vector used in calculations.
    double* _vec1;

    /// @var _vec1 double* A temporary 3D vector used in calculations.
    double* _vec2;

    /// @var _matr1 double* A temporary 3x3 matrix used in calculations. The
    ///     matrix is indexed row by row.
    double* _matr1;

    /// @var _position double* The particle positions in x, y and z coordinates.
    double* _position;

    /// @var _velocity double* The particle velocities in x, y and z components.
    double* _velocity;

    /// @var _velocity double* The particle velocities in x, y and z components.
    double* _velocity_halfs;

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
