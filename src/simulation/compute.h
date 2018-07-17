#pragma once

#include "kernel/kernel.h"
#include "data/neighbors.h"
#include <yaml-cpp/yaml.h>

class Compute {
public:
    /// Constructor.
    ///
    /// @param param YAML::Node& The parameter object
    /// @param kernel_d Kernel* The kernel used for density calculation
    /// @param kernel_p Kernel* The kernel used for pressure calculation
    /// @param kernel_v Kernel* The kernel used for viscosity calculation
    Compute(YAML::Node& param, Kernel* kernel_d, Kernel* kernel_p, Kernel* kernel_v);

    /// Destructor. Destroys the data fields properly, that were created
    /// during initialization.
    ~Compute();

    /// Calculates the current densities of the particles from their
    /// position and the kernel. Also returns the average density of
    /// the fluid.
    ///
    /// @return float The average density of the fluid
    void CalculateDensity();

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
    /// for an overall number of 3*N floats.
    /// 
    /// @return float* The particle positions
    float* GetPosition();

    /// Returns the particle densities as consecutive values for each particle.
    /// 
    /// @return float* The particle densities
    float* GetDensity();

    /// Returns the particle pressure as consecutive values for each particle.
    /// 
    /// @return float* The particle pressure
    float* GetPressure();

private:
    /// @var _param YAML::Node The parameter object containing the values
    /// of all necessary parameters.
    YAML::Node _param;

    /// @var _kernel_density Kernel* The kernel to use for calculations.
    Kernel* _kernel_density;

    /// @var _kernel_pressure Kernel* The kernel to use for calculations.
    Kernel* _kernel_pressure;

    /// @var _kernel_viscosity Kernel* The kernel to use for calculations.
    Kernel* _kernel_viscosity;

    /// @var _neighbors Neighbors* A class used to get the neighbors of a particle
    Neighbors* _neighbors;

    /// @var _isFirstStep bool Flag to indicate if we are doing the first time
    ///     step, which requires special handling.
    bool _isFirstStep;

    /// @var _dr float* A temporary 3D vector used in calculations.
    float* _dr;

    /// @var _fod float* A temporary 3D vector used in calculations.
    float* _fod;

    /// @var _matr1 float* A temporary 3x3 matrix used in calculations. The
    ///     matrix is indexed row by row.
    float* _matr1;

    /// @var _position float* The particle positions in x, y and z coordinates.
    float* _position;

    /// @var _velocity float* The particle velocities in x, y and z components.
    float* _velocity;

    /// @var _velocity float* The particle velocities in x, y and z components.
    float* _velocity_halfs;

    /// @var _force float* The sum of all forces acting on the particles in x,
    ///     y and z components.
    float* _force;

    /// @var _density float* The density of the fluid particles.
    float* _density;

    /// @var _pressure float* The pressure of the fluid particles themselves,
    ///     as opposed to the pressure acting on them by other particles.
    float* _pressure;
};
