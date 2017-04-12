#include "parameter.hpp"

#ifndef __INITIALIZATION_HPP
#define __INITIALIZATION_HPP

class Initialization {
public:
    /// Constructor.
    ///
    /// @param param Parameter* The parameter object
    Initialization(Parameter* param);

    /// Initializes the position of the particles.
    ///
    /// @param position double* The position data
    void InitPosition(double* position);

    /// Initializes the velocity of the particles.
    ///
    /// @param velocity double* The velocity data
    void InitVelocity(double* velocity);

    /// Initializes the density of the particles.
    ///
    /// @param density double* The density data
    void InitDensity(double* density);

    /// Initializes the force data of the particles.
    ///
    /// @param force double* The force data
    void InitForce(double* force);

    /// Initializes the pressure of the particles.
    ///
    /// @param pressure double* The pressure data
    void InitPressure(double* pressure);

private:
    /// @var _param Parameter* The object containing the parameters
    Parameter* _param;
};
#endif // __INITIALIZATION_HPP
