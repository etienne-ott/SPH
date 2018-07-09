#pragma once

#include <yaml-cpp/yaml.h>

class Initialization {
public:
    /// Constructor.
    ///
    /// @param param YAML::Node& The parameter object
    Initialization(YAML::Node& param);

    /// Initializes the position of the particles.
    ///
    /// @param position float* The position data
    /// @return int The number of created particles
    int InitPosition(float* position);

    /// Initializes the velocity of the particles.
    ///
    /// @param velocity float* The velocity data
    void InitVelocity(float* velocity);

    /// Initializes the density of the particles.
    ///
    /// @param density float* The density data
    void InitDensity(float* density);

    /// Initializes the force data of the particles.
    ///
    /// @param force float* The force data
    void InitForce(float* force);

    /// Initializes the pressure of the particles.
    ///
    /// @param pressure float* The pressure data
    void InitPressure(float* pressure);

private:
    /// @var _param YAML::Node The object containing the parameters
    YAML::Node _param;
};
