#include <cmath>
#include <cstdio>
#include "compute.hpp"
#include "initialization.hpp"

Compute::Compute(YAML::Node& param, Kernel* kernel) {
    _param = param;
    _kernel = kernel;

    _isFirstStep = true;
    _initPotNrg = 0.0;
    int N = param["N"].as<int>();

    _vec1 = new float[3];
    _vec2 = new float[3];
    _matr1 = new float[9];
    _position = new float[3 * N];
    _velocity_halfs = new float[3 * N];
    _velocity = new float[3 * N];
    _force = new float[3 * N];
    _density = new float[N];
    _pressure = new float[N];

    Initialization init = Initialization(_param);
    init.InitPosition(_position);
    init.InitVelocity(_velocity);
    init.InitPressure(_pressure);
    init.InitForce(_force);
    init.InitDensity(_density);

    // @todo Make height for potential energy dependant of domain
    // @todo Somehow approximate initial potential energy of pressure
    // and viscosity, depends on starting conditions
    for (int i = 0; i < N; i++) {
        _initPotNrg += param["g"].as<float>() * param["mass"].as<float>() * _position[i * 3 + 2];
    }
}

Compute::~Compute() {
    delete[] _position;
    delete[] _velocity;
    delete[] _velocity_halfs;
    delete[] _force;
    delete[] _density;
    delete[] _pressure;
    delete[] _vec1;
    delete[] _vec2;
    delete[] _matr1;
}
void Compute::CalculateDensity() {
    int N = _param["N"].as<int>();
    float mass = _param["mass"].as<float>();

    for (int i = 0; i < N; i++) {
        float sum = 0.0;
        float distance = 0.0;

        for (int j = 0; j < N; j++) {
            distance = pow(
                (_position[i * 3] - _position[j * 3]) * (_position[i * 3] - _position[j * 3])
                    + (_position[i * 3 + 1] - _position[j * 3 + 1]) * (_position[i * 3 + 1] - _position[j * 3 + 1])
                    + (_position[i * 3 + 2] - _position[j * 3 + 2]) * (_position[i * 3 + 2] - _position[j * 3 + 2]),
                0.5
            );
            sum += mass * _kernel->ValueOf(distance);
        }

        _density[i] = sum;
    }
}

void Compute::Timestep() {
    this->CalculateDensity();
    this->CalculatePressure();
    this->CalculateForces();
    this->VelocityIntegration(_isFirstStep);
    this->PositionIntegration();

    _isFirstStep = false;
    printf("\n");
}

void Compute::CalculatePressure() {
    float rho0 = _param["rho0"].as<float>(),
        k = _param["k"].as<float>(),
        gamma = _param["gamma"].as<float>();

    for (int i = 0; i < _param["N"].as<int>(); i++) {
        _pressure[i] = k * rho0
            * (pow(_density[i] / rho0, gamma) - 1)
            / gamma;
    }
}

void Compute::CalculateForces() {
    float distance = 0.0;
    float tmp = 0.0;
    float dvx = 0.0, dvy = 0.0, dvz = 0.0;
    float kinNrg = 0.0;
    float mass = _param["mass"].as<float>();
    float g = _param["g"].as<float>();
    float epsilon = _param["epsilon"].as<float>();
    float h = _param["h"].as<float>();
    float mu = _param["mu"].as<float>();

    for (int i = 0; i < _param["N"].as<int>(); i++) {
        _force[i * 3] = 0.0;
        _force[i * 3 + 1] = 0.0;
        _force[i * 3 + 2] = 0.0;

        kinNrg += 0.5 * mass * (_velocity[i * 3] * _velocity[i * 3]
            + _velocity[i * 3 + 1] * _velocity[i * 3 + 1]
            + _velocity[i * 3 + 2] * _velocity[i * 3 + 2]);

        // 1-body forces, currently only gravity
        _force[i * 3 + 2] += mass * g;

        // 2-body forces
        // @todo Use symmetry to reduce number of calculations by a factor of 2
        float* pressureForce = new float[3];
        pressureForce[0] = 0.0;
        pressureForce[1] = 0.0;
        pressureForce[2] = 0.0;
        float* viscosityForce = new float[3];
        viscosityForce[0] = 0.0;
        viscosityForce[1] = 0.0;
        viscosityForce[2] = 0.0;

        for (int j = 0; j < _param["N"].as<int>(); j++) {
            if (j == i) {
                continue;
            }

            // Calculate distance vector
            _vec1[0] = _position[i * 3] - _position[j * 3];
            _vec1[1] = _position[i * 3 + 1] - _position[j * 3 + 1];
            _vec1[2] = _position[i * 3 + 2] - _position[j * 3 + 2];
            distance = pow(_vec1[0] * _vec1[0] + _vec1[1] * _vec1[1] + _vec1[2] * _vec1[2], 0.5);

            // Pressure force
            _kernel->FOD(_vec1[0], _vec1[1], _vec1[2], distance, _vec2);
            tmp = (
                _pressure[i] / (_density[i] * _density[i])
                + _pressure[j] / (_density[j] * _density[j])
            );

            pressureForce[0] += tmp * _vec2[0];
            pressureForce[1] += tmp * _vec2[1];
            pressureForce[2] += tmp * _vec2[2];

            // Viscosity force
            _kernel->SOD(_vec1[0], _vec1[1], _vec1[2], distance, _matr1);
            dvx = _velocity[j * 3] - _velocity[i * 3];
            dvy = _velocity[j * 3 + 1] - _velocity[i * 3 + 1];
            dvz = _velocity[j * 3 + 2] - _velocity[i * 3 + 2];
            tmp = (dvx * _vec1[0] + dvy * _vec1[1] + dvz * _vec1[2])
                / (distance * distance + epsilon * h * h);

            viscosityForce[0] += tmp * _vec1[0] * _matr1[0] / distance
                + tmp * _vec1[1] * _matr1[1] / distance
                + tmp * _vec1[2] * _matr1[2] / distance;

            viscosityForce[1] += tmp * _vec1[0] * _matr1[3] / distance
                + tmp * _vec1[1] * _matr1[4] / distance
                + tmp * _vec1[2] * _matr1[5] / distance;

            viscosityForce[2] += tmp * _vec1[0] * _matr1[6] / distance
                + tmp * _vec1[1] * _matr1[7] / distance
                + tmp * _vec1[2] * _matr1[8] / distance;
        }

        _force[i * 3] -= _density[i] * pressureForce[0];
        _force[i * 3 + 1] -= _density[i] * pressureForce[1];
        _force[i * 3 + 2] -= _density[i] * pressureForce[2];

        _force[i * 3] += mu * viscosityForce[0];
        _force[i * 3 + 1] += mu * viscosityForce[1];
        _force[i * 3 + 2] += mu * viscosityForce[2];

        delete pressureForce;
        delete viscosityForce;
    }

    printf("Kinetic energy: %f; Potential energy (apprx): %f; ", kinNrg, _initPotNrg - kinNrg);
}

void Compute::VelocityIntegration(bool firstStep) {
    if (firstStep) {
        for (int i = 0; i < _param["N"].as<int>(); i++) {
            _velocity_halfs[i * 3] = _velocity[i * 3] +  0.5 * _param["dt"].as<float>() * _force[i * 3] / _density[i];
            _velocity_halfs[i * 3 + 1] = _velocity[i * 3 + 1] +  0.5 * _param["dt"].as<float>() * _force[i * 3 + 1] / _density[i];
            _velocity_halfs[i * 3 + 2] = _velocity[i * 3 + 2] +  0.5 * _param["dt"].as<float>() * _force[i * 3 + 2] / _density[i];
        }
        return;
    }

    for (int i = 0; i < _param["N"].as<int>(); i++) {
        _velocity_halfs[i * 3] += _param["dt"].as<float>() * _force[i * 3] / _density[i];
        _velocity_halfs[i * 3 + 1] += _param["dt"].as<float>() * _force[i * 3 + 1] / _density[i];
        _velocity_halfs[i * 3 + 2] += _param["dt"].as<float>() * _force[i * 3 + 2] / _density[i];
        _velocity[i * 3] = _velocity_halfs[i * 3] + 0.5 * _param["dt"].as<float>() * _force[i * 3] / _density[i];
        _velocity[i * 3 + 1] = _velocity_halfs[i * 3 + 1] + 0.5 * _param["dt"].as<float>() * _force[i * 3 + 1] / _density[i];
        _velocity[i * 3 + 2] = _velocity_halfs[i * 3 + 2] + 0.5 * _param["dt"].as<float>() * _force[i * 3 + 2] / _density[i];
    }
}

void Compute::PositionIntegration() {
    // Do _position integration (leap frog) with collision detection
    // against the standard rectangle spanned by (0,0,0)x(1,1,1), that is
    // open on the upper side (z > 1) and extended infinitely.
    // On collision we reflect the respective component and rescale the
    // velocity so the new absolute velocity is the old one multiplied
    // with a dampening factor.
    // @todo We do the reflection with the old velocity, but the path after
    // the reflection should be traversed with the dampened velocity
    // @todo The dampening is uniform now, but should actually dampen the
    // reflected component stronger, while still maintaining the correct
    // direction
    // @todo in fact, the reflection is too crude, since the kernel seemingly
    // extends into the wall, but should be "squished" against it
    float newval = 0.0;
    for (int i = 0; i < _param["N"].as<int>(); i++) {
        bool damp = false;

        // Reflection of x component
        newval = _position[i * 3] + _velocity_halfs[i * 3] * _param["dt"].as<float>();
        if (newval < 0.0) {
            damp = true;
            _position[i * 3] = -newval;
            _velocity_halfs[i * 3] = -_velocity_halfs[i * 3];
        } else if (newval > 1.0) {
            damp = true;
            _position[i * 3] = 2.0 - newval;
            _velocity_halfs[i * 3] = -_velocity_halfs[i * 3];
        } else {
            _position[i * 3] = newval;
        }

        // Reflection of y component
        newval = _position[i * 3 + 1] + _velocity_halfs[i * 3 + 1] * _param["dt"].as<float>();
        if (newval < 0.0) {
            damp = true;
            _position[i * 3 + 1] = -newval;
            _velocity_halfs[i * 3 + 1] = -_velocity_halfs[i * 3 + 1];
        } else if (newval > 1.0) {
            damp = true;
            _position[i * 3 + 1] = 2.0 - newval;
            _velocity_halfs[i * 3 + 1] = -_velocity_halfs[i * 3 + 1];
        } else {
            _position[i * 3 + 1] = newval;
        }

        // Reflection of z component
        newval = _position[i * 3 + 2] + _velocity_halfs[i * 3 + 2] * _param["dt"].as<float>();
        if (newval < 0.0) {
            damp = true;
            _position[i * 3 + 2] = -newval;
            _velocity_halfs[i * 3 + 2] = -_velocity_halfs[i * 3 + 2];
        } else {
            _position[i * 3 + 2] = newval;
        }

        // Dampening
        if (damp) {
            _velocity_halfs[i * 3] *= _param["dampening"].as<float>();
            _velocity_halfs[i * 3 + 1] *= _param["dampening"].as<float>();
            _velocity_halfs[i * 3 + 2] *= _param["dampening"].as<float>();
        }
    }
}

float* Compute::GetPosition() {
    return _position;
}

float* Compute::GetDensity() {
    return _density;
}