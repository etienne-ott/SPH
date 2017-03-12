#include <cmath>
#include "compute.hpp"
#include "random_pool.hpp"
#include "parameter.hpp"

Compute::Compute(Parameter* param, Kernel* kernel) {
    _param = param;
    _kernel = kernel;

    _isFirstStep = true;

    _vec1 = new double[3];
    _vec2 = new double[3];
    _matr1 = new double[9];
    _position = new double[3 * _param->N];
    _velocity_halfs = new double[3 * _param->N];
    _velocity = new double[3 * _param->N];
    _force = new double[3 * _param->N];
    _density = new double[_param->N];
    _pressure = new double[_param->N];

    RandomPool pool = RandomPool();

    for (int i = 0; i < _param->N; i++) {
        // We keep rerolling the position until we are within a radius of 0.4
        // of the position (0.5,0.5,1.0)
        // @todo Do this more cleverly with spherical coordinates
        double distance = 1e3;
        while (distance > 0.4) {
            _position[i * 3] = pool.NextDouble();
            _position[i * 3 + 1] = pool.NextDouble();
            _position[i * 3 + 2] = pool.NextDouble();
            distance = sqrt((0.5 - _position[i * 3]) * (0.5 - _position[i * 3])
                + (0.5 - _position[i * 3 + 1]) * (0.5 - _position[i * 3 + 1])
                + (1.0 - _position[i * 3 + 2]) * (1.0 - _position[i * 3 + 2]));
        }

        _velocity[i * 3] = pool.NextDouble(0.0, 0.1);
        _velocity[i * 3 + 1] = pool.NextDouble(0.0, 0.1);
        _velocity[i * 3 + 2] = pool.NextDouble(0.0, 0.1);

        _force[i * 3] = 0.0;
        _force[i * 3 + 1] = 0.0;
        _force[i * 3 + 2] = 0.0;

        _density[i] = 0.0;
        _pressure[i] = 0.0;
    }

    this->CalculateDensity();
    _param->NormalizeMass(_density, _param->N);
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

double Compute::CalculateDensity() {
    double avgsum = 0.0;

    for (int i = 0; i < _param->N; i++) {
        double sum = 0.0;
        double distance = 0.0;

        for (int j = 0; j < _param->N; j++) {
            distance = pow(
                (_position[i * 3] - _position[j * 3]) * (_position[i * 3] - _position[j * 3])
                    + (_position[i * 3 + 1] - _position[j * 3 + 1]) * (_position[i * 3 + 1] - _position[j * 3 + 1])
                    + (_position[i * 3 + 2] - _position[j * 3 + 2]) * (_position[i * 3 + 2] - _position[j * 3 + 2]),
                0.5
            );
            sum += _param->mass * _kernel->Function(distance);
        }

        _density[i] = sum;
        avgsum += sum;
    }

    return avgsum / _param->N;
}

void Compute::Timestep() {
    // Calculate density and update smoothing length. h should be set to
    // avg^(-1/d), where avg is the average density of the fluid and d is the
    // number of dimensions.
    double h = pow(this->CalculateDensity(), -1.0 / 3.0);
    _kernel->SetH(h);

    this->CalculatePressure();
    this->CalculateForces();
    this->VelocityIntegration(_isFirstStep);
    this->PositionIntegration();

    _isFirstStep = false;
}

void Compute::CalculatePressure() {
    for (int i = 0; i < _param->N; i++) {
        _pressure[i] = _param->k * (pow(_density[i] / _param->rho0, 7.0) - 1);
    }
}

void Compute::CalculateForces() {
    double distance = 0.0;
    double tmp = 0.0;

    for (int i = 0; i < _param->N; i++) {
        _force[i * 3] = 0.0;
        _force[i * 3 + 1] = 0.0;
        _force[i * 3 + 2] = 0.0;

        // 1-body forces, currently only gravity
        _force[i * 3 + 2] += _density[i] * _param->mass * _param->g * _param->FSGravity;

        // 2-body forces
        // @todo Use symmetry to reduce number of calculations by a factor of 2
        for (int j = 0; j < _param->N; j++) {
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
            tmp = 0.5 * (_pressure[i] + _pressure[j]) * (_param->mass / _density[j]);

            _force[i * 3] -= _param->FSPressure * tmp * _vec2[0];
            _force[i * 3 + 1] -= _param->FSPressure * tmp * _vec2[1];
            _force[i * 3 + 2] -= _param->FSPressure * tmp * _vec2[2];

            // Viscosity force
            _kernel->SOD(_vec1[0], _vec1[1], _vec1[2], distance, _matr1);
            tmp = _param->FSViscosity * _param->mu * _param->mass / _density[j];
            _force[i * 3] += tmp * (
                (_velocity[j * 3] - _velocity[i * 3]) * _matr1[0]
                + (_velocity[j * 3 + 1] - _velocity[i * 3 + 1]) * _matr1[1]
                + (_velocity[j * 3 + 2] - _velocity[i * 3 + 2]) * _matr1[2]
                );

            _force[i * 3 + 1] += tmp * (
                (_velocity[j * 3] - _velocity[i * 3]) * _matr1[3]
                + (_velocity[j * 3 + 1] - _velocity[i * 3 + 1]) * _matr1[4]
                + (_velocity[j * 3 + 2] - _velocity[i * 3 + 2]) * _matr1[5]
                );

            _force[i * 3 + 2] += tmp * (
                (_velocity[j * 3] - _velocity[i * 3]) * _matr1[6]
                + (_velocity[j * 3 + 1] - _velocity[i * 3 + 1]) * _matr1[7]
                + (_velocity[j * 3 + 2] - _velocity[i * 3 + 2]) * _matr1[8]
                );
        }
    }
}

void Compute::VelocityIntegration(bool firstStep) {
    if (firstStep) {
        for (int i = 0; i < _param->N; i++) {
            _velocity_halfs[i * 3] = _velocity[i * 3] +  0.5 + _param->dt * _force[i * 3] / _density[i];
            _velocity_halfs[i * 3 + 1] = _velocity[i * 3 + 1] +  0.5 + _param->dt * _force[i * 3 + 1] / _density[i];
            _velocity_halfs[i * 3 + 2] = _velocity[i * 3 + 2] +  0.5 + _param->dt * _force[i * 3 + 2] / _density[i];
        }
        return;
    }

    for (int i = 0; i < _param->N; i++) {
        _velocity_halfs[i * 3] += _param->dt * _force[i * 3] / _density[i];
        _velocity_halfs[i * 3 + 1] += _param->dt * _force[i * 3 + 1] / _density[i];
        _velocity_halfs[i * 3 + 2] += _param->dt * _force[i * 3 + 2] / _density[i];
        _velocity[i * 3] = _velocity_halfs[i * 3] + 0.5 * _param->dt * _force[i * 3] / _density[i];
        _velocity[i * 3 + 1] = _velocity_halfs[i * 3 + 1] + 0.5 * _param->dt * _force[i * 3 + 1] / _density[i];
        _velocity[i * 3 + 2] = _velocity_halfs[i * 3 + 2] + 0.5 * _param->dt * _force[i * 3 + 2] / _density[i];
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
    double newval = 0.0;
    for (int i = 0; i < _param->N; i++) {
        bool damp = false;

        // Reflection of x component
        newval = _position[i * 3] + _velocity_halfs[i * 3] * _param->dt;
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
        newval = _position[i * 3 + 1] + _velocity_halfs[i * 3 + 1] * _param->dt;
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
        newval = _position[i * 3 + 2] + _velocity_halfs[i * 3 + 2] * _param->dt;
        if (newval < 0.0) {
            damp = true;
            _position[i * 3 + 2] = -newval;
            _velocity_halfs[i * 3 + 2] = -_velocity_halfs[i * 3 + 2];
        } else {
            _position[i * 3 + 2] = newval;
        }

        // Dampening
        if (damp) {
            _velocity_halfs[i * 3] *= _param->dampening;
            _velocity_halfs[i * 3 + 1] *= _param->dampening;
            _velocity_halfs[i * 3 + 2] *= _param->dampening;
        }
    }
}

double* Compute::GetPosition() {
    return _position;
}

double* Compute::GetDensity() {
    return _density;
}