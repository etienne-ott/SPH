#include <cmath>
#include "compute.hpp"
#include "random_pool.hpp"

Compute::Compute(int N, Kernel* kernel) {
    _N = N;
    _kernel = kernel;

    _vec1 = new double[3];
    _vec2 = new double[3];
    _matr1 = new double[9];
    _position = new double[3 * N];
    _velocity = new double[3 * N];
    _force = new double[3 * N];
    _density = new double[N];
    _pressure = new double[N];

    RandomPool pool = RandomPool();

    for (int i = 0; i < N; i++) {
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
}

Compute::~Compute() {
    delete[] _position;
    delete[] _velocity;
    delete[] _force;
    delete[] _density;
    delete[] _pressure;
    delete[] _vec1;
    delete[] _vec2;
    delete[] _matr1;
}

void Compute::CalculateDensity() {
    double rho0 = 1000.0;
    double mass = rho0 / _N;

    for (int i = 0; i < _N; i++) {
        double sum = 0.0;
        double distance = 0.0;
        for (int j = 0; j < _N; j++) {
            distance = pow(
                (_position[j * 3] - _position[i * 3]) * (_position[j * 3] - _position[i * 3])
                    + (_position[j * 3 + 1] - _position[i * 3 + 1]) * (_position[j * 3 + 1] - _position[i * 3 + 1])
                    + (_position[j * 3 + 2] - _position[i * 3 + 2]) * (_position[j * 3 + 2] - _position[i * 3 + 2]),
                0.5
            );
            sum += mass * _kernel->Function(distance);
        }
        _density[i] = sum;
    }
}

void Compute::Timestep() {
    double k = 500.0;
    double g = 9.81;
    double dt = 0.1;
    double rho0 = 1000.0;
    double mass = rho0 / _N;
    double mu = 0.1;
    double dampening = 0.9;
    // @todo Force scaling should not be necessary
    double FSPressure = 0.00003;
    double FSGravity = 0.00001;
    double FSViscosity = 0.01;

    this->CalculateDensity();

    // Calculate pressure
    for (int i = 0; i < _N; i++) {
        _pressure[i] = k * (pow(_density[i] / rho0, 7.0) - 1);
    }

    // Calculate forces on particle i
    double distance = 0.0;
    double tmp = 0.0;

    for (int i = 0; i < _N; i++) {
        _force[i * 3] = 0.0;
        _force[i * 3 + 1] = 0.0;
        _force[i * 3 + 2] = 0.0;

        for (int j = 0; j < _N; j++) {
            if (j == i) {
                continue;
            }

            _vec1[0] = _position[i * 3] - _position[j * 3];
            _vec1[1] = _position[i * 3 + 1] - _position[j * 3 + 1];
            _vec1[2] = _position[i * 3 + 2] - _position[j * 3 + 2];
            distance = pow(_vec1[0] * _vec1[0] + _vec1[1] * _vec1[1] + _vec1[2] * _vec1[2], 0.5);

            // pressure
            _kernel->FOD(_vec1[0], _vec1[1], _vec1[2], distance, _vec2);
            tmp = 0.5 * (_pressure[i] + _pressure[j]) * (mass / _density[j]);

            _force[i * 3] -= FSPressure * tmp * _vec2[0];
            _force[i * 3 + 1] -= FSPressure * tmp * _vec2[1];
            _force[i * 3 + 2] -= FSPressure * tmp * _vec2[2];

            // viscosity
            _kernel->SOD(_vec1[0], _vec1[1], _vec1[2], distance, _matr1);
            tmp = FSViscosity * mu * mass / _density[j];
            _force[i * 3] -= tmp * (
                (_velocity[j * 3] - _velocity[i * 3]) * _matr1[0]
                + (_velocity[j * 3 + 1] - _velocity[i * 3 + 1]) * _matr1[1]
                + (_velocity[j * 3 + 2] - _velocity[i * 3 + 2]) * _matr1[2]
                );

            _force[i * 3 + 1] -= tmp * (
                (_velocity[j * 3] - _velocity[i * 3]) * _matr1[3]
                + (_velocity[j * 3 + 1] - _velocity[i * 3 + 1]) * _matr1[4]
                + (_velocity[j * 3 + 2] - _velocity[i * 3 + 2]) * _matr1[5]
                );

            _force[i * 3 + 2] -= tmp * (
                (_velocity[j * 3] - _velocity[i * 3]) * _matr1[6]
                + (_velocity[j * 3 + 1] - _velocity[i * 3 + 1]) * _matr1[7]
                + (_velocity[j * 3 + 2] - _velocity[i * 3 + 2]) * _matr1[8]
                );
        }

        // gravity
        _force[i * 3 + 2] -= _density[i] * mass * g * FSGravity;
    }

    // Do _velocity integration (explicit euler)
    for (int i = 0; i < _N; i++) {
        _velocity[i * 3] += dt * _force[i * 3] / mass;
        _velocity[i * 3 + 1] += dt * _force[i * 3 + 1] / mass;
        _velocity[i * 3 + 2] += dt * _force[i * 3 + 2] / mass;
    }

    // Do _position integration (explicit euler) with collision detection
    // against the standard rectangle spanned by (0,0,0)x(1,1,1), that is
    // open on the upper side (z > 1).
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
    for (int i = 0; i < _N; i++) {
        bool damp = false;

        // Reflection of x component
        newval = _position[i * 3] + _velocity[i * 3] * dt;
        if (newval < 0.0) {
            damp = true;
            _position[i * 3] = -newval;
            _velocity[i * 3] = -_velocity[i * 3];
        } else if (newval > 1.0) {
            damp = true;
            _position[i * 3] = 2.0 - newval;
            _velocity[i * 3] = -_velocity[i * 3];
        } else {
            _position[i * 3] = newval;
        }

        // Reflection of y component
        newval = _position[i * 3 + 1] + _velocity[i * 3 + 1] * dt;
        if (newval < 0.0) {
            damp = true;
            _position[i * 3 + 1] = -newval;
            _velocity[i * 3 + 1] = -_velocity[i * 3 + 1];
        } else if (newval > 1.0) {
            damp = true;
            _position[i * 3 + 1] = 2.0 - newval;
            _velocity[i * 3 + 1] = -_velocity[i * 3 + 1];
        } else {
            _position[i * 3 + 1] = newval;
        }

        // Reflection of z component
        newval = _position[i * 3 + 2] + _velocity[i * 3 + 2] * dt;
        if (newval < 0.0) {
            damp = true;
            _position[i * 3 + 2] = -newval;
            _velocity[i * 3 + 2] = -_velocity[i * 3 + 2];
        } else {
            _position[i * 3 + 2] = newval;
        }

        // Dampening
        if (damp) {
            _velocity[i * 3] *= dampening;
            _velocity[i * 3 + 1] *= dampening;
            _velocity[i * 3 + 2] *= dampening;
        }
    }
}

double* Compute::GetPosition() {
    return _position;
}

double* Compute::GetDensity() {
    return _density;
}