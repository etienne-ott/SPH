#include <cmath>
#include "compute.hpp"
#include "random_pool.hpp"

Compute::Compute(int N, Kernel* kernel) {
    _N = N;
    _kernel = kernel;

    _vec1 = new double[3];
    _vec2 = new double[3];
    _position = new double[3 * N];
    _velocity = new double[3 * N];
    _force = new double[3 * N];
    _density = new double[N];
    _pressure = new double[N];

    RandomPool pool = RandomPool();

    for (int i = 0; i < N; i++) {
        _position[i * 3] = pool.NextDouble();
        _position[i * 3 + 1] = pool.NextDouble();
        _position[i * 3 + 2] = pool.NextDouble();

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
}

void Compute::Timestep() {
    double k = 500.0;
    double g = 9.81;
    double dt = 0.1;
    double mass = 1.0 / _N;
    double forceScale = 0.000003;

    // Calculate density
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

    // Calculate pressure
    for (int i = 0; i < _N; i++) {
        _pressure[i] = k * (pow(_density[i] / 1.0, 7.0) - 1);
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

            _force[i * 3] -= tmp * _vec2[0];
            _force[i * 3 + 1] -= tmp * _vec2[1];
            _force[i * 3 + 2] -= tmp * _vec2[2];
        }

        // gravity
        _force[i * 3 + 2] -= _density[i] * mass * g;
    }

    // Do _velocity integration (explicit euler)
    for (int i = 0; i < _N; i++) {
        _velocity[i * 3] += forceScale * dt * _force[i * 3] / mass;
        _velocity[i * 3 + 1] += forceScale * dt * _force[i * 3 + 1] / mass;
        _velocity[i * 3 + 2] += forceScale * dt * _force[i * 3 + 2] / mass;
    }

    // Do _position integration (explicit euler)
    for (int i = 0; i < _N; i++) {
        _position[i * 3] += _velocity[i * 3] * dt;
        _position[i * 3 + 1] += _velocity[i * 3 + 1] * dt;
        _position[i * 3 + 2] += _velocity[i * 3 + 2] * dt;
    }
}

double* Compute::GetPosition() {
    return _position;
}

double* Compute::GetDensity() {
    return _density;
}