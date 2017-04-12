#include "initialization.hpp"
#include "random_pool.hpp"

Initialization::Initialization(Parameter* param) {
    _param = param;
}

void Initialization::InitPosition(double* position) {
    double r = _param->h / _param->kappa, x = 0, y = 0, z = 0,
        offsetX = r - 0.5 * r, offsetY = r - 0.5 * r, offsetZ = r - 0.5 * r;
    bool oddRow = true, oddLayer = true;
    RandomPool pool = RandomPool();

    for (int i = 0; i < _param->N; i++) {
        position[i * 3] = offsetX + x + pool.NextDouble(0.0, 0.00001);
        position[i * 3 + 1] = offsetY + y + pool.NextDouble(0.0, 0.00001);
        position[i * 3 + 2] = offsetZ + z + pool.NextDouble(0.0, 0.00001);

        x += r;
        if (offsetX + x > 1.0 - r * 0.25) {
            oddRow = !oddRow;
            x = 0;
            offsetX = oddLayer ? (oddRow ? 0.5 * r : r) : (oddRow ? r : 0.5 * r);
            y += sqrt(3) * r * 0.5;
        }

        if (offsetY + y > 1.0 - r * 0.25) {
            oddLayer = !oddLayer;
            oddRow = true;
            x = 0;
            y = 0;
            offsetX = oddLayer ? 0.5 * r : r;
            offsetY = oddLayer ? 0.5 * r: (sqrt(3) - 1) * r * 0.5;
            z += sqrt(3) * r * 0.5;
        }
    }
}

void Initialization::InitVelocity(double* velocity) {
    RandomPool pool = RandomPool();

    for (int i = 0; i < _param->N; i++) {
        velocity[i * 3] = pool.NextDouble(0.0, 0.001);
        velocity[i * 3 + 1] = pool.NextDouble(0.0, 0.001);
        velocity[i * 3 + 2] = pool.NextDouble(0.0, 0.001);
    }
}

void Initialization::InitDensity(double* density) {
    for (int i = 0; i < _param->N; i++) {
        density[i] = 0.0;
    }
}

void Initialization::InitForce(double* force) {
    for (int i = 0; i < _param->N; i++) {
        force[i * 3] = 0.0;
        force[i * 3 + 1] = 0.0;
        force[i * 3 + 2] = 0.0;
    }
}

void Initialization::InitPressure(double* pressure) {
    for (int i = 0; i < _param->N; i++) {
        pressure[i] = 0.0;
    }
}