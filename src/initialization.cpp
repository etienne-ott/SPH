#include "random_pool.h"
#include "initialization.h"

Initialization::Initialization(YAML::Node& param) {
    _param = param;
}

void Initialization::InitPosition(float* position) {
    RandomPool pool = RandomPool();

    for (int i = 0; i < _param["N"].as<int>(); i++) {
        position[i * 3] = pool.NextFloat(0.5f, 1.f);
        position[i * 3 + 1] = pool.NextFloat(0.5f, 1.f);
        position[i * 3 + 2] = pool.NextFloat(0.5f, 1.f);
    }
}

void Initialization::InitVelocity(float* velocity) {
    RandomPool pool = RandomPool();

    for (int i = 0; i < _param["N"].as<int>(); i++) {
        velocity[i * 3] = pool.NextFloat(0.0, 0.001);
        velocity[i * 3 + 1] = pool.NextFloat(0.0, 0.001);
        velocity[i * 3 + 2] = pool.NextFloat(0.0, 0.001);
    }
}

void Initialization::InitDensity(float* density) {
    float rho0 = _param["rho0"].as<float>();
    for (int i = 0; i < _param["N"].as<int>(); i++) {
        density[i] = rho0;
    }
}

void Initialization::InitForce(float* force) {
    for (int i = 0; i < _param["N"].as<int>(); i++) {
        force[i * 3] = 0.0;
        force[i * 3 + 1] = 0.0;
        force[i * 3 + 2] = 0.0;
    }
}

void Initialization::InitPressure(float* pressure) {
    for (int i = 0; i < _param["N"].as<int>(); i++) {
        pressure[i] = 0.0;
    }
}