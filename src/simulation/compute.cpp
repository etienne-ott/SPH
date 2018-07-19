#include <cmath>
#include <cstdio>
#include "simulation/compute.h"
#include "simulation/initialization.h"
#include "util/misc_math.h"
#include <string>

Compute::Compute(YAML::Node& param, Kernel* kernel_d, Kernel* kernel_p, Kernel* kernel_v) {
    _isFirstStep = true;
    int N = param["N"].as<int>();
    float h = param["h"].as<float>();

    _param = param;
    _kernel_density = kernel_d;
    _kernel_pressure = kernel_p;
    _kernel_viscosity = kernel_v;

    // We don't necessarily create all N particles,
    // so we need to reduce N to the actual number created
    Initialization init = Initialization(_param);
    float* tmp = new float[3 * N];

    int nrCreated = init.InitPosition(tmp);
    if (nrCreated < N) {
        printf("We wanted %d particles but only created %d\n", N, nrCreated);
    }
    _param["N"] = std::to_string(nrCreated);
    N = nrCreated;

    // copy over positions
    _position = new float[3 * N];
    for (int i = 0; i < 3 * N; i++) {
        _position[i] = tmp[i];
    }
    delete[] tmp;

    _dr = new float[3];
    _fod = new float[3];
    _matr1 = new float[9];
    _velocity_halfs = new float[3 * N];
    _velocity = new float[3 * N];
    _force = new float[3 * N];
    _density = new float[N];
    _pressure = new float[N];

    float* tmp2 = new float[6];
    tmp2[0] = param["bbox_x_lower"].as<float>();
    tmp2[1] = param["bbox_y_lower"].as<float>();
    tmp2[2] = param["bbox_z_lower"].as<float>();
    tmp2[3] = param["bbox_x_upper"].as<float>();
    tmp2[4] = param["bbox_y_upper"].as<float>();
    tmp2[5] = param["bbox_z_upper"].as<float>();

    _neighbors = new Neighbors(h, N, tmp2);
    delete[] tmp2;

    init.InitVelocity(_velocity);
    init.InitPressure(_pressure);
    init.InitForce(_force);
    init.InitDensity(_density);
}

Compute::~Compute() {
    delete[] _position;
    delete[] _velocity;
    delete[] _velocity_halfs;
    delete[] _force;
    delete[] _density;
    delete[] _pressure;
    delete[] _dr;
    delete[] _fod;
    delete[] _matr1;
    delete _neighbors;
}
void Compute::CalculateDensity() {
    int N = _param["N"].as<int>();
    float mass = _param["mass"].as<float>();
    float h = _param["h"].as<float>();

    for (int i = 0; i < N; i++) {
        float sum = 0.0;
        float distance = 0.0;

        std::vector<int> candidates = _neighbors->getNeighbors(i);
        for (uint k = 0; k < candidates.size(); k++) {
            int j = candidates.at(k);

            distance = fastSqrt2(
                (_position[i * 3] - _position[j * 3]) * (_position[i * 3] - _position[j * 3])
                    + (_position[i * 3 + 1] - _position[j * 3 + 1]) * (_position[i * 3 + 1] - _position[j * 3 + 1])
                    + (_position[i * 3 + 2] - _position[j * 3 + 2]) * (_position[i * 3 + 2] - _position[j * 3 + 2])
            );
            if (distance <= h) sum += mass * _kernel_density->ValueOf(distance);
        }

        _density[i] = sum;
    }
}

void Compute::Timestep() {
    _neighbors->sortParticlesIntoGrid(_position);

    this->CalculateDensity();
    this->CalculatePressure();
    this->CalculateForces();
    this->VelocityIntegration(_isFirstStep);
    this->PositionIntegration();

    _isFirstStep = false;
}

void Compute::CalculatePressure() {
    float rho0 = _param["rho0"].as<float>(),
        k = _param["k"].as<float>(),
        gamma = _param["gamma"].as<float>(),
        k_mod = k * rho0 / gamma;
    std::string model = _param["pressure_model"].as<std::string>();

    if (model == "P_GAMMA_ELASTIC") {
        for (int i = 0; i < _param["N"].as<int>(); i++) {
            _pressure[i] = (float)(k_mod * (pow(_density[i] / rho0, gamma) - 1.f));
            _pressure[i] = _pressure[i] * (_pressure[i] > 0);
        }
    } else if (model == "P_DIFFERENCE") {
        for (int i = 0; i < _param["N"].as<int>(); i++) {
            _pressure[i] = k * (_density[i] - rho0);
            _pressure[i] = _pressure[i] * (_pressure[i] > 0);
        }
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

    // reset forces
    for (int i = 0; i < _param["N"].as<int>(); i++) {
        _force[i * 3] = 0.0;
        _force[i * 3 + 1] = 0.0;
        _force[i * 3 + 2] = 0.0;
    }

    // now iterate over the particles and calculate the forces. we only do so
    // for other particles in the neighborhood with j > i and then apply the
    // forces to both particles in opposite directions. this uses the force
    // symmetry to improve performance
    int ix = 0, iy = 1, iz = 2;
    for (int i = 0; i < _param["N"].as<int>(); i++) {
        // calculate kinetic energy for debugging purposes
        kinNrg += 0.5 * mass * (_velocity[ix] * _velocity[ix]
            + _velocity[iy] * _velocity[iy]
            + _velocity[iz] * _velocity[iz]);

        // 1-body forces, currently only gravity
        _force[iy] += mass * g;

        // 2-body forces
        std::vector<int> candidates = _neighbors->getNeighbors(i);
        int jx, jy, jz;

        for (uint k = 0; k < candidates.size(); k++) {
            int j = candidates.at(k);

            if (j <= i) {
                continue;
            }

            jx = j * 3;
            jy = j * 3 + 1;
            jz = j * 3 + 2;

            // Calculate distance vector
            _dr[0] = _position[ix] - _position[jx];
            _dr[1] = _position[iy] - _position[jy];
            _dr[2] = _position[iz] - _position[jz];
            distance = fastSqrt2(_dr[0] * _dr[0] + _dr[1] * _dr[1] + _dr[2] * _dr[2]);

            // Pressure force
            _kernel_pressure->FOD(_dr[0], _dr[1], _dr[2], distance, _fod);
            tmp = mass * mass * (_pressure[i] / (_density[i] * _density[i])
                + _pressure[j] / (_density[j] * _density[j]));

            _force[ix] -= tmp * _fod[0];
            _force[iy] -= tmp * _fod[1];
            _force[iz] -= tmp * _fod[2];

            _force[jx] += tmp * _fod[0];
            _force[jy] += tmp * _fod[1];
            _force[jz] += tmp * _fod[2];

            // Viscosity force
            _kernel_viscosity->FOD(_dr[0], _dr[1], _dr[2], distance, _fod);
            dvx = _velocity[ix] - _velocity[jx];
            dvy = _velocity[iy] - _velocity[jy];
            dvz = _velocity[iz] - _velocity[jz];
            tmp = 2.f * mass * mass * mu / _density[j] / (
                _dr[0] * _dr[0] + _dr[1] * _dr[1] + _dr[2] * _dr[2]
                + epsilon * h * h
            );

            _force[ix] += tmp * dvx * (_dr[0] * _fod[0]);
            _force[iy] += tmp * dvy * (_dr[1] * _fod[1]);
            _force[iz] += tmp * dvz * (_dr[2] * _fod[2]);

            _force[jx] -= tmp * dvx * (_dr[0] * _fod[0]);
            _force[jy] -= tmp * dvy * (_dr[1] * _fod[1]);
            _force[jz] -= tmp * dvz * (_dr[2] * _fod[2]);
        }

        ix += 3; iy += 3; iz += 3;
    }

    printf("Kinetic energy: %f;", kinNrg);
}

void Compute::VelocityIntegration(bool firstStep) {
    float inv_mass = 1.f / _param["mass"].as<float>();
    float dt = _param["dt"].as<float>();
    float factor1, factor2;

    if (firstStep) {
        int ix = 0, iy = 1, iz = 2;
        factor1 = dt * inv_mass * 0.5;

        for (int i = 0; i < _param["N"].as<int>(); i++) {
            _velocity_halfs[ix] = _velocity[ix] +  factor1 * _force[ix];
            _velocity_halfs[iy] = _velocity[iy] +  factor1 * _force[iy];
            _velocity_halfs[iz] = _velocity[iz] +  factor1 * _force[iz];
            ix += 3; iy += 3; iz += 3;
        }
        return;
    }

    int ix = 0, iy = 1, iz = 2;
    factor1 = 0.5f * dt * inv_mass;
    factor2 = dt * inv_mass;

    for (int i = 0; i < _param["N"].as<int>(); i++) {
        _velocity_halfs[ix] += factor2 * _force[ix];
        _velocity_halfs[iy] += factor2 * _force[iy];
        _velocity_halfs[iz] += factor2 * _force[iz];

        _velocity[ix] = _velocity_halfs[ix] + factor1 * _force[ix];
        _velocity[iy] = _velocity_halfs[iy] + factor1 * _force[iy];
        _velocity[iz] = _velocity_halfs[iz] + factor1 * _force[iz];

        ix += 3; iy += 3; iz += 3;
    }
}

void Compute::PositionIntegration() {
    float dt = _param["dt"].as<float>();
    float lx = _param["bbox_x_lower"].as<float>();
    float ly = _param["bbox_y_lower"].as<float>();
    float lz = _param["bbox_z_lower"].as<float>();
    float ux = _param["bbox_x_upper"].as<float>();
    float uy = _param["bbox_y_upper"].as<float>();
    float uz = _param["bbox_z_upper"].as<float>();

    // Do _position integration (leap frog) with collision detection
    // against the standard cube spanned by (0,0,0)x(1,1,1).
    // On collision we reflect the respective component and rescale the
    // velocity so the new absolute velocity is the old one multiplied
    // with a dampening factor, except for the upper z-bound.
    // @todo We do the reflection with the old velocity, but the path after
    // the reflection should be traversed with the dampened velocity
    // @todo The dampening is uniform now, but should actually dampen the
    // reflected component stronger, while still maintaining the correct
    // direction
    // @todo in fact, the reflection is too crude, since the kernel seemingly
    // extends into the wall, but should be "squished" against it
    float newval = 0.0;
    float dampening = _param["dampening"].as<float>();
    int ix = 0, iy = 1, iz = 2;

    for (int i = 0; i < _param["N"].as<int>(); i++) {
        bool damp = false;

        // Reflection of x component
        newval = _position[ix] + _velocity_halfs[ix] * dt;
        if (newval < lx) {
            damp = true;
            _position[ix] = lx + fabs(lx - newval);
            _velocity_halfs[ix] = -_velocity_halfs[ix];
        } else if (newval > ux) {
            damp = true;
            _position[ix] = ux - fabs(newval - ux);
            _velocity_halfs[ix] = -_velocity_halfs[ix];
        } else {
            _position[ix] = newval;
        }

        // Reflection of y component
        newval = _position[iy] + _velocity_halfs[iy] * dt;
        if (newval < ly) {
            damp = true;
            _position[iy] = ly + fabs(ly - newval);
            _velocity_halfs[iy] = -_velocity_halfs[iy];
        } else if (newval > uy) {
            damp = true;
            _position[iy] = uy - fabs(newval - uy);
            _velocity_halfs[iy] = -_velocity_halfs[iy];
        } else {
            _position[iy] = newval;
        }

        // Reflection of z component
        newval = _position[iz] + _velocity_halfs[iz] * dt;
        if (newval < lz) {
            damp = true;
            _position[iz] = lz + fabs(lz - newval);
            _velocity_halfs[iz] = -_velocity_halfs[iz];
        } else if (newval > uz) {
            damp = true;
            _position[iz] = uz - fabs(newval - uz);
            _velocity_halfs[iz] = -_velocity_halfs[iz];
        } else {
            _position[iz] = newval;
        }

        // Dampening
        if (damp) {
            _velocity_halfs[ix] *= dampening;
            _velocity_halfs[iy] *= dampening;
            _velocity_halfs[iz] *= dampening;
        }

        ix += 3; iy += 3; iz += 3;
    }
}

float* Compute::GetPosition() {
    return _position;
}

float* Compute::GetDensity() {
    return _density;
}

float* Compute::GetPressure() {
    return _pressure;
}
