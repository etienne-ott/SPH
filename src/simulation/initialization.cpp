#include "util/random_pool.h"
#include "simulation/initialization.h"
#include "distribution/distributionEnums.h"
#include "distribution/domain.h"
#include "distribution/fastPoissonDisk.h"
#include "distribution/goldenSet.h"
#include "distribution/halton.h"
#include "distribution/hammersley.h"
#include "distribution/spherePacking.h"
#include "distribution/volumeGrid.h"
#include "distribution/whiteNoise.h"
#include <string>

Initialization::Initialization(YAML::Node& param) {
    _param = param;
}

int Initialization::InitPosition(float* position) {
    Domain dom = Domain(DomainType::Cube, _param);

    if (_param["domain_type"].as<std::string>() == "ellipsoid") {
        dom = Domain(DomainType::Sphere, _param);
    } else if (_param["domain_type"].as<std::string>() == "mesh") {
        Mesh m = Mesh();
        m.loadMeshFromOBJFile(_param["mesh_file"].as<std::string>());
        dom.setMesh(&m);
    }

    int nrCreated = 0;
    int N = _param["N"].as<int>();

    switch (_param["distribution_type"].as<int>()) {
        case 1: {
            VolumeGrid distr = VolumeGrid(N);
            nrCreated = distr.createPoints(position, &dom);
            break;
        }

        case 2: {
            SpherePacking distr = SpherePacking(N, _param["particle_size"].as<float>());
            nrCreated = distr.createPoints(position, &dom);
            break;
        }

        case 3: {
            WhiteNoise distr = WhiteNoise(_param["seed"].as<long>(), N);
            nrCreated = distr.createPoints(position, &dom);
            break;
        }

        case 4: {
            FastPoissonDisk distr = FastPoissonDisk(
                _param["seed"].as<long>(),
                N,
                _param["disk_radius"].as<float>(),
                _param["disk_tries"].as<int>()
            );
            nrCreated = distr.createPoints(position, &dom);
            break;
        }

        case 5: {
            Hammersley distr = Hammersley(
                N,
                _param["pb_1"].as<int>(),
                _param["pb_2"].as<int>()
            );
            nrCreated = distr.createPoints(position, &dom);
            break;
        }

        case 6: {
            Halton distr = Halton(
                N,
                _param["pb_1"].as<int>(),
                _param["pb_2"].as<int>(),
                _param["pb_3"].as<int>()
            );
            nrCreated = distr.createPoints(position, &dom);
            break;
        }
    }

    return nrCreated;
}

void Initialization::InitVelocity(float* velocity) {
    RandomPool pool = RandomPool();

    for (int i = 0; i < _param["N"].as<int>(); i++) {
        velocity[i * 3] = pool.nextFloat(0.0, 0.001);
        velocity[i * 3 + 1] = pool.nextFloat(0.0, 0.001);
        velocity[i * 3 + 2] = pool.nextFloat(0.0, 0.001);
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