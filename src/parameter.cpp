#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "parameter.hpp"

using namespace std;

Parameter::Parameter() {
    N = 50;
    h = 0.3;
    Rix = 400;
    Riy = 600;
    Ro = 50;
    tend = 500.0;
    k = 500.0;
    gamma = 7.0;
    g = 9.81;
    dt = 0.1;
    rho0 = 1000.0;
    mass = 1.0;
    mu = 0.1;
    epsilon = 0.001;
    kappa = 1.0;
    dampening = 0.9;
    // @todo Force scaling should not be necessary
    FSPressure = 0.00003;
    FSGravity = 0.00001;
    FSViscosity = 0.01;
}

void Parameter::Load(const char *filename) {
    ifstream fin(filename);
    string line;

    while (getline(fin, line)) {
        if (
            line.find("#") != string::npos
            || line.length() == 0
        ) {
            continue;
        }

        istringstream sin(line.substr(line.find("=") + 1));

        // The following checks need to be in a specific order. The reason
        // for this is that looking for the string "h" might be true for
        // the input string "h = 0.1", but also for "rho0 = 1000". Therefore
        // we need to check the strings, that contain other, shorter parameter
        // names, first, then the short ones.
        // @todo Do this more cleverly, so the order doesn't matter and the
        // exact parameter name is required
        if (line.find("epsilon") != string::npos)
            sin >> epsilon;
        else if (line.find("kappa") != string::npos)
            sin >> kappa;
        else if (line.find("dampening") != string::npos)
            sin >> dampening;
        else if (line.find("fspressure") != string::npos)
            sin >> FSPressure;
        else if (line.find("fsgravity") != string::npos)
            sin >> FSGravity;
        else if (line.find("fsviscosity") != string::npos)
            sin >> FSViscosity;
        else if (line.find("gamma") != string::npos)
            sin >> gamma;
        else if (line.find("Rix") != string::npos)
            sin >> Rix;
        else if (line.find("Riy") != string::npos)
            sin >> Riy;
        else if (line.find("Ro") != string::npos)
            sin >> Ro;
        else if (line.find("rho0") != string::npos)
            sin >> rho0;
        else if (line.find("tend") != string::npos)
            sin >> tend;
        else if (line.find("dt") != string::npos)
            sin >> dt;
        else if (line.find("N") != string::npos)
            sin >> N;
        else if (line.find("h") != string::npos)
            sin >> h;
        else if (line.find("k") != string::npos)
            sin >> k;
        else if (line.find("g") != string::npos)
            sin >> g;
        else if (line.find("mu") != string::npos)
            sin >> mu;
        else {
            printf("Unknown parameter in line: %s\n", line.c_str());
        }
    }
}

double Parameter::NormalizeMass(double* density, int N) {
    double sumRho = 0.0, sumRhoSq = 0.0;
    for (int i = 0; i < N; i++) {
        sumRho += density[i];
        sumRhoSq += density[i] * density[i];
    }
    mass = rho0 * sumRho / sumRhoSq;
    return mass;
}