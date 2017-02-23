#include <iostream>
#include <cmath>
#include "renderer.hpp"
#include "kernel.hpp"
#include "random_pool.hpp"
#include "vtk.hpp"

#define WRITE_VTK_OUTPUT false

void initFields(int N, double* position, double* velocity, double* force, double* density, double* pressure) {
    RandomPool pool = RandomPool();

    for (int i = 0; i < N; i++) {
        position[i * 3] = pool.NextDouble();
        position[i * 3 + 1] = pool.NextDouble();
        position[i * 3 + 2] = pool.NextDouble();

        velocity[i * 3] = pool.NextDouble(0.0, 0.1);
        velocity[i * 3 + 1] = pool.NextDouble(0.0, 0.1);
        velocity[i * 3 + 2] = pool.NextDouble(0.0, 0.1);

        force[i * 3] = 0.0;
        force[i * 3 + 1] = 0.0;
        force[i * 3 + 2] = 0.0;

        density[i] = 0.0;
        pressure[i] = 0.0;
    }
}

void renderPositions(double* pos, Renderer* r, int N, int R) {
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < R; j++) {
            r->SetPixelRGB(i, j, 255, 255, 255);
        }
    }

    for (int i = 0; i < N; i++) {
        int iposx = (int)(R * pos[i * 3]),
            iposz = (int)(R * pos[i * 3 + 2]),
            yrange = (int)(10 * pos[i * 3 + 1]);
        for (int k = iposx - yrange < 0 ? 0 : iposx - yrange; k < iposx + yrange && k < R; k++) {
            for (int l = iposz - yrange < 0 ? 0 : iposz - yrange; l < iposz + yrange && l < R; l++) {
                r->SetPixelRGB(k, l, 0, 0, 0);
            }
        }
    }
}

int main() {
    int N = 50;
    int R = 400;
    double h = 0.3;
    double k = 500.0;
    double g = 9.81;
    double t = 0.0;
    double tend = 50.0;
    double dt = 0.1;
    double mass = 1.0 / N;
    double forceScale = 0.000003;

    Renderer r = Renderer();
    r.Init(R, R);

    double* vec1 = new double[3];
    double* vec2 = new double[3];
    double* position = new double[3 * N];
    double* velocity = new double[3 * N];
    double* force = new double[3 * N];
    double* density = new double[N];
    double* pressure = new double[N];

    Kernel kernel = Kernel(h, N);

    VTK vtk = VTK("VTK/", &kernel, mass, 20);

    initFields(N, position, velocity, force, density, pressure);

    renderPositions(position, &r, N, R);
    r.Render();

    while (t < tend) {
        // Calculate density
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            double distance = 0.0;
            for (int j = 0; j < N; j++) {
                distance = pow(
                    (position[j * 3] - position[i * 3]) * (position[j * 3] - position[i * 3])
                        + (position[j * 3 + 1] - position[i * 3 + 1]) * (position[j * 3 + 1] - position[i * 3 + 1])
                        + (position[j * 3 + 2] - position[i * 3 + 2]) * (position[j * 3 + 2] - position[i * 3 + 2]),
                    0.5
                );
                sum += mass * kernel.Function(distance);
            }
            density[i] = sum;
        }

        // Calculate pressure
        for (int i = 0; i < N; i++) {
            pressure[i] = k * (pow(density[i] / 1.0, 7.0) - 1);
        }

        // Calculate forces on particle i
        double distance = 0.0;
        double tmp = 0.0;

        for (int i = 0; i < N; i++) {
            force[i * 3] = 0.0;
            force[i * 3 + 1] = 0.0;
            force[i * 3 + 2] = 0.0;

            for (int j = 0; j < N; j++) {
                if (j == i) {
                    continue;
                }

                vec1[0] = position[i * 3] - position[j * 3];
                vec1[1] = position[i * 3 + 1] - position[j * 3 + 1];
                vec1[2] = position[i * 3 + 2] - position[j * 3 + 2];
                distance = pow(vec1[0] * vec1[0] + vec1[1] * vec1[1] + vec1[2] * vec1[2], 0.5);

                // pressure
                kernel.FOD(vec1[0], vec1[1], vec1[2], distance, vec2);
                tmp = 0.5 * (pressure[i] + pressure[j]) * (mass / density[j]);

                force[i * 3] -= tmp * vec2[0];
                force[i * 3 + 1] -= tmp * vec2[1];
                force[i * 3 + 2] -= tmp * vec2[2];
            }

            // gravity
            force[i * 3 + 2] -= density[i] * mass * g;
        }

        // Do velocity integration (explicit euler)
        for (int i = 0; i < N; i++) {
            velocity[i * 3] += forceScale * dt * force[i * 3] / mass;
            velocity[i * 3 + 1] += forceScale * dt * force[i * 3 + 1] / mass;
            velocity[i * 3 + 2] += forceScale * dt * force[i * 3 + 2] / mass;
        }

        // Do position integration (explicit euler)
        for (int i = 0; i < N; i++) {
            position[i * 3] += velocity[i * 3] * dt;
            position[i * 3 + 1] += velocity[i * 3 + 1] * dt;
            position[i * 3 + 2] += velocity[i * 3 + 2] * dt;
        }

        if (WRITE_VTK_OUTPUT) {
            vtk.WriteDensity(density, position, N);
        }

        renderPositions(position, &r, N, R);
        r.Render();
        SDL_Delay(30);

        t += dt;
    }

    bool running = true;
    SDL_Event event;

    while (running) {
        SDL_Delay(30);
        if (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            }
        }
    }

    delete[] position;
    delete[] velocity;
    delete[] force;
    delete[] density;
    delete[] pressure;
    delete[] vec1;
    delete[] vec2;

    return 0;
}