#include <iostream>
#include <random>
#include <chrono>
#include <cmath>
#include "renderer.hpp"

double rand(std::default_random_engine& engine) {
    return (double)engine() / engine.max();
}

double kernel(double r, double h, int N) {
    double v = r / h;
    if (v >= 0.0 && v < 1.0) {
        v = (4.0 - 6.0 * r * r + 3.0 * r * r * r);
    } else if (v >= 1.0 && v < 2.0) {
        v = pow((2.0 - r), 3.0); 
    } else {
        return 0.0;
    }
    return v / (4.0 * h * h * h * N);
}

void printField(double* field, int len, int dim) {
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < dim; j++) {
            printf("%f ", field[i * dim + j]);
        }
        printf("\n");
    }
}

void initFields(int N, double* position, double* velocity, double* density) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    for (int i = 0; i < N; i++) {
        position[i * 3] = rand(generator);
        position[i * 3 + 1] = rand(generator);
        position[i * 3 + 2] = rand(generator);

        velocity[i * 3] = rand(generator) * 0.05 - 0.025;
        velocity[i * 3 + 1] = rand(generator) * 0.05 - 0.025;
        velocity[i * 3 + 2] = rand(generator) * 0.05 - 0.025;

        density[i] = 1.0 + rand(generator) * 0.2;
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
            iposy = (int)(R * pos[i * 3 + 1]),
            zrange = (int)(10 * pos[i * 3 + 2]);
        for (int k = iposx - zrange < 0 ? 0 : iposx - zrange; k < iposx + zrange && k < R; k++) {
            for (int l = iposy - zrange < 0 ? 0 : iposy - zrange; l < iposy + zrange && l < R; l++) {
                r->SetPixelRGB(k, l, 0, 0, 0);
            }
        }
    }
}

int main() {
    int N = 20;
    int R = 200;
    double h = 0.3;
    double t = 0.0;
    double tend = 50.0;
    double dt = 0.1;
    double mass = 1.0 / N;

    Renderer r = Renderer();
    r.Init(R, R);

    double* position = new double[3 * N];
    double* velocity = new double[3 * N];
    double* density = new double[N];

    initFields(N, position, velocity, density); 

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
                sum += mass * kernel(distance, h, N);
            }
            density[i] = sum;
        }

        // Do position integration (explicit euler)
        for (int i = 0; i < N; i++) {
            position[i * 3] += velocity[i * 3] * dt;
            position[i * 3 + 1] += velocity[i * 3 + 1] * dt;
            position[i * 3 + 2] += velocity[i * 3 + 2] * dt;
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
    delete[] density;

    return 0;
}