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

void gradKernel(double rx, double ry, double rz, double r, double h, int N, double* ret) {
    double v = r / h, ir = 1.0 / r;

    if (v >= 0.0 && v < 1.0) {
        ret[0] = -12.0 * rx * ir + 9.0 * rx * pow(ir, 2.0);
        ret[1] = -12.0 * ry * ir + 9.0 * ry * pow(ir, 2.0);
        ret[2] = -12.0 * rz * ir + 9.0 * rz * pow(ir, 2.0);
    } else if (v >= 1.0 && v < 2.0) {
        ret[0] = -3.0 * pow(2.0 - r, 2.0) * rx * ir;
        ret[1] = -3.0 * pow(2.0 - r, 2.0) * ry * ir;
        ret[2] = -3.0 * pow(2.0 - r, 2.0) * rz * ir;
    } else {
        ret[0] = 0.0;
        ret[1] = 0.0;
        ret[2] = 0.0;
    }

    ret[0] /= (4.0 * h * h * h * N);
    ret[1] /= (4.0 * h * h * h * N);
    ret[2] /= (4.0 * h * h * h * N);
}

void writeDensity(double* density, double* position, double h, int N, double mass, int i) {
    char* filename = new char[50];
    sprintf(filename, "VTK/field_%i.vts", i);

    FILE* handle = fopen(filename, "w");
    delete filename;

    int size = 20;
    double ds = 1.0 / size;

    fprintf(handle, "<?xml version=\"1.0\"?>\n");
    fprintf(handle, "<VTKFile type=\"StructuredGrid\">\n");
    fprintf(handle, "<StructuredGrid WholeExtent=\"0 %i 0 %i 0 %i \">\n", size, size, size);
    fprintf(handle, "<Piece Extent=\"0 %i 0 %i 0 %i \">\n", size, size, size);
    fprintf(handle, "<Points>\n");
    fprintf(handle, "<DataArray type=\"Float64\" format=\"ascii\" NumberOfComponents=\"3\">\n");

    for (int z = 0; z <= size; ++z) {
        for (int y = 0; y <= size; ++y) {
            for (int x = 0; x <= size; ++x) {
                fprintf(handle, "%le %le %le\n",
                    (double)x * ds,
                    (double)y * ds,
                    (double)z * ds
                );
            }
        }
    }

    fprintf(handle, "</DataArray>\n");
    fprintf(handle, "</Points>\n");
    fprintf(handle, "<PointData>\n");

    fprintf(handle,
    "<DataArray Name=\"%s\" type=\"Float64\" format=\"ascii\">\n", "density");

    for (int z = 0; z <= size; ++z) {
        for (int y = 0; y <= size; ++y) {
            for (int x = 0; x <= size; ++x) {

                double sum = 0.0;
                double distance = 0.0;

                for (int j = 0; j < N; j++) {
                    distance = pow(
                        (position[j * 3] - x*ds) * (position[j * 3] - x*ds)
                            + (position[j * 3 + 1] - y*ds) * (position[j * 3 + 1] - y*ds)
                            + (position[j * 3 + 2] - z*ds) * (position[j * 3 + 2] - z*ds),
                        0.5
                    );
                    sum += mass * density[j] * kernel(distance, h, N);
                }

                fprintf(handle, "%le ", (sum > 0.5 ? 1.0 : 0.0));
            }

            fprintf(handle, "\n");
        }

        fprintf(handle, "\n");
    }

    fprintf(handle, "</DataArray>\n");
    fprintf(handle, "</PointData>\n");
    fprintf(handle, "</Piece>\n");
    fprintf(handle, "</StructuredGrid>\n");
    fprintf(handle, "</VTKFile>\n");

    fclose(handle);
}

void printField(double* field, int len, int dim) {
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < dim; j++) {
            printf("%f ", field[i * dim + j]);
        }
        printf("\n");
    }
}

void initFields(int N, double* position, double* velocity, double* force, double* density, double* pressure) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    for (int i = 0; i < N; i++) {
        position[i * 3] = rand(generator);
        position[i * 3 + 1] = rand(generator);
        position[i * 3 + 2] = rand(generator);

        velocity[i * 3] = rand(generator) * 0.05 - 0.025;
        velocity[i * 3 + 1] = rand(generator) * 0.05 - 0.025;
        velocity[i * 3 + 2] = rand(generator) * 0.05 - 0.025;

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
    double k = 500.0;
    double t = 0.0;
    double tend = 5.0;
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
                sum += mass * kernel(distance, h, N);
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

                gradKernel(vec1[0], vec1[1], vec1[2], distance, h, N, vec2);
                tmp = 0.5 * (pressure[i] + pressure[j]) * (mass / density[j]);

                force[i * 3] -= tmp * vec2[0];
                force[i * 3 + 1] -= tmp * vec2[1];
                force[i * 3 + 2] -= tmp * vec2[2];
            }
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

        writeDensity(density, position, h, N, mass, (int)(t / dt));

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