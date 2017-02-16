#include <iostream>
#include <random>
#include <chrono>
#include <cmath>

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

int main() {
    int N = 20;
    double h = 0.3;
    double mass = 1.0 / N;

    double* position = new double[3 * N];
    double* velocity = new double[3 * N];
    double* density = new double[N];

    initFields(N, position, velocity, density); 

    printField(velocity, N, 3);
    printf("\n");

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

    printField(density, N, 1);

    printf("All finished\n");

    delete[] position;
    delete[] velocity;
    delete[] density;

    return 0;
}