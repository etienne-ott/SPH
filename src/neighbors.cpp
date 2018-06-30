#include "neighbors.hpp"
#include <cmath>
#include <iostream>

Neighbors::Neighbors(float h, int N) {
    this->h = h;
    this->N = N;
    this->M = std::ceil(1.f / h);

    this->grid = std::vector<std::vector<int>>();
    this->currentList = std::vector<int>();
    this->indices = new int[3 * N];

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < M; k++) {
                this->grid.push_back(std::vector<int>());
            }
        }
    }
}

Neighbors::~Neighbors() {
    delete[] this->indices;
}

void Neighbors::sortParticlesIntoGrid(float* positions) {
    // clear old data
    for (int i = 0; i < this->grid.size(); i++) {
        this->grid.at(i).clear();
    }

    // sort particles into grid
    float invh = 1.f / h;

    for (int i = 0; i < N; i++) {
        int idx_x = std::floor(positions[i * 3] * invh);
        int idx_y = std::floor(positions[i * 3 + 1] * invh);
        int idx_z = std::floor(positions[i * 3 + 2] * invh);

        this->indices[i * 3] = idx_x;
        this->indices[i * 3 + 1] = idx_y;
        this->indices[i * 3 + 2] = idx_z;

        this->grid.at(idx_z * M * M + idx_y * M + idx_x).push_back(i);
    }
}

std::vector<int> Neighbors::getNeighbors(int idx) {
    this->currentList.clear();

    for (int i = std::max(0, this->indices[idx * 3] - 1); i <= std::min(M - 1, this->indices[idx * 3] + 1); i++) {
        for (int j = std::max(0, this->indices[idx * 3 + 1] - 1); j <= std::min(M - 1, this->indices[idx * 3 + 1] + 1); j++) {
            for (int k = std::max(0, this->indices[idx * 3 + 2] - 1); k <= std::min(M - 1, this->indices[idx * 3 + 2] + 1); k++) {
                for (int l = 0; l < this->grid.at(k * M * M + j * M + i).size(); l++) {
                    this->currentList.push_back(this->grid.at(k * M * M + j * M + i).at(l));
                }
            }
        }
    }

    return this->currentList;
}