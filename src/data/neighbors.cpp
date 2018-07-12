#include "data/neighbors.h"
#include <cmath>
#include <iostream>

Neighbors::Neighbors(float h, int N, float* bbox) {
    this->h = h;
    this->N = N;
    this->size_x = std::ceil((bbox[3] - bbox[0]) / h);
    this->size_y = std::ceil((bbox[4] - bbox[1]) / h);
    this->size_z = std::ceil((bbox[5] - bbox[2]) / h);

    this->grid = std::vector<std::vector<int>>();
    this->currentList = std::vector<int>();
    this->indices = new int[3 * N];

    for (int i = 0; i < size_x; i++) {
        for (int j = 0; j < size_y; j++) {
            for (int k = 0; k < size_z; k++) {
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

        this->grid.at(idx_z * size_y * size_x + idx_y * size_x + idx_x).push_back(i);
    }
}

std::vector<int> Neighbors::getNeighbors(int idx) {
    this->currentList.clear();

    for (int i = std::max(0, this->indices[idx * 3] - 1); i <= std::min(size_x - 1, this->indices[idx * 3] + 1); i++) {
        for (int j = std::max(0, this->indices[idx * 3 + 1] - 1); j <= std::min(size_y - 1, this->indices[idx * 3 + 1] + 1); j++) {
            for (int k = std::max(0, this->indices[idx * 3 + 2] - 1); k <= std::min(size_z - 1, this->indices[idx * 3 + 2] + 1); k++) {
                for (int l = 0; l < this->grid.at(k * size_y * size_x + j * size_x + i).size(); l++) {
                    this->currentList.push_back(this->grid.at(k * size_y * size_x + j * size_x + i).at(l));
                }
            }
        }
    }

    return this->currentList;
}