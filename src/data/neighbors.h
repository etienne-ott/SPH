#pragma once

#include <vector>

class Neighbors {
private:
    std::vector<std::vector<int>> grid;
    int* indices;
    float h;
    int N;
    int size_x;
    int size_y;
    int size_z;
    float lx, ly, lz;
    std::vector<int> currentList;

public:
    Neighbors(float h, int N, float* bbox);

    ~Neighbors();

    void sortParticlesIntoGrid(float* positions);

    std::vector<int> getNeighbors(int idx);
};
