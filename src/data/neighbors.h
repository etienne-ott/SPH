#pragma once

#include <vector>
#include "util/parallel_bounds.h"

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

public:
    Neighbors(float h, int N, float* bbox);

    ~Neighbors();

    void sortParticlesIntoGrid(float* positions, ParallelBounds& pBounds);

    void getNeighbors(int idx, std::vector<int>& list);
};
