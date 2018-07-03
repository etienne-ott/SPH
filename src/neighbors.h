#pragma once

#include <vector>

class Neighbors {
private:
    std::vector<std::vector<int>> grid;
    int* indices;
    float h;
    int N;
    int M;
    std::vector<int> currentList;

public:
    Neighbors(float h, int N);

    ~Neighbors();

    void sortParticlesIntoGrid(float* positions);

    std::vector<int> getNeighbors(int idx);
};
