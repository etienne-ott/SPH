#include "util/random_pool.h"
#include <random>
#include <chrono>

RandomPool::RandomPool()
    : RandomPool((long)std::chrono::system_clock::now().time_since_epoch().count())
    {}

RandomPool::RandomPool(long seed) {
    _seed = seed;
    _generator = std::default_random_engine(seed);

    // Draw and discard the first ten random numbers from the generator. See
    // the accompanying documentation for why this is necessary. Essentially,
    // these numbers are not equally distributed, which is what we want.
    for (int i = 0; i < 10; i++) {
        _generator();
    }
}

float RandomPool::nextFloat() {
    return (float)_generator() / _generator.max();
}

float RandomPool::nextFloat(float mean, float scale) {
    return mean + scale * ((float)_generator() / (float)_generator.max() - 0.5f);
}

double RandomPool::nextDouble() {
    return (double)_generator() / _generator.max();
}

double RandomPool::nextDouble(double mean, double scale) {
    return mean + scale * ((double)_generator() / (double)_generator.max() - 0.5);
}

int RandomPool::nextInt(int min, int max) {
    return (int) round(min - 0.5 + ((double)_generator() / _generator.max()) * (max - min + 1));
}