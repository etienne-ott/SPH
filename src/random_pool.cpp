#include "random_pool.h"
#include <random>
#include <chrono>

RandomPool::RandomPool()
    : RandomPool((long)std::chrono::system_clock::now().time_since_epoch().count())
    {}

RandomPool::RandomPool(long seed) {
    _seed = seed;
    _generator = std::default_random_engine(seed);
}

float RandomPool::NextFloat() {
    return (float)_generator() / _generator.max();
}

float RandomPool::NextFloat(float mean, float scale) {
    return mean + scale * ((float)_generator() / _generator.max() - 0.5);
}