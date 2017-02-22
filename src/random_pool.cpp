#include "random_pool.hpp"
#include <random>
#include <chrono>

RandomPool::RandomPool()
    : RandomPool((long)std::chrono::system_clock::now().time_since_epoch().count())
    {}

RandomPool::RandomPool(long seed) {
    _seed = seed;
    _generator = std::default_random_engine(seed);
}

double RandomPool::NextDouble() {
    return (double)_generator() / _generator.max();
}

double RandomPool::NextDouble(double mean, double scale) {
    return mean + scale * ((double)_generator() / _generator.max() - 0.5);
}