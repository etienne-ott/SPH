#include <random>

#ifndef __RANDOM_POOL_HPP
#define __RANDOM_POOL_HPP

class RandomPool {
public:
    /// Constructor. Initializes the random pool with the current unix
    /// time as seed.
    RandomPool();

    /// Constructor. Initializes the random pool with the given seed.
    RandomPool(long seed);

    /// Returns a randomly selected double value between 0.0 and 1.0
    ///
    /// @return double A random double between 0.0 and 1.0
    double NextDouble();

    /// Returns a randomly selected double value between
    /// (mean - 0.5 * scale) and (mean + 0.5 * scale).
    ///
    /// @param mean double The mean value of the distribution
    /// @param scale double The spread of the distribution
    /// @return double A random double in a range determined by the mean and
    ///     scale parameters.
    double NextDouble(double mean, double scale);

private:
    /// @var _seed long The seed used by the random number generator.
    long _seed;

    /// @var _generator std::default_random_engine The random engine used.
    std::default_random_engine _generator;
};
#endif //__RANDOM_POOL_HPP
