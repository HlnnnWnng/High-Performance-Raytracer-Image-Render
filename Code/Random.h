// Random.h
#ifndef RANDOM_H
#define RANDOM_H

#include <random>

class Random {
public:
    static std::mt19937& getGenerator() {
        static std::mt19937 generator(std::random_device{}());
        return generator;
    }

    static double uniformDouble() {
        static std::uniform_real_distribution<double> distribution(0.0, 1.0);
        return distribution(getGenerator());
    }

    static double uniformDoubleBetween(double min = 0.0, double max = 1.0) {
        static thread_local std::mt19937 generator(std::random_device{}());
        std::uniform_real_distribution<double> distribution(min, max);
        return distribution(generator);
    }
};

#endif // RANDOM_H
