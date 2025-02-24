#include "Vec3.h"
#include <random>

// Generate a random Vec3 with components in [0, 1)
Vec3 Vec3::random() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return Vec3(distribution(generator), distribution(generator), distribution(generator));
}

// Generate a random Vec3 with components in [min, max)
Vec3 Vec3::random(double min, double max) {
    static std::uniform_real_distribution<double> distribution(min, max);
    static std::mt19937 generator;
    return Vec3(distribution(generator), distribution(generator), distribution(generator));
}

// Scalar multiplication from the left
Vec3 operator*(double t, const Vec3& v) {
    return Vec3(t * v.x, t * v.y, t * v.z);
}

Vec3 Vec3::hadamardProduct(const Vec3& a, const Vec3& b) {
    return Vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}

// Element-wise division
Vec3 Vec3::operator/(const Vec3& other) const {
    return Vec3(x / other.x, y / other.y, z / other.z);
}

