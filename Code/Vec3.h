#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

class Vec3 {
public:
    double x, y, z;

    // Constructors
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double val) : x(val), y(val), z(val) {} // Add this constructor
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    // Operator overloading
    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(double t) const { return Vec3(x * t, y * t, z * t); }
    Vec3 operator/(double t) const { return (*this) * (1 / t); }

    // Compound assignment operators
    Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    Vec3& operator-=(const Vec3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
    Vec3& operator*=(double t) { x *= t; y *= t; z *= t; return *this; }
    Vec3& operator/=(double t) { return *this *= 1 / t; }

    // Dot product
    double dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }

    // Cross product
    Vec3 cross(const Vec3& v) const { 
        return Vec3(
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        ); 
    }

    // Length and normalization
    double length() const { return std::sqrt(lengthSquared()); }
    double lengthSquared() const { return x * x + y * y + z * z; }
    Vec3 normalized() const { return *this / length(); }

    // Static utility functions
    static Vec3 random();
    static Vec3 random(double min, double max);

    // Output stream operator for debugging
    friend std::ostream& operator<<(std::ostream& os, const Vec3& v) {
        os << "[" << v.x << ", " << v.y << ", " << v.z << "]";
        return os;
    }

    static Vec3 hadamardProduct(const Vec3& a, const Vec3& b);

    // Overload the division operator for element-wise division
    Vec3 operator/(const Vec3& other) const;

    bool isZero() const {
        const double EPSILON = 1e-8;
        return std::abs(x) < EPSILON && std::abs(y) < EPSILON && std::abs(z) < EPSILON;
    }

};

// Utility functions for Vec3
Vec3 operator*(double t, const Vec3& v);

#endif // VEC3_H
