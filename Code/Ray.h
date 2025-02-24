#ifndef RAY_H
#define RAY_H

#include "Vec3.h"

class Ray {
public:
    Vec3 origin;    // Ray origin
    Vec3 direction; // Ray direction (should be normalized)

    // Constructor
    Ray(const Vec3& origin_, const Vec3& direction_)
        : origin(origin_), direction(direction_.normalized()) {}

    // Compute a point along the ray at distance t
    Vec3 at(double t) const {
        return origin + direction * t;
    }
};

#endif // RAY_H
