// AABB.h
#ifndef AABB_H
#define AABB_H

#include "Vec3.h"
#include "Ray.h"
#include <algorithm> // For std::fmin and std::fmax
#include <cmath>     // For fabs

class AABB {
public:
    // Constructors
    AABB();
    AABB(const Vec3& a, const Vec3& b);

    // Accessors
    Vec3 min() const;
    Vec3 max() const;

    // Member Functions
    bool hit(const Ray& ray, double tMin, double tMax) const;

    // Static Utility Function
    static AABB surroundingBox(const AABB& box1, const AABB& box2);

private:
    Vec3 _min;
    Vec3 _max;
};

#endif // AABB_H
