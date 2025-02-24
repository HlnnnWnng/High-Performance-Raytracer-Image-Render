// AABB.cpp
#include "AABB.h"

// Default Constructor
AABB::AABB() {}

// Parameterized Constructor
AABB::AABB(const Vec3& a, const Vec3& b) : _min(a), _max(b) {}

// Accessor for Minimum Point
Vec3 AABB::min() const { 
    return _min; 
}

// Accessor for Maximum Point
Vec3 AABB::max() const { 
    return _max; 
}

// Hit Function Implementation
bool AABB::hit(const Ray& ray, double tmin, double tmax) const {
    double ox = ray.origin.x;
    double oy = ray.origin.y;
    double oz = ray.origin.z;
    double dx = ray.direction.x;
    double dy = ray.direction.y;
    double dz = ray.direction.z;

    double tx_min, ty_min, tz_min;
    double tx_max, ty_max, tz_max;

    double x0 = _min.x, y0 = _min.y, z0 = _min.z;
    double x1 = _max.x, y1 = _max.y, z1 = _max.z;

    // X slab
    if (fabs(dx) < 1e-6) {
        if (ox < x0 || ox > x1)
            return false;
        tx_min = -INFINITY;
        tx_max = INFINITY;
    } else {
        if (dx >= 0) {
            tx_min = (x0 - ox) / dx;
            tx_max = (x1 - ox) / dx;
        } else {
            tx_min = (x1 - ox) / dx;
            tx_max = (x0 - ox) / dx;
        }
    }

    // Y slab
    if (fabs(dy) < 1e-6) {
        if (oy < y0 || oy > y1)
            return false;
        ty_min = -INFINITY;
        ty_max = INFINITY;
    } else {
        if (dy >= 0) {
            ty_min = (y0 - oy) / dy;
            ty_max = (y1 - oy) / dy;
        } else {
            ty_min = (y1 - oy) / dy;
            ty_max = (y0 - oy) / dy;
        }
    }

    // Z slab
    if (fabs(dz) < 1e-6) {
        if (oz < z0 || oz > z1)
            return false;
        tz_min = -INFINITY;
        tz_max = INFINITY;
    } else {
        if (dz >= 0) {
            tz_min = (z0 - oz) / dz;
            tz_max = (z1 - oz) / dz;
        } else {
            tz_min = (z1 - oz) / dz;
            tz_max = (z0 - oz) / dz;
        }
    }

    // Compute the largest t_min and smallest t_max
    double t0 = std::max({tx_min, ty_min, tz_min});
    double t1 = std::min({tx_max, ty_max, tz_max});

    return t0 < t1 && t1 > tmin && t0 < tmax;
}

// Static Function to Compute Surrounding Box
AABB AABB::surroundingBox(const AABB& box1, const AABB& box2) {
    Vec3 small(
        std::fmin(box1.min().x, box2.min().x),
        std::fmin(box1.min().y, box2.min().y),
        std::fmin(box1.min().z, box2.min().z)
    );
    Vec3 big(
        std::fmax(box1.max().x, box2.max().x),
        std::fmax(box1.max().y, box2.max().y),
        std::fmax(box1.max().z, box2.max().z)
    );
    return AABB(small, big);
}
