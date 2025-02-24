#include "TriangleAreaLight.h"
#include "Random.h"
#include <cmath>

TriangleAreaLight::TriangleAreaLight(const Vec3& v0_, const Vec3& v1_, const Vec3& v2_, const Vec3& emission_)
    : v0(v0_), v1(v1_), v2(v2_), emission(emission_) {
    normal = (v1 - v0).cross(v2 - v0).normalized();
}

Vec3 TriangleAreaLight::getPosition() const {
    // Return the centroid of the triangle
    return (v0 + v1 + v2) / 3.0;
}

Vec3 TriangleAreaLight::getIntensity() const {
    // Not meaningful for area lights; return zero or total power if needed
    return Vec3(0.0, 0.0, 0.0);
}

double TriangleAreaLight::area() const {
    return 0.5 * (v1 - v0).cross(v2 - v0).length();
}

