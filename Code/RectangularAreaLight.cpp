#include "RectangularAreaLight.h"
#include "Random.h"

RectangularAreaLight::RectangularAreaLight(const Vec3& position_, const Vec3& u_, const Vec3& v_, const Vec3& radiance_)
    : position(position_), u(u_), v(v_), radiance(radiance_) {
    normal = u.cross(v).normalized();
}

Vec3 RectangularAreaLight::getPosition() const {
    return position;
}

Vec3 RectangularAreaLight::getIntensity() const {
    // Not meaningful for area lights; can return zero or total power if needed
    return Vec3(0.0, 0.0, 0.0);
}

double RectangularAreaLight::area() const {
    return u.cross(v).length() * 4.0; // Since u and v are half-vectors
}

