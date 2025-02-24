// PointLight.cpp
#include "PointLight.h"

PointLight::PointLight(const Vec3& position_, const Vec3& intensity_)
    : position(position_), intensity(intensity_) {}

Vec3 PointLight::getPosition() const {
    return position;
}

Vec3 PointLight::getIntensity() const {
    return intensity;
}
