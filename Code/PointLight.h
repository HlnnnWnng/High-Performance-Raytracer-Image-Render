// PointLight.h
#ifndef POINTLIGHT_H
#define POINTLIGHT_H

#include "Light.h"

class PointLight : public Light {
public:
    Vec3 position;
    Vec3 intensity;

    PointLight(const Vec3& position_, const Vec3& intensity_);

    virtual Vec3 getPosition() const override;
    virtual Vec3 getIntensity() const override;

    // For PointLight, sampling always returns the fixed position and intensity
    bool sample(const Vec3& referencePoint, Vec3& lightPoint, Vec3& lightNormal, double& pdf) const override {
        lightPoint = position;
        lightNormal = Vec3(0.0, 0.0, 0.0); // Not used for point lights
        pdf = 1.0; // Delta distribution
        return true;
    }

    Vec3 getRadiance() const override {
        return intensity;
    }
};

#endif // POINTLIGHT_H
