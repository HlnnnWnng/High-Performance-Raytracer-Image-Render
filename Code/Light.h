// Light.h
#ifndef LIGHT_H
#define LIGHT_H

#include "Vec3.h"

class Light {
public:
    virtual ~Light() = default;
    virtual Vec3 getPosition() const = 0;
    virtual Vec3 getIntensity() const = 0;

        // Virtual function to sample a point on the light
    virtual bool sample(const Vec3& referencePoint, Vec3& lightPoint, Vec3& lightNormal, double& pdf) const = 0;

    // Virtual function to get radiance or intensity
    virtual Vec3 getRadiance() const = 0;

};

#endif // LIGHT_H
