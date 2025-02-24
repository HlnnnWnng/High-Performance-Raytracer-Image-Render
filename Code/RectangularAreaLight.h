#ifndef RECTANGULARAREALIGHT_H
#define RECTANGULARAREALIGHT_H

#include "Light.h"
#include "Vec3.h"
#include "Random.h"

class RectangularAreaLight : public Light {
public:
    Vec3 position; // Center of the rectangle
    Vec3 u;        // Vector along one side (half-width)
    Vec3 v;        // Vector along another side (half-height)
    Vec3 normal;   // Normal vector
    Vec3 radiance; // Emitted radiance (watts per steradian per square meter)

    RectangularAreaLight(const Vec3& position_, const Vec3& u_, const Vec3& v_, const Vec3& radiance_);

    virtual Vec3 getPosition() const override;
    virtual Vec3 getIntensity() const override; // Not used for area lights


    // Get the area of the rectangle
    double area() const;

    bool sample(const Vec3& referencePoint, Vec3& lightPoint, Vec3& lightNormal, double& pdf) const override {
        // Uniformly sample a point on the rectangle using two random numbers
        double r1 = Random::uniformDouble();
        double r2 = Random::uniformDouble();

        lightPoint = position + u * r1 + v * r2;
        lightNormal = normal;
        pdf = 1.0 / area();
        return true;
    }

    Vec3 getRadiance() const override {
        return radiance;
    }
};

#endif // RECTANGULARAREALIGHT_H
