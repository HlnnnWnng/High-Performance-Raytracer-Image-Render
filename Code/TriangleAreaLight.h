#ifndef TRIANGLEAREALIGHT_H
#define TRIANGLEAREALIGHT_H

#include "Light.h"
#include "Vec3.h"
#include "Random.h"
#include <memory>

class TriangleAreaLight : public Light {
public:
    Vec3 v0, v1, v2;      // Triangle vertices
    Vec3 emission;        // Emission color (radiance)

    TriangleAreaLight(const Vec3& v0_, const Vec3& v1_, const Vec3& v2_, const Vec3& emission_);

    virtual Vec3 getPosition() const override;
    virtual Vec3 getIntensity() const override; // Not used for area lights

    // Get the area of the triangle
    double area() const;

    bool sample(const Vec3& referencePoint, Vec3& lightPoint, Vec3& lightNormal, double& pdf) const override {
        // Uniformly sample a point on the triangle using barycentric coordinates
        double u = Random::uniformDouble();
        double v = Random::uniformDouble();

        if (u + v > 1.0) {
            u = 1.0 - u;
            v = 1.0 - v;
        }

        lightPoint = v0 + u * (v1 - v0) + v * (v2 - v0);
        lightNormal = normal;
        pdf = 1.0 / area();
        return true;
    }

    Vec3 getRadiance() const override {
        return emission;
    }

private:
    Vec3 normal; // Triangle normal
};

#endif // TRIANGLEAREALIGHT_H
