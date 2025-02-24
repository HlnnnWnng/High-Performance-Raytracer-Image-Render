#ifndef SPHERE_H
#define SPHERE_H

#include "Shape.h"
#include "Vec3.h"
#include "Material.h"

class Sphere : public Shape {
public:
    Vec3 center;
    double radius;
    Material material;

    Sphere(const Vec3& center_, double radius_, const Material& material_);

    virtual bool intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const override;
    virtual bool boundingBox(AABB& outputBox) const override;
    std::string getType() const override{
        return "Sphere";
    }
};

#endif // SPHERE_H
