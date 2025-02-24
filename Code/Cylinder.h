#ifndef CYLINDER_H
#define CYLINDER_H

#include "Shape.h"
#include "Vec3.h"
#include "Material.h"

class Cylinder : public Shape {
public:
    Vec3 center;
    Vec3 axis;
    double radius;
    double height;
    Material material;

    Cylinder(const Vec3& center_, const Vec3& axis_, double radius_, double height_, const Material& material_);

    virtual bool intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const override;
    virtual bool boundingBox(AABB& outputBox) const override;
    std::string getType() const override{
        return "Cylinder";
    }

private:
    Vec3 uAxis() const;
    Vec3 vAxis() const;
};

#endif // CYLINDER_H
