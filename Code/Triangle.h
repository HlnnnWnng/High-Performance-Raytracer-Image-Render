#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Shape.h"
#include "Vec3.h"
#include "Vec2.h"
#include "Material.h"

class Triangle : public Shape {
public:
    Vec3 v0, v1, v2;
    Vec2 uv0, uv1, uv2;
    Vec3 normal0, normal1, normal2;
    Material material;
    AABB bbox;

    Triangle(const Vec3& v0_, const Vec3& v1_, const Vec3& v2_,
             const Vec2& uv0_, const Vec2& uv1_, const Vec2& uv2_,
             const Material& material_);

             // New Constructor with Normals
    Triangle(const Vec3& v0, const Vec3& v1, const Vec3& v2,
             const Vec2& uv0, const Vec2& uv1, const Vec2& uv2,
             const Vec3& normal0, const Vec3& normal1, const Vec3& normal2,
             const Material& material_);

    virtual bool intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const override;
    virtual bool boundingBox(AABB& outputBox) const override;
    
    std::string getType() const override{
        return "Triangle";
    }
    // Getter for bounding box (used in BVH)
    AABB boundingMeshBox() const;
};

#endif // TRIANGLE_H


