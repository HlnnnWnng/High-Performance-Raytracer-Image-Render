// Shape.h
#ifndef SHAPE_H
#define SHAPE_H

#include "Ray.h"
#include "Vec3.h"
#include "Material.h"
#include "AABB.h"

struct HitRecord {
    Vec3 point;
    Vec3 normal;
    Material material;

    double t;              // Ray parameter at intersection
    bool front_face;      // Indicates if the intersection is on the front face
    double u, v;           // UV coordinates for texture mapping

    inline void set_face_normal(const Ray& r, const Vec3& outward_normal) {
        front_face = r.direction.dot(outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};



class Shape {
public:
    // Virtual Destructor for Proper Cleanup of Derived Classes
    virtual ~Shape() {}

    // Pure Virtual Functions
    virtual bool intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const = 0;
    virtual bool boundingBox(AABB& outputBox) const = 0;
    // Add this method to identify the shape type
    virtual std::string getType() const = 0;
};

#endif // SHAPE_H





