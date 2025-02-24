#include "Triangle.h"
#include <algorithm> 

Triangle::Triangle(const Vec3& v0_, const Vec3& v1_, const Vec3& v2_,
                   const Vec2& uv0_, const Vec2& uv1_, const Vec2& uv2_,
                   const Material& material_)
    : v0(v0_), v1(v1_), v2(v2_),
      uv0(uv0_), uv1(uv1_), uv2(uv2_),
      material(material_) {
    // Precompute the face normal
    Vec3 faceNormal = (v1 - v0).cross(v2 - v0).normalized();
    normal0 = normal1 = normal2 = faceNormal;

    // Compute bounding box during construction
    double minX = std::min({v0.x, v1.x, v2.x});
    double minY = std::min({v0.y, v1.y, v2.y});
    double minZ = std::min({v0.z, v1.z, v2.z});
    double maxX = std::max({v0.x, v1.x, v2.x});
    double maxY = std::max({v0.y, v1.y, v2.y});
    double maxZ = std::max({v0.z, v1.z, v2.z});


    // Expand the box slightly to account for numerical errors
    const double epsilon = 1e-5;
    bbox = AABB(
        Vec3(minX - epsilon, minY - epsilon, minZ - epsilon),
        Vec3(maxX + epsilon, maxY + epsilon, maxZ + epsilon)
    );
}

Triangle::Triangle(const Vec3& v0_, const Vec3& v1_, const Vec3& v2_,
                 const Vec2& uv0_, const Vec2& uv1_, const Vec2& uv2_,
                 const Vec3& normal0_, const Vec3& normal1_, const Vec3& normal2_,
                 const Material& material_)
    : v0(v0_), v1(v1_), v2(v2_),
      uv0(uv0_), uv1(uv1_), uv2(uv2_),
      normal0(normal0_), normal1(normal1_), normal2(normal2_),
      material(material_) {
    // Compute bounding box during construction
    double minX = std::min({v0.x, v1.x, v2.x});
    double minY = std::min({v0.y, v1.y, v2.y});
    double minZ = std::min({v0.z, v1.z, v2.z});
    double maxX = std::max({v0.x, v1.x, v2.x});
    double maxY = std::max({v0.y, v1.y, v2.y});
    double maxZ = std::max({v0.z, v1.z, v2.z});

    // Expand the box slightly to account for numerical errors
    const double epsilon = 1e-5;
    bbox = AABB(
        Vec3(minX - epsilon, minY - epsilon, minZ - epsilon),
        Vec3(maxX + epsilon, maxY + epsilon, maxZ + epsilon)
    );
}

AABB Triangle::boundingMeshBox() const {
    return bbox;
}

bool Triangle::intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const {
    const double EPSILON = 1e-8;
    Vec3 edge1 = v1 - v0;
    Vec3 edge2 = v2 - v0;
    Vec3 h = ray.direction.cross(edge2);
    double a = edge1.dot(h);

    if (fabs(a) < EPSILON)
        return false; // Ray is parallel to the triangle

    double f = 1.0 / a;
    Vec3 s = ray.origin - v0;
    double u = f * s.dot(h);

    if (u < 0.0 || u > 1.0)
        return false;

    Vec3 q = s.cross(edge1);
    double v = f * ray.direction.dot(q);

    if (v < 0.0 || u + v > 1.0)
        return false;

    double t = f * edge2.dot(q);

    if (t < tMin || t > tMax)
        return false;

    // Compute barycentric coordinate w
    double w = 1.0 - u - v;

    // Intersection occurs
    hitRecord.t = t;
    hitRecord.point = ray.at(t);

    // **Interpolate the normal**
    Vec3 interpolatedNormal = (normal0 * w) + (normal1 * u) + (normal2 * v);
    interpolatedNormal = interpolatedNormal.normalized();

    // Ensure correct normal direction
    if (interpolatedNormal.dot(ray.direction) > 0) {
        interpolatedNormal = -interpolatedNormal;
    }

    hitRecord.normal = interpolatedNormal;
    hitRecord.material = material;

    // Interpolate UV coordinates
    double texU = w * uv0.x + u * uv1.x + v * uv2.x;
    double texV = w * uv0.y + u * uv1.y + v * uv2.y;

    hitRecord.u = texU;
    hitRecord.v = texV;

    return true;
}

bool Triangle::boundingBox(AABB& outputBox) const {
    Vec3 minPoint(
        std::min({ v0.x, v1.x, v2.x }),
        std::min({ v0.y, v1.y, v2.y }),
        std::min({ v0.z, v1.z, v2.z })
    );
    Vec3 maxPoint(
        std::max({ v0.x, v1.x, v2.x }),
        std::max({ v0.y, v1.y, v2.y }),
        std::max({ v0.z, v1.z, v2.z })
    );
    // Expand the bounding box slightly to account for floating-point inaccuracies
    double epsilon = 1e-4;
    minPoint = minPoint - Vec3(epsilon, epsilon, epsilon);
    maxPoint = maxPoint + Vec3(epsilon, epsilon, epsilon);
    outputBox = AABB(minPoint, maxPoint);
    return true;
}

