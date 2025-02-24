#include "Sphere.h"

Sphere::Sphere(const Vec3& center_, double radius_, const Material& material_)
    : center(center_), radius(radius_), material(material_) {}

bool Sphere::intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const {
    Vec3 oc = ray.origin - center;
    double a = ray.direction.dot(ray.direction);
    double half_b = oc.dot(ray.direction);
    double c = oc.dot(oc) - radius * radius;
    double discriminant = half_b * half_b - a * c;

    if (discriminant > 0) {
        double sqrtD = std::sqrt(discriminant);
        double root = (-half_b - sqrtD) / a;
        if (root < tMax && root > tMin) {
            hitRecord.t = root;
            hitRecord.point = ray.at(hitRecord.t);
            hitRecord.normal = (hitRecord.point - center) / radius;
            hitRecord.material = material;

            // Compute UV coordinates
            Vec3 p_c = (hitRecord.point - center) / radius;
            double theta = atan2(p_c.z, p_c.x) + M_PI;
            double phi = acos(p_c.y);

            double u = theta / (2 * M_PI);
            double v = phi / M_PI;

            hitRecord.u = u;
            hitRecord.v = v;

            return true;
        }
        root = (-half_b + sqrtD) / a;
        if (root < tMax && root > tMin) {
            hitRecord.t = root;
            hitRecord.point = ray.at(hitRecord.t);
            hitRecord.normal = (hitRecord.point - center) / radius;
            hitRecord.material = material;

            // Compute UV coordinates
            Vec3 p_c = (hitRecord.point - center) / radius;
            double theta = atan2(p_c.z, p_c.x) + M_PI;
            double phi = acos(p_c.y);

            double u = theta / (2 * M_PI);
            double v = phi / M_PI;

            hitRecord.u = u;
            hitRecord.v = v;

            return true;
        }
    }
    return false;
}

bool Sphere::boundingBox(AABB& outputBox) const {
    Vec3 minPoint = center - Vec3(radius, radius, radius);
    Vec3 maxPoint = center + Vec3(radius, radius, radius);
    outputBox = AABB(minPoint, maxPoint);
    return true;
}

