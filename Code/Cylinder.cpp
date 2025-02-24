// Cylinder.cpp
#include "Cylinder.h"
#include <cmath>
#include <limits>
#include <algorithm> // For std::min and std::max

Cylinder::Cylinder(const Vec3& center_, const Vec3& axis_, double radius_, double height_, const Material& material_)
    : center(center_), axis(axis_.normalized()), radius(radius_), height(height_), material(material_) {}

bool Cylinder::intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const {
    // Ray direction and origin
    Vec3 d = ray.direction;
    Vec3 o = ray.origin;

    // Cylinder's axis and center
    Vec3 ca = axis;
    Vec3 cc = center;

    // Vector from cylinder center to ray origin
    Vec3 oc = o - cc;

    // Components perpendicular to the cylinder's axis
    Vec3 d_ca = d - ca * d.dot(ca);
    Vec3 oc_ca = oc - ca * oc.dot(ca);

    double a = d_ca.dot(d_ca);
    double b = 2.0 * d_ca.dot(oc_ca);
    double c = oc_ca.dot(oc_ca) - radius * radius;

    double discriminant = b * b - 4.0 * a * c;

    if (discriminant < 0.0) {
        return false; // No intersection with infinite cylinder
    }

    double sqrtDiscriminant = std::sqrt(discriminant);

    // Find the nearest t value within [tMin, tMax] that intersects the finite cylinder
    double t0 = (-b - sqrtDiscriminant) / (2.0 * a);
    double t1 = (-b + sqrtDiscriminant) / (2.0 * a);

    double tNear = std::numeric_limits<double>::max();
    bool hit = false;
    bool capHit = false; // Flag to indicate if we hit a cap


    // Check both t0 and t1 for valid intersection within the cylinder's height
    for (double tCandidate : {t0, t1}) {
        if (tCandidate < tMin || tCandidate > tMax) {
            continue;
        }

        Vec3 point = ray.at(tCandidate);
        double y = (point - cc).dot(ca);
        // Height ranges from -height to +height (total height is 2 * height)
        if (y >= -height && y <= height) {
            if (tCandidate < tNear) {
                tNear = tCandidate;
                hit = true;
                capHit = false;

                // Populate HitRecord
                hitRecord.t = tNear;
                hitRecord.point = point;

                // Normal calculation
                Vec3 projection = cc + ca * y;
                Vec3 normal = (point - projection).normalized();

                hitRecord.normal = normal;
                hitRecord.material = material;

                // Compute UV coordinates for lateral surface
                // Compute u
                Vec3 dp = point - cc;
                Vec3 p_perp = dp - ca * y;
                double theta = atan2(p_perp.z, p_perp.x);
                if (theta < 0.0) {
                    theta += 2 * M_PI;
                }
                double u = theta / (2 * M_PI);

                // Compute v
                double v = (y + height) / (2 * height); // Map y from [-height, +height] to [0,1]
                v = std::clamp(v, 0.0, 1.0);

                // Assign UV coordinates
                hitRecord.u = u;
                hitRecord.v = v;
            }
        }
    }
    // Check for intersection with the caps (top and bottom disks)
    double denom = d.dot(ca);
    if (std::abs(denom) > 1e-6) { // Avoid division by zero
        // Top cap at y = +height
        double tCapTop = ((cc + ca * height - o).dot(ca)) / denom;
        if (tCapTop >= tMin && tCapTop <= tMax) {
            Vec3 p = ray.at(tCapTop);
            if ((p - (cc + ca * height)).lengthSquared() <= radius * radius) {
                if (tCapTop < tNear) {
                    tNear = tCapTop;
                    hit = true;
                    capHit = true;

                    // Populate HitRecord
                    hitRecord.t = tNear;
                    hitRecord.point = p;
                    hitRecord.normal = ca; // Normal points in the direction of the axis
                    hitRecord.material = material;

                    // Compute UV coordinates for top cap
                    Vec3 localP = p - (cc + ca * height);
                    double x = localP.dot(uAxis());
                    double z = localP.dot(vAxis());
                    double u = (x / (2 * radius)) + 0.5;
                    double v = (z / (2 * radius)) + 0.5;
                    hitRecord.u = u;
                    hitRecord.v = v;
                }
            }
        }

        // Bottom cap at y = -height
        double tCapBottom = ((cc - ca * height - o).dot(ca)) / denom;
        if (tCapBottom >= tMin && tCapBottom <= tMax) {
            Vec3 p = ray.at(tCapBottom);
            if ((p - (cc - ca * height)).lengthSquared() <= radius * radius) {
                if (tCapBottom < tNear) {
                    tNear = tCapBottom;
                    hit = true;
                    capHit = true;

                    // Populate HitRecord
                    hitRecord.t = tNear;
                    hitRecord.point = p;
                    hitRecord.normal = -ca; // Normal points opposite to the axis
                    hitRecord.material = material;

                    // Compute UV coordinates for bottom cap
                    Vec3 localP = p - (cc - ca * height);
                    double x = localP.dot(uAxis());
                    double z = localP.dot(vAxis());
                    double u = (x / (2 * radius)) + 0.5;
                    double v = (z / (2 * radius)) + 0.5;
                    hitRecord.u = u;
                    hitRecord.v = v;
                }
            }
        }
    }

    return hit;
}

// Helper functions to compute orthogonal axes to the cylinder's axis
Vec3 Cylinder::uAxis() const {
    Vec3 arbitrary = std::abs(axis.x) > 0.9 ? Vec3(0, 1, 0) : Vec3(1, 0, 0);
    return axis.cross(arbitrary).normalized();
}

Vec3 Cylinder::vAxis() const {
    return axis.cross(uAxis()).normalized();
}

bool Cylinder::boundingBox(AABB& outputBox) const {
    // Compute the bounding box of the cylinder
    Vec3 extent = axis * height;
    Vec3 minPoint = center - extent - Vec3(radius, radius, radius);
    Vec3 maxPoint = center + extent + Vec3(radius, radius, radius);

    // Component-wise min and max
    double minX = std::fmin(minPoint.x, maxPoint.x);
    double minY = std::fmin(minPoint.y, maxPoint.y);
    double minZ = std::fmin(minPoint.z, maxPoint.z);

    double maxX = std::fmax(minPoint.x, maxPoint.x);
    double maxY = std::fmax(minPoint.y, maxPoint.y);
    double maxZ = std::fmax(minPoint.z, maxPoint.z);

    outputBox = AABB(Vec3(minX, minY, minZ), Vec3(maxX, maxY, maxZ));
    return true;
}