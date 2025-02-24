#include "MeshBVHNode.h"
#include <algorithm>
#include <cstdlib> // For rand()

MeshBVHNode::MeshBVHNode() : left(nullptr), right(nullptr), triangle(nullptr) {}

MeshBVHNode::MeshBVHNode(const std::vector<std::shared_ptr<Triangle>>& triangles, int start, int end) {
    int objectSpan = end - start;

    if (objectSpan == 1) {
        triangle = triangles[start];
        box = triangle->boundingMeshBox();
        left = right = nullptr;
    } else {
        // Choose axis to sort along: 0 = x, 1 = y, 2 = z
        int axis = rand() % 3;

        auto comparator = [axis](const std::shared_ptr<Triangle>& a, const std::shared_ptr<Triangle>& b) {
            double aMin = 0.0;
            double bMin = 0.0;

            if (axis == 0) {
                aMin = a->boundingMeshBox().min().x;
                bMin = b->boundingMeshBox().min().x;
            } else if (axis == 1) {
                aMin = a->boundingMeshBox().min().y;
                bMin = b->boundingMeshBox().min().y;
            } else if (axis == 2) {
                aMin = a->boundingMeshBox().min().z;
                bMin = b->boundingMeshBox().min().z;
            }

            return aMin < bMin;
        };

        // Sort triangles based on the chosen axis
        std::vector<std::shared_ptr<Triangle>> sortedTriangles(triangles.begin() + start, triangles.begin() + end);
        std::sort(sortedTriangles.begin(), sortedTriangles.end(), comparator);

        int mid = objectSpan / 2;
        left = std::make_shared<MeshBVHNode>(sortedTriangles, 0, mid);
        right = std::make_shared<MeshBVHNode>(sortedTriangles, mid, objectSpan);

        box = AABB::surroundingBox(left->box, right->box);
        triangle = nullptr;
    }
}

bool MeshBVHNode::intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const {
    if (!box.hit(ray, tMin, tMax))
        return false;

    if (triangle) {
        return triangle->intersect(ray, tMin, tMax, hitRecord);
    }

    bool hitLeft = left && left->intersect(ray, tMin, tMax, hitRecord);
    bool hitRight = right && right->intersect(ray, tMin, hitLeft ? hitRecord.t : tMax, hitRecord);

    return hitLeft || hitRight;
}
