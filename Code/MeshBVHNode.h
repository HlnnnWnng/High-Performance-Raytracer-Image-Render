#ifndef MESH_BVH_NODE_H
#define MESH_BVH_NODE_H

#include "AABB.h"
#include "Ray.h"
#include "Scene.h"
#include "Triangle.h"
#include <memory>
#include <vector>

class MeshBVHNode {
public:
    MeshBVHNode();
    MeshBVHNode(const std::vector<std::shared_ptr<Triangle>>& triangles, int start, int end);

    bool intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const;

private:
    std::shared_ptr<MeshBVHNode> left;
    std::shared_ptr<MeshBVHNode> right;
    std::shared_ptr<Triangle> triangle; // Non-nullptr for leaf nodes
    AABB box;
};

#endif // MESH_BVH_NODE_H
