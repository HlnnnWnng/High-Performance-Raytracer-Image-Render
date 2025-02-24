// BVHNode.h
#ifndef BVHNODE_H
#define BVHNODE_H

#include "Shape.h"    // Assuming Shape is the base class similar to Hittable
#include "AABB.h"
#include <vector>
#include <memory>
#include <algorithm> // For std::sort
#include <cstdlib>   // For rand()

class BVHNode : public Shape {
public:
    // Constructors
    BVHNode();
    BVHNode(std::vector<std::shared_ptr<Shape>>& objects, int start, int end);

    // Overridden Functions
    virtual bool intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const override;
    virtual bool boundingBox(AABB& outputBox) const override;
    std::string getType() const override{
        return "BVHNode";
    }

private:
    // Child Nodes
    std::shared_ptr<Shape> left;
    std::shared_ptr<Shape> right;

    // Bounding Box
    AABB box;

    // Comparator Function
    static int boxCompare(const std::shared_ptr<Shape>& a, const std::shared_ptr<Shape>& b, int axis);
};

#endif // BVHNODE_H
