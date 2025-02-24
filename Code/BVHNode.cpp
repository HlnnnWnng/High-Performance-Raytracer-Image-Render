// BVHNode.cpp
#include "BVHNode.h"

// Default Constructor
BVHNode::BVHNode() : left(nullptr), right(nullptr), box() {}

// Comparator Function Implementation
int BVHNode::boxCompare(const std::shared_ptr<Shape>& a, const std::shared_ptr<Shape>& b, int axis) {
    AABB boxA, boxB;

    // Compare along the specified axis
    if (axis == 0) { // X-axis
        if (boxA.min().x < boxB.min().x)
            return -1;
        else
            return 1;
    } else if (axis == 1) { // Y-axis
        if (boxA.min().y < boxB.min().y)
            return -1;
        else
            return 1;
    } else { // Z-axis
        if (boxA.min().z < boxB.min().z)
            return -1;
        else
            return 1;
    }
}


// Parameterized Constructor
BVHNode::BVHNode(std::vector<std::shared_ptr<Shape>>& objects, int start, int end) {
    int objectSpan = end - start;

    if (objectSpan == 1) {
        left = right = objects[start];
    } else if (objectSpan == 2) {
        if (boxCompare(objects[start], objects[start + 1], 0) < 0) { // 0 for x-axis
            left = objects[start];
            right = objects[start + 1];
        } else {

            left = objects[start + 1];
            right = objects[start];
        }
    } else {
        // Choose axis to sort along: 0 = x, 1 = y, 2 = z
        int axis = rand() % 3;

        // Define comparator based on the chosen axis
        auto comparator = [&](const std::shared_ptr<Shape>& a, const std::shared_ptr<Shape>& b) -> bool {
            return boxCompare(a, b, axis) < 0;
        };

        // Sort objects based on the chosen axis
        std::sort(objects.begin() + start, objects.begin() + end, comparator);

        // Recursively build BVH
        int mid = start + objectSpan / 2;
        left = std::make_shared<BVHNode>(objects, start, mid);
        right = std::make_shared<BVHNode>(objects, mid, end);
    }

    // Compute the bounding box for this node
    AABB boxLeft, boxRight;
    if (!left->boundingBox(boxLeft) || !right->boundingBox(boxRight)) {
        // Handle error appropriately, e.g., throw an exception
        // For simplicity, we'll assume all shapes have a bounding box
    }

    box = AABB::surroundingBox(boxLeft, boxRight);
}

// Overridden intersect Function Implementation
bool BVHNode::intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const {
    if (!box.hit(ray, tMin, tMax))
        return false;
    HitRecord leftRecord, rightRecord;
    bool hitLeft = left->intersect(ray, tMin, tMax, leftRecord);
    bool hitRight = right->intersect(ray, tMin, hitLeft ? leftRecord.t : tMax, rightRecord);

    if (hitLeft && hitRight) {
        if (leftRecord.t < rightRecord.t)
            hitRecord = leftRecord;
        else
            hitRecord = rightRecord;
        return true;
    } else if (hitLeft) {
        hitRecord = leftRecord;
        return true;
    } else if (hitRight) {
        hitRecord = rightRecord;
        return true;
    } else {
        return false;
    }
}

// Overridden boundingBox Function Implementation
bool BVHNode::boundingBox(AABB& outputBox) const {
    outputBox = box;
    return true;
}



