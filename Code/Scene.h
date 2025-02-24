// Scene.h
#ifndef SCENE_H
#define SCENE_H

#include "Shape.h"
#include "BVHNode.h"
#include "Light.h" // Assuming you have a Light class
#include <vector>
#include <memory>

class Scene {
public:
    // Constructors
    Scene();

    // Methods to Add Shapes and Lights
    void addShape(const std::shared_ptr<Shape>& shape);
    void addLight(const std::shared_ptr<Light>& light);

    // Method to Build BVH
    void buildBVH();

    // Intersection Method
    bool intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const;

    // Getter and Setter for Background Color
    void setBackgroundColor(const Vec3& color);
    Vec3 getBackgroundColor() const;

    // Getter and Setter for Ambient Light
    void setAmbientLight(const Vec3& ambient);
    Vec3 getAmbientLight() const;

    // Getter for Lights
    const std::vector<std::shared_ptr<Light>>& getLights() const;

private:
    std::vector<std::shared_ptr<Shape>> shapes;
    std::shared_ptr<BVHNode> bvhRoot;
    std::vector<std::shared_ptr<Light>> lights;

    // Background and Ambient Colors
    Vec3 backgroundColor;
    Vec3 ambientLight;
};

#endif // SCENE_H