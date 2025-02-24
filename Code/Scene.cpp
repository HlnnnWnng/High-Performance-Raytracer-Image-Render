// Scene.cpp
#include "Scene.h"
#include <cfloat> // For DBL_MAX

// Constructor
Scene::Scene() : bvhRoot(nullptr), backgroundColor(0.0, 0.0, 0.0), ambientLight(0.5, 0.5, 0.5) {}

// Add a Shape to the Scene
void Scene::addShape(const std::shared_ptr<Shape>& shape) {
    shapes.push_back(shape);
}

// Add a Light to the Scene
void Scene::addLight(const std::shared_ptr<Light>& light) {
    lights.push_back(light);
}

// Build the BVH Tree
void Scene::buildBVH() {
    if (shapes.empty()) return;

    bvhRoot = std::make_shared<BVHNode>(shapes, 0, shapes.size());
}

// Intersection Method
bool Scene::intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const {
    if (bvhRoot) {
        return bvhRoot->intersect(ray, tMin, tMax, hitRecord);
    } else {
        // Fallback: Iterate Over All Shapes
        bool hitAnything = false;
        double closestSoFar = tMax;
        HitRecord tempRecord;

        for (const auto& shape : shapes) {
            if (shape->intersect(ray, tMin, closestSoFar, tempRecord)) {
                hitAnything = true;
                closestSoFar = tempRecord.t;
                hitRecord = tempRecord;
            }
        }
        return hitAnything;
    }
}

// Setter for Background Color
void Scene::setBackgroundColor(const Vec3& color) {
    backgroundColor = color;
}

// Getter for Background Color
Vec3 Scene::getBackgroundColor() const {
    return backgroundColor;
}

// Setter for Ambient Light
void Scene::setAmbientLight(const Vec3& ambient) {
    ambientLight = ambient;
}

// Getter for Ambient Light
Vec3 Scene::getAmbientLight() const {
    return ambientLight;
}

// Getter for Lights
const std::vector<std::shared_ptr<Light>>& Scene::getLights() const {
    return lights;
}


