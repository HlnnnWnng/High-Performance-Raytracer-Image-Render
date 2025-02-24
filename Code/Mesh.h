#ifndef MESH_H
#define MESH_H

#include "Shape.h"
#include "Material.h"
#include "Triangle.h"

#include "Light.h"
#include "TriangleAreaLight.h"
#include "Vec3.h"
#include "Vec2.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include "MeshBVHNode.h"

struct FaceVertex {
    int v;   // Vertex index
    int vt;  // Texture coordinate index
    int vn;  // Normal index
};

struct Face {
    std::vector<FaceVertex> vertices;
    std::string materialName;
};

class Mesh : public Shape {
public:
    std::vector<Vec3> vertices;           // List of vertices
    std::vector<Vec2> texCoords;          // List of texture coordinates
    std::vector<Vec3> normals;            // List of normals
    std::vector<Face> faces;              // List of faces (with material names)
    std::unordered_map<std::string, Material> materialMap; // Materials loaded from MTL
    Material defaultMaterial;             // Default material

    // Constructor
    Mesh(const std::string& filename, const Material& defaultMaterial_, Scene& scene);

    // Retrieve area lights associated with this mesh
    const std::vector<std::shared_ptr<Light>>& getAreaLights() const { return areaLights; }

    // Intersection method
    virtual bool intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const override;
    virtual bool boundingBox(AABB& outputBox) const override;
    std::string getType() const override{
        return "Mesh";
    }

private:
    void loadMTLFile(const std::string& mtlFilename, const std::string& objFilePath);
    std::vector<std::shared_ptr<Triangle>> triangles;
    std::vector<std::shared_ptr<Light>> areaLights; // Store area lights here
    // Add a member to store the bounding box
    AABB boundingBox_;

    std::shared_ptr<MeshBVHNode> bvhRoot;

    void buildBVH();

    // Function to compute the bounding box during loading
    void computeBoundingBox();
};

#endif // MESH_H


