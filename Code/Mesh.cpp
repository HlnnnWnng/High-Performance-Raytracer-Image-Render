// Mesh.cpp
#include "Mesh.h"
#include "Triangle.h"
#include <filesystem>
#include <algorithm> // For std::min and std::max
#define USE_BVH

Mesh::Mesh(const std::string& filename, const Material& defaultMaterial_, Scene& scene)
    : defaultMaterial(defaultMaterial_) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open OBJ file " << filename << "\n";
        return;
    }

    // Get the directory of the OBJ file to resolve relative paths
    std::filesystem::path objFilePath(filename);
    std::string objDirectory = objFilePath.parent_path().string();

    std::string line;
    std::string currentMaterialName = ""; // Default material name
    while (getline(file, line)) {
        std::istringstream ss(line);
        std::string prefix;
        ss >> prefix;

        if (prefix == "v") {
            double x, y, z;
            ss >> x >> y >> z;
            vertices.emplace_back(x, y, z);
        } else if (prefix == "vt") {
            double u, v;
            ss >> u >> v;
            texCoords.emplace_back(u, v);
        } else if (prefix == "vn") {
            double x, y, z;
            ss >> x >> y >> z;
            normals.emplace_back(x, y, z);
        } else if (prefix == "f") {
            std::vector<FaceVertex> faceVertices;
            std::string vertexStr;
            while (ss >> vertexStr) {
                FaceVertex fv = { -1, -1, -1 };
                size_t firstSlash = vertexStr.find('/');
                size_t secondSlash = vertexStr.find('/', firstSlash + 1);

                if (firstSlash == std::string::npos) {
                    // Format: v
                    fv.v = std::stoi(vertexStr) - 1;
                } else if (secondSlash == std::string::npos) {
                    // Format: v/vt
                    fv.v = std::stoi(vertexStr.substr(0, firstSlash)) - 1;
                    fv.vt = std::stoi(vertexStr.substr(firstSlash + 1)) - 1;
                } else {
                    if (firstSlash + 1 == secondSlash) {
                        // Format: v//vn
                        fv.v = std::stoi(vertexStr.substr(0, firstSlash)) - 1;
                        fv.vn = std::stoi(vertexStr.substr(secondSlash + 1)) - 1;
                    } else {
                        // Format: v/vt/vn
                        fv.v = std::stoi(vertexStr.substr(0, firstSlash)) - 1;
                        fv.vt = std::stoi(vertexStr.substr(firstSlash + 1, secondSlash - firstSlash - 1)) - 1;
                        fv.vn = std::stoi(vertexStr.substr(secondSlash + 1)) - 1;
                    }
                }

                faceVertices.push_back(fv);
            }
            Face face;
            face.vertices = faceVertices;
            face.materialName = currentMaterialName; // Assign current material
            faces.push_back(face);
        } else if (prefix == "mtllib") {
            // Load material library
            std::string mtlFilename;
            ss >> mtlFilename;
            // Resolve relative path
            std::string fullPath = objDirectory + "/" + mtlFilename;
            loadMTLFile(fullPath, objDirectory);
        } else if (prefix == "usemtl") {
            // Set current material
            ss >> currentMaterialName;
        }
        // Ignore other prefixes
    }

    file.close();

    // Extract triangles and store them
    for (const auto& face : faces) {
        if (face.vertices.size() < 3)
            continue;

        // Retrieve material for this face
        Material faceMaterial = defaultMaterial;
        if (!face.materialName.empty()) {
            auto it = materialMap.find(face.materialName);
            if (it != materialMap.end()) {
                faceMaterial = it->second;
            }
        }

        // Triangulate the face (in case it's a polygon with more than 3 vertices)
        // Inside the loop where you create triangles
        for (size_t i = 1; i < face.vertices.size() - 1; ++i) {
            // Indices for the triangle
            const FaceVertex& fv0 = face.vertices[0];
            const FaceVertex& fv1 = face.vertices[i];
            const FaceVertex& fv2 = face.vertices[i + 1];

            // Vertex positions
            Vec3 vertex0 = vertices[fv0.v];
            Vec3 vertex1 = vertices[fv1.v];
            Vec3 vertex2 = vertices[fv2.v];

            // Texture coordinates
            Vec2 tex0(0.0, 0.0);
            Vec2 tex1(0.0, 0.0);
            Vec2 tex2(0.0, 0.0);

            if (fv0.vt >= 0 && fv0.vt < static_cast<int>(texCoords.size())) {
                tex0 = texCoords[fv0.vt];
            }
            if (fv1.vt >= 0 && fv1.vt < static_cast<int>(texCoords.size())) {
                tex1 = texCoords[fv1.vt];
            }
            if (fv2.vt >= 0 && fv2.vt < static_cast<int>(texCoords.size())) {
                tex2 = texCoords[fv2.vt];
            }

            // Compute face normal once
            Vec3 faceNormal = (vertex1 - vertex0).cross(vertex2 - vertex0).normalized();

            // Normals
            Vec3 normal0, normal1, normal2;

            if (fv0.vn >= 0 && fv0.vn < static_cast<int>(normals.size())) {
                normal0 = normals[fv0.vn];
            } else {
                normal0 = faceNormal;
            }

            if (fv1.vn >= 0 && fv1.vn < static_cast<int>(normals.size())) {
                normal1 = normals[fv1.vn];
            } else {
                normal1 = faceNormal;
            }

            if (fv2.vn >= 0 && fv2.vn < static_cast<int>(normals.size())) {
                normal2 = normals[fv2.vn];
            } else {
                normal2 = faceNormal;
            }

            // Create a triangle and store it
            auto triangle = std::make_shared<Triangle>(
                vertex0, vertex1, vertex2,
                tex0, tex1, tex2,
                normal0, normal1, normal2,
                faceMaterial
            );

            triangles.push_back(triangle);

            // If the material is emissive, add it as an area light
            if (faceMaterial.emission.lengthSquared() > 0.0) {
                auto areaLight = std::make_shared<TriangleAreaLight>(vertex0, vertex1, vertex2, faceMaterial.emission);
                areaLights.emplace_back(areaLight);
            }
        }
    }

    computeBoundingBox();

    #ifdef USE_BVH
    // Build the internal BVH
    buildBVH();
    #endif
}

void Mesh::loadMTLFile(const std::string& mtlFilename, const std::string& objFilePath) {
    std::ifstream mtlFile(mtlFilename);
    if (!mtlFile.is_open()) {
        std::cerr << "Error: Could not open MTL file " << mtlFilename << "\n";
        return;
    }

    std::string line;
    Material currentMaterial = defaultMaterial; // Start with default material
    std::string currentMaterialName;

    while (getline(mtlFile, line)) {
        std::istringstream ss(line);
        std::string prefix;
        ss >> prefix;

        if (prefix == "newmtl") {
            if (!currentMaterialName.empty()) {
                // Save the previous material
                materialMap[currentMaterialName] = currentMaterial;
            }
            ss >> currentMaterialName;
            currentMaterial = defaultMaterial; // Reset material to default
        } else if (prefix == "Kd") {
            double r, g, b;
            ss >> r >> g >> b;
            currentMaterial.diffuseColor = Vec3(r, g, b);
            currentMaterial.kd = (r + g + b) / 3.0;
        } else if (prefix == "Ks") {
            double r, g, b;
            ss >> r >> g >> b;
            currentMaterial.specularColor = Vec3(r, g, b);
        } else if (prefix == "Ns") {
            double ns;
            ss >> ns;
            currentMaterial.specularExponent = ns;
        } else if (prefix == "Ke") {
            double r, g, b;
            ss >> r >> g >> b;
            currentMaterial.emission = Vec3(r, g, b);
        } else if (prefix == "Tf") {
            double r, g, b;
            ss >> r >> g >> b;
            currentMaterial.transmissionFilter = Vec3(r, g, b);
        }else if (prefix == "Ni") {
            double ni;
            ss >> ni;
            currentMaterial.refractiveIndex = ni;
        } else if (prefix == "d") {
            double dissolveValue;
            ss >> dissolveValue;
            currentMaterial.transparency = 1- dissolveValue;
            if (dissolveValue < 1.0) {
                currentMaterial.kt = 1.0 - dissolveValue;
            } else {
                currentMaterial.kt = 0.0;
            }
        } else if (prefix == "illum") {
            int illumModel;
            ss >> illumModel;
        // Set flags based on illum model
            if (illumModel == 3 || illumModel == 5 || illumModel == 7 || illumModel == 8) {
                currentMaterial.isReflective = true;
            }
            if (illumModel == 4 || illumModel == 6 || illumModel == 7 || illumModel == 9) {
                currentMaterial.isRefractive = true;
            }

        } else if (prefix == "map_Kd") {
            std::string textureFilename;
            ss >> textureFilename;
            // Resolve relative path
            std::string texturePath = objFilePath + "/" + textureFilename;
            currentMaterial.texture = std::make_shared<Texture>(texturePath);
        }
        // Parse PBR parameters
        else if (prefix == "Pr") {
            double roughnessValue;
            ss >> roughnessValue;
            double roughnessClamped = std::max(roughnessValue, 0.001);
            currentMaterial.roughness = roughnessClamped;
            currentMaterial.specularExponent = pow((1.0 - roughnessClamped), 2) * 1000.0;
        } else if (prefix == "Pm") {
            double metallicValue;
            ss >> metallicValue;
            currentMaterial.metallic = metallicValue;
            currentMaterial.reflectivity = metallicValue;
        } else if (prefix == "sigma_a") {
            double sigma_a_value;
            ss >> sigma_a_value;
            currentMaterial.sigma_a = sigma_a_value;
        } else if (prefix == "sigma_s") {
            double sigma_s_value;
            ss >> sigma_s_value;
            currentMaterial.sigma_s = sigma_s_value;
        } else if (prefix == "g") {
            double g_value;
            ss >> g_value;
            currentMaterial.g = g_value;
        }

    }

    // Save the last material
    if (!currentMaterialName.empty()) {
        materialMap[currentMaterialName] = currentMaterial;
    }

    mtlFile.close();
}

bool Mesh::intersect(const Ray& ray, double tMin, double tMax, HitRecord& hitRecord) const {
    // First, check if the ray intersects the mesh's bounding box
    if (!boundingBox_.hit(ray, tMin, tMax))
        return false;

    // If the BVH is built, use it 
    if (bvhRoot) {
        return bvhRoot->intersect(ray, tMin, tMax, hitRecord);
    } else {
        // Fallback: Iterate through all triangles (should not reach here if BVH is built)
        bool hitAnything = false;
        double closestSoFar = tMax;

        for (const auto& triangle : triangles) {
            HitRecord tempRecord;
            if (triangle->intersect(ray, tMin, closestSoFar, tempRecord)) {
                hitAnything = true;
                closestSoFar = tempRecord.t;
                hitRecord = tempRecord;
            }
        }

        return hitAnything;
    }
}

void Mesh::computeBoundingBox() {
    if (vertices.empty()) {
        // Set an empty bounding box
        boundingBox_ = AABB();
        return;
    }

    // Initialize min and max points
    Vec3 minPoint(std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::max());
    Vec3 maxPoint(std::numeric_limits<double>::lowest(),
                  std::numeric_limits<double>::lowest(),
                  std::numeric_limits<double>::lowest());

    // Iterate over all vertices to find min and max
    for (const auto& vertex : vertices) {
        minPoint.x = std::min(minPoint.x, vertex.x);
        minPoint.y = std::min(minPoint.y, vertex.y);
        minPoint.z = std::min(minPoint.z, vertex.z);

        maxPoint.x = std::max(maxPoint.x, vertex.x);
        maxPoint.y = std::max(maxPoint.y, vertex.y);
        maxPoint.z = std::max(maxPoint.z, vertex.z);
    }

    boundingBox_ = AABB(minPoint, maxPoint);
}

bool Mesh::boundingBox(AABB& outputBox) const {
    outputBox = boundingBox_;
    return true;
}

void Mesh::buildBVH() {
    if (triangles.empty()) {
        bvhRoot = nullptr;
    } else {
        bvhRoot = std::make_shared<MeshBVHNode>(triangles, 0, triangles.size());
    }
}
