// raytracer.cpp
#include "raytracer.h"
#include "Camera.h"
#include "Vec3.h"
#include "Ray.h"
#include "color.h"
#include "Scene.h"
#include "Sphere.h"
#include "Triangle.h"
#include "Cylinder.h"
#include "Mesh.h"
#include "Random.h"
#include "Vec3.h"
#include "json.hpp" // Include the JSON library

#include "Material.h"
#include "Light.h"
#include "PointLight.h"
#include "RectangularAreaLight.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <memory>
#include <random> 
#include <cstdlib> // For std::stoi


// For convenience
using json = nlohmann::json;
Vec3 computeColor(const Ray& ray, const Scene& scene, int depth, int maxBounces, std::string renderMode);


// Function to write a PPM image file
void writePPM(const std::string& filename, int width, int height, const std::string& pixelData) {
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Error: Could not open the file " << filename << " for writing.\n";
        return;
    }
    ofs << "P3\n" << width << " " << height << "\n255\n";
    ofs << pixelData;
    ofs.close();
}

// Function to read the scene from a JSON file
bool loadScene(const std::string& filename, Scene& scene, Camera& camera, std::string& renderMode, int& maxBounces) {
    std::ifstream inFile(filename);
    if (!inFile) {
        std::cerr << "Error: Could not open the scene file " << filename << "\n";
        return false;
    }

    json jsonData;
    try {
        inFile >> jsonData;
    } catch (const std::exception& e) {
        std::cerr << "Error: Failed to parse JSON file - " << e.what() << "\n";
        return false;
    }

    // Read render mode and max bounces
    renderMode = jsonData.value("rendermode", "phong");
    maxBounces = jsonData.value("nbounces", 1);

    // Read camera parameters
    auto jsonCamera = jsonData["camera"];
    Vec3 cameraPosition(jsonCamera["position"][0], jsonCamera["position"][1], jsonCamera["position"][2]);
    Vec3 lookAt(jsonCamera["lookAt"][0], jsonCamera["lookAt"][1], jsonCamera["lookAt"][2]);
    Vec3 upVector(jsonCamera["upVector"][0], jsonCamera["upVector"][1], jsonCamera["upVector"][2]);
    double exposure = jsonCamera.value("exposure", 1.0);
    double fov = jsonCamera["fov"];
    int imageWidth = jsonCamera["width"];
    int imageHeight = jsonCamera["height"];

    // Set up the camera

    // Read aperture and focus distance
    double aperture = jsonCamera.value("aperture", 0.0);
    double focusDistance = jsonCamera.value("focusdistance", (cameraPosition - lookAt).length());

    double aspectRatio = static_cast<double>(imageWidth) / imageHeight;

    // Update camera initialization
    camera = Camera(cameraPosition, lookAt, upVector, fov, aspectRatio, aperture, focusDistance);
    camera.setImageDimensions(imageWidth, imageHeight);
    camera.exposure = exposure;



    // Read background color
    if (jsonData["scene"].contains("backgroundcolor")) {
        auto bgColorArray = jsonData["scene"]["backgroundcolor"];
        scene.setBackgroundColor(Vec3(bgColorArray[0], bgColorArray[1], bgColorArray[2]));
    } else {
        scene.setBackgroundColor(Vec3(0.0, 0.0, 0.0)); // Default background color
    }

    // Read ambient light
    if (jsonData["scene"].contains("ambientlight")) {
        auto ambientArray = jsonData["scene"]["ambientlight"];
        scene.setAmbientLight(Vec3(ambientArray[0], ambientArray[1], ambientArray[2]));
    } else {
        // Default ambient light
        scene.setAmbientLight(Vec3(0.5, 0.5, 0.5)); // Default Ambient Light
    }

    // Read lights
    if (jsonData["scene"].contains("lightsources")) {
        for (const auto& lightData : jsonData["scene"]["lightsources"]) {
            std::string type = lightData["type"];
            if (type == "pointlight") {
                Vec3 position(lightData["position"][0], lightData["position"][1], lightData["position"][2]);
                Vec3 intensity(lightData["intensity"][0], lightData["intensity"][1], lightData["intensity"][2]);
                scene.addLight(std::make_shared<PointLight>(position, intensity));
            } else if (type == "arealight") {
                Vec3 position(lightData["position"][0], lightData["position"][1], lightData["position"][2]);
                Vec3 u(lightData["u"][0], lightData["u"][1], lightData["u"][2]);
                Vec3 v(lightData["v"][0], lightData["v"][1], lightData["v"][2]);
                Vec3 intensity(lightData["intensity"][0], lightData["intensity"][1], lightData["intensity"][2]);
                scene.addLight(std::make_shared<RectangularAreaLight>(position, u, v, intensity));
            } else {
                std::cerr << "Warning: Unknown light type '" << type << "'\n";
            }
        }
    }


    // Read shapes
    for (const auto& shapeData : jsonData["scene"]["shapes"]) {
        std::string type = shapeData["type"];
        Material material;

        // Parse material if present
        if (shapeData.contains("material")) {
            const auto& materialData = shapeData["material"];
            double ks = materialData["ks"];
            double kd = materialData["kd"];
            double specularExponent = materialData["specularexponent"];
            Vec3 diffuseColor(materialData["diffusecolor"][0], materialData["diffusecolor"][1], materialData["diffusecolor"][2]);
            Vec3 specularColor(materialData["specularcolor"][0], materialData["specularcolor"][1], materialData["specularcolor"][2]);
            bool isReflective = materialData["isreflective"];
            double reflectivity = materialData["reflectivity"];
            bool isRefractive = materialData["isrefractive"];
            double refractiveIndex = materialData["refractiveindex"];
            double dissolve = materialData.value("dissolve", 1.0); // Default to 1.0 if not specified
            double metallic = materialData.value("metallic", 0.0);  // Default to 0.0 if not specified
            double roughness = materialData.value("roughness", 0.5); // Default to 0.5 if not specified
            

            material =  Material(ks, kd, specularExponent, diffuseColor, specularColor,
                    isReflective, reflectivity, isRefractive, refractiveIndex);

            // Load texture if specified
            if (materialData.contains("texture")) {
                std::string textureFilename = materialData["texture"];
                material.texture = std::make_shared<Texture>(textureFilename);
            }

            if (materialData.contains("scaleU")) {
                material.scaleU = materialData["scaleU"];
            }
            if (materialData.contains("scaleV")) {
                material.scaleV = materialData["scaleV"];
            }
        } else {
            // Default material if none provided
            material = Material();
        }

        // Create shape based on type
        if (type == "sphere") {
            Vec3 center(shapeData["center"][0], shapeData["center"][1], shapeData["center"][2]);
            double radius = shapeData["radius"];
            scene.addShape(std::make_shared<Sphere>(center, radius, material));
        } else if (type == "triangle") {
            Vec3 v0(shapeData["v0"][0], shapeData["v0"][1], shapeData["v0"][2]);
            Vec3 v1(shapeData["v1"][0], shapeData["v1"][1], shapeData["v1"][2]);
            Vec3 v2(shapeData["v2"][0], shapeData["v2"][1], shapeData["v2"][2]);

            // Parse UV coordinates
            Vec2 uv0(0.0, 0.0); // Default UVs if not specified
            Vec2 uv1(1.0, 0.0);
            Vec2 uv2(0.0, 1.0);

            if (shapeData.contains("uv0")) {
                uv0 = Vec2(shapeData["uv0"][0], shapeData["uv0"][1]);
            }
            if (shapeData.contains("uv1")) {
                uv1 = Vec2(shapeData["uv1"][0], shapeData["uv1"][1]);
            }
            if (shapeData.contains("uv2")) {
                uv2 = Vec2(shapeData["uv2"][0], shapeData["uv2"][1]);
            }
        

            scene.addShape(std::make_shared<Triangle>(v0, v1, v2, uv0, uv1, uv2, material));

        } else if (type == "cylinder") {
            Vec3 center(shapeData["center"][0], shapeData["center"][1], shapeData["center"][2]);
            Vec3 axis(shapeData["axis"][0], shapeData["axis"][1], shapeData["axis"][2]);
            double radius = shapeData["radius"];
            double height = shapeData["height"];
            scene.addShape(std::make_shared<Cylinder>(center, axis, radius, height, material));
        }else if (type == "mesh") {
            std::string objFilename = shapeData["filename"];
            // Create the mesh and add it to the scene
            auto mesh = std::make_shared<Mesh>(objFilename, material, scene);
            scene.addShape(mesh);
                
            // Retrieve and add area lights from the mesh to the scene's light list
            for (const auto& areaLight : mesh->getAreaLights()) {
                scene.addLight(areaLight);
            }
        }else {
                std::cerr << "Warning: Unknown shape type '" << type << "'\n";
        }
    }

    return true;
}

Material parseMaterial(const json& materialData) {
    double ks = materialData["ks"];
    double kd = materialData["kd"];
    double specularExponent = materialData["specularexponent"];
    Vec3 diffuseColor(materialData["diffusecolor"][0], materialData["diffusecolor"][1], materialData["diffusecolor"][2]);
    Vec3 specularColor(materialData["specularcolor"][0], materialData["specularcolor"][1], materialData["specularcolor"][2]);
    bool isReflective = materialData["isreflective"];
    double reflectivity = materialData["reflectivity"];
    bool isRefractive = materialData["isrefractive"];
    double refractiveIndex = materialData["refractiveindex"];

    return Material(ks, kd, specularExponent, diffuseColor, specularColor,
                    isReflective, reflectivity, isRefractive, refractiveIndex);
}

double fresnel(const Vec3& I, const Vec3& N, double eta_i, double eta_t) {
    double cosi = std::clamp(I.dot(N), -1.0, 1.0);
    double etai = eta_i, etat = eta_t;
    if (cosi > 0) {
        std::swap(etai, etat);
    }
    // Compute sin of transmission angle using Snell's law
    double sint = etai / etat * sqrt(std::max(0.0, 1 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1.0) {
        return 1.0;
    } else {
        double cost = sqrt(std::max(0.0, 1 - sint * sint));
        cosi = fabs(cosi);
        double Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        double Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        return (Rs * Rs + Rp * Rp) / 2.0;
    }
}

Vec3 handleReflection(const Ray& ray, const HitRecord& hitRecord, const Scene& scene, int depth, int maxBounces, double rrProbability, const std::string& renderMode) {
    Vec3 reflectionDir = ray.direction - 2 * ray.direction.dot(hitRecord.normal) * hitRecord.normal;
    reflectionDir = reflectionDir.normalized();
    Ray reflectionRay(hitRecord.point + hitRecord.normal * 0.001, reflectionDir);
    Vec3 reflectionColor = computeColor(reflectionRay, scene, depth + 1, maxBounces, renderMode);
    return reflectionColor * hitRecord.material.reflectivity;
}

Vec3 handleRefraction(const Ray& ray, const HitRecord& hitRecord, const Scene& scene, int depth, int maxBounces, double rrProbability, const std::string& renderMode) {
    double eta_i = 1.0; // Refractive index of air
    double eta_t = hitRecord.material.refractiveIndex;

    Vec3 normal = hitRecord.normal;
    double cosi = ray.direction.dot(normal);
    if (cosi > 0.0) {
        // Ray is inside the object
        normal = -normal;
        cosi = -cosi;
        std::swap(eta_i, eta_t);
    }

    double fresnelCoefficient = fresnel(ray.direction, normal, eta_i, eta_t);
    double reflectionCoefficient = fresnelCoefficient;
    double refractionCoefficient = 1.0 - reflectionCoefficient;

    Vec3 color(0.0, 0.0, 0.0);

    // Reflection
    if (reflectionCoefficient > 0.0) {
        Vec3 reflectionDir = ray.direction - 2 * ray.direction.dot(hitRecord.normal) * hitRecord.normal;
        reflectionDir = reflectionDir.normalized();
        Ray reflectionRay(hitRecord.point + hitRecord.normal * 0.001, reflectionDir);
        Vec3 reflectionColor = computeColor(reflectionRay, scene, depth + 1, maxBounces, renderMode);
        color += reflectionColor * reflectionCoefficient;
    }

    // Refraction
    double etaRatio = eta_i / eta_t;
    double cosTheta = -ray.direction.dot(normal);
    double k = 1.0 - etaRatio * etaRatio * (1.0 - cosTheta * cosTheta);
    if (k >= 0.0) {
        // Refraction occurs
        Vec3 refractionDir = (ray.direction * etaRatio + normal * (etaRatio * cosTheta - sqrt(k))).normalized();
        Ray refractionRay(hitRecord.point - normal * 0.001, refractionDir);
        Vec3 refractionColor = computeColor(refractionRay, scene, depth + 1, maxBounces, renderMode);
        // Adjust the transmission based on roughness
        double specularTransmissionFactor = 1.0 - hitRecord.material.roughness;
        color += Vec3::hadamardProduct(hitRecord.material.transmissionFilter, refractionColor) * refractionCoefficient * specularTransmissionFactor;
    } else {
        // Total internal reflection
        // Already handled in reflection above
    }

    return color;
}

// Henyey-Greenstein phase function sampling
Vec3 samplePhaseFunction(const Vec3& normal, double g) {
    double u1 = Random::uniformDouble();
    double u2 = Random::uniformDouble();

    double cosTheta;
    if (fabs(g) < 1e-3) {
        cosTheta = 1.0 - 2.0 * u1;
    } else {
        double sqrTerm = (1.0 - g * g) / (1.0 - g + 2.0 * g * u1);
        cosTheta = (1.0 + g * g - sqrTerm * sqrTerm) / (2.0 * g);
    }
    double sinTheta = sqrt(1.0 - cosTheta * cosTheta);
    double phi = 2.0 * M_PI * u2;

    // Spherical coordinates to Cartesian
    double x = sinTheta * cos(phi);
    double y = sinTheta * sin(phi);
    double z = cosTheta;

    // Orthonormal basis
    Vec3 w = normal;
    Vec3 helper = (fabs(w.x) > 0.1) ? Vec3(0.0, 1.0, 0.0) : Vec3(1.0, 0.0, 0.0);
    Vec3 u = w.cross(helper).normalized();
    Vec3 v = w.cross(u);

    // Transform to world coordinates
    Vec3 direction = (u * x + v * y + w * z).normalized();
    return direction;
}


Vec3 computeColor(const Ray& ray, const Scene& scene, int depth, int maxBounces, std::string renderMode) {

    if (depth >= maxBounces) {
        return Vec3(0.0, 0.0, 0.0); // Return black if maximum recursion depth reached
    }
    
    HitRecord hitRecord;
    if (scene.intersect(ray, 0.001, std::numeric_limits<double>::infinity(), hitRecord)) {
        if (renderMode == "pathtracing") {

            // Initialize the color
            Vec3 color(0.0, 0.0, 0.0);

            // Handle emission
            // **Handle emission**
            if (hitRecord.material.emission.lengthSquared() > 0.0) {
                color += hitRecord.material.emission;
                // Do not return immediately; allow further lighting computations
            }

             double rrProbability = 1.0;
            if (depth >= 3) {
                rrProbability = 0.8; // Adjust as needed
                if (Random::uniformDouble() > rrProbability) {
                    return Vec3(0.0, 0.0, 0.0);
                }
            }

            // Compute the diffuse color at the hit point
            Vec3 diffuseColor = hitRecord.material.diffuseColor;
            if (hitRecord.material.texture) {
                double scaledU = hitRecord.u * hitRecord.material.scaleU;
                double scaledV = hitRecord.v * hitRecord.material.scaleV;
                diffuseColor = hitRecord.material.texture->getColor(scaledU, scaledV);
            }

            int lightSamplesPerAreaLight = 64;
            // **Direct Illumination**
            // For each light source, compute the direct illumination
            for (const auto& light : scene.getLights()) {
                Vec3 directLight(0.0, 0.0, 0.0);

                Vec3 lightPoint, lightNormal;
                double pdf;

                // Sample the light
                if (light->sample(hitRecord.point, lightPoint, lightNormal, pdf)) {
                    Vec3 lightDir = lightPoint - hitRecord.point;
                    double distanceSquared = lightDir.lengthSquared();
                    double distance = sqrt(distanceSquared);
                    lightDir /= distance; // Normalize

                    // Shadow Ray Origin
                    Vec3 shadowOrigin = hitRecord.point + hitRecord.normal * 0.001;

                    // Shadow Ray
                    Ray shadowRay(shadowOrigin, lightDir);
                    HitRecord tempRecord;

                    // Check for occlusion
                    if (scene.intersect(shadowRay, 0.001, distance - 0.001, tempRecord)) {
                        continue; // Occluded, skip this light sample
                    }


                   double cosTheta = std::max(0.0, hitRecord.normal.dot(lightDir));
                    double cosThetaLight = (lightNormal.lengthSquared() > 0.0) ? std::max(0.0, lightNormal.dot(-lightDir)) : 1.0; // For point lights, lightNormal might be zero

                    Vec3 radiance = light->getRadiance();

                    Vec3 brdf = (diffuseColor / M_PI) * hitRecord.material.kd;

                    double G = (cosTheta * cosThetaLight) / distanceSquared;

                    Vec3 sampleContribution = Vec3::hadamardProduct(brdf, radiance) * G / pdf;

                    directLight += sampleContribution;

                }

                // Accumulate direct lighting
                color += directLight;
            }

            // **Indirect Illumination**
            // Sample a direction using BRDF, passing the diffuseColor
            Vec3 outDir;
            double pdf;
            Vec3 brdf = hitRecord.material.sampleBRDF(hitRecord.normal, -ray.direction.normalized(), outDir, pdf, diffuseColor);

            if (pdf > 0.0 && brdf.lengthSquared() > 0.0) {
                double cosTheta = std::max(0.0, hitRecord.normal.dot(outDir));
                Ray scattered(hitRecord.point + hitRecord.normal * 0.001, outDir);
                Vec3 incomingRadiance = computeColor(scattered, scene, depth + 1, maxBounces, renderMode);

                double denominator = pdf * rrProbability;
                if (denominator > 0.0) {
                    Vec3 indirect = Vec3::hadamardProduct(brdf, incomingRadiance) * (cosTheta / denominator);
                    color += indirect;
                }
            }

            if (hitRecord.material.transmissionFilter.lengthSquared() > 0.0 && hitRecord.material.transparency > 0.0) {
                Vec3 scatteredDir;
                double pdfDiffuse;
                Vec3 brdfDiffuse = hitRecord.material.sampleDiffuseTransmission(hitRecord.normal, -ray.direction.normalized(), scatteredDir, pdfDiffuse);

                if (pdfDiffuse > 0.0 && brdfDiffuse.lengthSquared() > 0.0) {
                    double cosTheta = std::max(0.0, hitRecord.normal.dot(scatteredDir));
                    Vec3 normal = hitRecord.normal;
                    Vec3 scatteredOrigin = hitRecord.point + normal * 0.001;
                    Ray scatteredRay(scatteredOrigin, scatteredDir);
                    Vec3 incomingScatteredRadiance = computeColor(scatteredRay, scene, depth + 1, maxBounces, renderMode);

                    // Apply Transmission Filter and Scale by Transparency
                    Vec3 filteredRadiance = Vec3::hadamardProduct(brdfDiffuse, incomingScatteredRadiance);
                    color += Vec3::hadamardProduct(filteredRadiance, hitRecord.material.transmissionFilter) * (cosTheta / pdfDiffuse) * hitRecord.material.transparency;
                }
            }

            
            // Handle BSSDRF**
            if (hitRecord.material.sigma_s > 0.0) {
                double sigma_t = hitRecord.material.sigma_a + hitRecord.material.sigma_s;
                double u = Random::uniformDouble();
                u = std::max(1e-6, u); // Prevent log(0)
                double scatterDistance = -log(1.0 - u) / sigma_t;

                // **Compute the Scattered Origin from the Hit Point**
                Vec3 scatteredOrigin = hitRecord.point + ray.direction.normalized() * scatterDistance;

                // **Sample a New Direction Using the Phase Function**
                Vec3 scatteredDir = samplePhaseFunction(hitRecord.normal, hitRecord.material.g);

                Ray scatteredRay(scatteredOrigin + scatteredDir * 0.001, scatteredDir);
                
                // **Recursive Call to Compute Incoming Radiance**
                Vec3 incomingRadiance = computeColor(scatteredRay, scene, depth + 1, maxBounces, renderMode);

                // **Apply Beer-Lambert Law for Absorption**
                Vec3 absorption = Vec3(
                    exp(-hitRecord.material.sigma_a * scatterDistance),
                    exp(-hitRecord.material.sigma_a * scatterDistance),
                    exp(-hitRecord.material.sigma_a * scatterDistance)
                );

                // **Scattering Contribution**
                double scatteringProb = hitRecord.material.sigma_s / sigma_t;
                color += Vec3::hadamardProduct(absorption, incomingRadiance) * scatteringProb;
            }

            // **Handle Specular Transmission (Specular Transparency)**
            if (hitRecord.material.isRefractive) {
                // Handle specular transmission via refraction
                Vec3 specularTransmission = handleRefraction(ray, hitRecord, scene, depth, maxBounces, rrProbability, renderMode);
                color += specularTransmission;
            }

            // **Handle Reflective Materials**
            if (hitRecord.material.isReflective && hitRecord.material.reflectivity > 0.0) {
                Vec3 reflection = handleReflection(ray, hitRecord, scene, depth, maxBounces, rrProbability, renderMode);
                color += reflection;
            }

            return color;

        }else if (renderMode == "phong"){
            Vec3 diffuseColor = hitRecord.material.diffuseColor;
            /*if (hitRecord.material.texture) {
                diffuseColor = hitRecord.material.texture->getColor(hitRecord.u, hitRecord.v);
            }*/

            if (hitRecord.material.texture) {
                double scaledU = hitRecord.u * hitRecord.material.scaleU;
                double scaledV = hitRecord.v * hitRecord.material.scaleV;
                diffuseColor = hitRecord.material.texture->getColor(scaledU, scaledV);
                /*std::cout << "Material specularColor: (" 
                    << diffuseColor.x << ", " 
                    << diffuseColor.y << ", " 
                    << diffuseColor.z << ")\n";*/
            }

            Vec3 viewDir = -ray.direction.normalized();
            Vec3 ambient = Vec3::hadamardProduct(diffuseColor, scene.getAmbientLight()) * hitRecord.material.kd;
            Vec3 color = ambient;

            for (const auto& light : scene.getLights()) {
                Vec3 lightDir = (light->getPosition() - hitRecord.point).normalized();

                // Shadow check
                Ray shadowRay(hitRecord.point + hitRecord.normal * 0.001, lightDir);
                HitRecord tempRecord;
                if (scene.intersect(shadowRay, 0.001, (light->getPosition() - hitRecord.point).length(), tempRecord)) {
                    continue; // In shadow, skip this light
                }

                // Diffuse component
                double diff = std::max(hitRecord.normal.dot(lightDir), 0.0);
                double diffuseFactor = hitRecord.material.kd * diff;
                Vec3 diffuse = diffuseColor * diffuseFactor;
                diffuse = Vec3::hadamardProduct(diffuse, light->getIntensity());

                // Specular component
                Vec3 halfDir = (lightDir + viewDir).normalized();
                double specAngle = std::max(hitRecord.normal.dot(halfDir), 0.0);
                double specularFactor = hitRecord.material.ks * pow(specAngle, hitRecord.material.specularExponent);
                Vec3 specular = hitRecord.material.specularColor * specularFactor;
                specular = Vec3::hadamardProduct(specular, light->getIntensity());

                color += diffuse + specular;
            }
            // Add emission
            color += hitRecord.material.emission;

            // **Handle Transparency**
            if (hitRecord.material.transmissionFilter.lengthSquared() > 0.0 && hitRecord.material.transparency > 0.0 && depth < maxBounces) {
                // Compute the transmitted color by modulating with transmissionFilter
                Vec3 transmittedColor = hitRecord.material.transmissionFilter * hitRecord.material.transparency;
                Ray transmittedRay(hitRecord.point + ray.direction.normalized() * 0.001, ray.direction.normalized());
                Vec3 incomingTransmittedRadiance = computeColor(transmittedRay, scene, depth + 1, maxBounces, renderMode);
                color += Vec3::hadamardProduct(transmittedColor, incomingTransmittedRadiance);
            }


            Vec3 reflectionColor(0.0, 0.0, 0.0);
            Vec3 refractionColor(0.0, 0.0, 0.0);

            // Check if the material is reflective or refractive
            if (depth < maxBounces) {
                // Handle refractive materials
                if (hitRecord.material.isRefractive) {
                    double eta_i = 1.0; // Refractive index of air
                    double eta_t = hitRecord.material.refractiveIndex;

                    Vec3 normal = hitRecord.normal;
                    double cosi = ray.direction.dot(normal);
                    if (cosi > 0.0) {
                        // Ray is inside the object
                        normal = -normal;
                        cosi = -cosi;
                        std::swap(eta_i, eta_t);
                    }

                    double fresnelCoefficient = fresnel(ray.direction, normal, eta_i, eta_t);
                    double reflectionCoefficient = fresnelCoefficient;
                    double refractionCoefficient = 1.0 - reflectionCoefficient;

                    // Reflection
                    if (reflectionCoefficient > 0.0) {
                        Vec3 reflectionDir = ray.direction - 2 * ray.direction.dot(hitRecord.normal) * hitRecord.normal;
                        reflectionDir = reflectionDir.normalized();
                        Ray reflectionRay(hitRecord.point + hitRecord.normal * 0.001, reflectionDir);
                        reflectionColor = computeColor(reflectionRay, scene, depth + 1, maxBounces, renderMode) * reflectionCoefficient;
                    }

                    // Refraction
                    double etaRatio = eta_i / eta_t;
                    double cosTheta = -ray.direction.dot(normal);
                    double k = 1.0 - etaRatio * etaRatio * (1.0 - cosi * cosi);
                    if (k >= 0.0) {
                        // Refraction occurs
                        Vec3 refractionDir = (ray.direction * etaRatio + normal * (etaRatio * cosi - sqrt(k))).normalized();
                        Ray refractionRay(hitRecord.point - normal * 0.001, refractionDir);
                        refractionColor = computeColor(refractionRay, scene, depth + 1, maxBounces, renderMode) * refractionCoefficient;
                    } else {
                        // Total internal reflection
                        Vec3 reflectionDir = ray.direction - 2 * ray.direction.dot(hitRecord.normal) * hitRecord.normal;
                        reflectionDir = reflectionDir.normalized();
                        Ray reflectionRay(hitRecord.point + hitRecord.normal * 0.001, reflectionDir);
                        reflectionColor += computeColor(reflectionRay, scene, depth + 1, maxBounces, renderMode) * reflectionCoefficient;
                    }

                    // Combine reflection and refraction
                    color += reflectionColor + refractionColor;
                }
                // Handle reflective but non-refractive materials
                else if (hitRecord.material.isReflective && hitRecord.material.reflectivity > 0.0) {
                    Vec3 reflectionDir = ray.direction - 2 * ray.direction.dot(hitRecord.normal) * hitRecord.normal;
                    reflectionDir = reflectionDir.normalized();
                    Ray reflectionRay(hitRecord.point + hitRecord.normal * 0.001, reflectionDir);
                    reflectionColor = computeColor(reflectionRay, scene, depth + 1, maxBounces, renderMode);

                    // Combine reflection with current color using material's reflectivity
                    color = color * (1.0 - hitRecord.material.reflectivity) + reflectionColor * hitRecord.material.reflectivity;
                }
            }
            return color;

        }else if (renderMode == "binary"){
            return Vec3(1.0f, 0.0f, 0.0f);
            //return hitRecord.material.diffuseColor;
        }else{
            // Default to Phong shading if unknown render mode
            std::cerr << "Warning: Unknown render mode '" << renderMode << "'. Defaulting to 'phong'.\n";
            // (Optionally, you can set a default color or handle differently)
            return computeColor(ray, scene, depth, maxBounces, "phong");
        }
        // Add ambient light if necessary (not included in this example)
    } else {
        if (renderMode == "phong"){
            return scene.getBackgroundColor();
        }else if (renderMode == "binary"){
            return Vec3(0.0f, 0.0f, 0.0f);
        }else if (renderMode == "pathtracing"){
            return scene.getBackgroundColor();
        }else {
            // Default return value for unknown render modes
            return Vec3(0.0, 0.0, 0.0);
        }
    }
}


Vec3 applyExposure(const Vec3& color, double exposure) {
    // Apply exposure scaling
    Vec3 adjustedColor = Vec3(
        1.0 - exp(-color.x * exposure),
        1.0 - exp(-color.y * exposure),
        1.0 - exp(-color.z * exposure)
    );
    return adjustedColor;
}


// Define whether to use BVH
#define USE_BVH

int main(int argc, char* argv[]) {
    // Check if the correct number of arguments are provided
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <path_to_json_file> <samples_per_pixel><output_file_name>\n";
        return 1;
    }

    // Parse arguments
    std::string jsonFilePath = argv[1];
    int samplesPerPixel = std::stoi(argv[2]); // Convert the second argument to an integer
    std::string outputFileName = argv[3];

    // Validate samplesPerPixel
    if (samplesPerPixel <= 0) {
        std::cerr << "Error: Samples per pixel must be a positive integer.\n";
        return 1;
    }

    // Scene setup
    Scene scene;
    Camera camera;
    std::string renderMode;
    int maxBounces;

    // Load the scene from the JSON file
    if (!loadScene(jsonFilePath, scene, camera, renderMode, maxBounces)) {
        std::cerr << "Failed to load scene from file: " << jsonFilePath << "\n";
        return 1;
    }

    // **Start Timing Here**
    auto startTime = std::chrono::high_resolution_clock::now();

    // Build the BVH conditionally
#ifdef USE_BVH
    auto bvhStart = std::chrono::high_resolution_clock::now();
    scene.buildBVH();
    auto bvhEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> bvhElapsed = bvhEnd - bvhStart;
    std::cout << "BVH Construction time: " << bvhElapsed.count() << " seconds.\n";
#endif

    int imageWidth = camera.imageWidth;
    int imageHeight = camera.imageHeight;

    // Random number generator setup
    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    // Prepare to collect pixel data
    std::stringstream pixelDataStream;

    // Render the scene
    for (int j = imageHeight - 1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < imageWidth; ++i) {
            Vec3 pixelColor(0.0, 0.0, 0.0);
            int sqrtSamples = static_cast<int>(std::sqrt(samplesPerPixel));
            for (int s = 0; s < samplesPerPixel; ++s) {
                int x = s % sqrtSamples;
                int y = s / sqrtSamples;
                double u = (i + (x + distribution(generator)) / sqrtSamples) / (imageWidth - 1);
                double v = (j + (y + distribution(generator)) / sqrtSamples) / (imageHeight - 1);

                Ray ray = camera.getRay(u, v);
                Vec3 sampleColor = computeColor(ray, scene, 0, maxBounces, renderMode);
                pixelColor += sampleColor;
            }

            // Average the color over the number of samples
            pixelColor /= static_cast<double>(samplesPerPixel);

            // Apply exposure adjustment
            Vec3 adjustedColor = applyExposure(pixelColor, camera.exposure);

            // Write the pixel color to the stream
            writeColor(pixelDataStream, adjustedColor);
        }
    }
    std::cerr << "\nDone.\n";

    // **End Timing Here**
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = endTime - startTime;

    // Output the rendering time
#ifdef USE_BVH
    std::cout << "Rendering time with BVH: " << elapsed.count() << " seconds.\n";
#else
    std::cout << "Rendering time without BVH: " << elapsed.count() << " seconds.\n";
#endif

    // Write the image to a file
    writePPM(outputFileName, imageWidth, imageHeight, pixelDataStream.str());

    return 0;
}