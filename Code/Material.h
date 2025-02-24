#ifndef MATERIAL_H
#define MATERIAL_H

#include "Vec3.h"
#include "Texture.h"
#include <memory>

class Material {
public:
    double ks;                 // Specular reflection coefficient
    double kd;                 // Diffuse reflection coefficient
    double specularExponent;   // Specular exponent
    Vec3 diffuseColor;         // Diffuse color
    Vec3 specularColor;        // Specular color
    bool isReflective;
    double reflectivity;
    bool isRefractive;
    double refractiveIndex;
    double kt;
    double transparency;
    Vec3 emission;
    Vec3 transmissionFilter;   // Transmission filter color (Tf)
    double sigma_a = 0;
    double sigma_s = 0;
    double g = 0;

    double scaleU = 1.0;
    double scaleV = 1.0;

    double roughness = 0.5;
    double metallic = 0.0;

    std::shared_ptr<Texture> texture;

    // Constructors
    Material();
    Material(double ks_, double kd_, double specularExponent_,
             const Vec3& diffuseColor_, const Vec3& specularColor_,
             bool isReflective_, double reflectivity_,
             bool isRefractive_, double refractiveIndex_);

    // BRDF Sampling Method
    Vec3 sampleBRDF(const Vec3& normal, const Vec3& viewDir, Vec3& outDir,
                    double& pdf, const Vec3& diffuseColor) const;

    Vec3 sampleDiffuseTransmission(const Vec3& normal, const Vec3& viewDir,
                                   Vec3& outDir, double& pdf) const;         
    bool refract(const Vec3& incident, const Vec3& normal, Vec3& refracted, double eta_t) const;
    Vec3 reflect(const Vec3& incident, const Vec3& normal) const;      
    Vec3 evaluateBRDF(const Vec3& normal, const Vec3& viewDir, const Vec3& lightDir, const Vec3& diffuseColor) const;

private:
    // Sampling functions
    Vec3 sampleCosineHemisphere(const Vec3& normal) const;
    Vec3 samplePhongLobe(const Vec3& reflectionDir, double exponent) const;

    // GGX sampling functions
    Vec3 sampleGGX(const Vec3& normal, const Vec3& viewDir, double roughness) const;
    double ggxNormalDistribution(const Vec3& h, const Vec3& n, double roughness) const;
    double ggxGeometryFunction(const Vec3& v, const Vec3& l, const Vec3& n, double roughness) const;

    // Helper functions
    void sampleDiffuse(const Vec3& normal, Vec3& outDir, double& pdf) const;
    void sampleSpecular(const Vec3& normal, const Vec3& viewDir,
                        Vec3& outDir, double& pdf) const;


};

#endif // MATERIAL_H
