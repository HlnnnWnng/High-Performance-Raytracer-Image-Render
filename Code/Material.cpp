#include "Material.h"
#include "Random.h"
#include <cmath>
#include <algorithm>

Material::Material()
    : ks(0.0), kd(0.0), specularExponent(0.0),
      diffuseColor(0.0, 0.0, 0.0), specularColor(0.0, 0.0, 0.0),
      isReflective(false), reflectivity(0.0),
      isRefractive(false), refractiveIndex(1.0) {}

Material::Material(double ks_, double kd_, double specularExponent_,
                   const Vec3& diffuseColor_, const Vec3& specularColor_,
                   bool isReflective_, double reflectivity_,
                   bool isRefractive_, double refractiveIndex_)
    : ks(ks_), kd(kd_), specularExponent(specularExponent_),
      diffuseColor(diffuseColor_), specularColor(specularColor_),
      isReflective(isReflective_), reflectivity(reflectivity_),
      isRefractive(isRefractive_), refractiveIndex(refractiveIndex_) {}



Vec3 fresnelSchlick(const Vec3& F0, double cosTheta) {
    return F0 + (Vec3(1.0, 1.0, 1.0) - F0) * pow(1.0 - cosTheta, 5.0);
}


void orthonormalBasis(const Vec3& n, Vec3& tangent, Vec3& bitangent) {
    if (fabs(n.x) > fabs(n.z)) {
        tangent = Vec3(-n.y, n.x, 0.0).normalized();
    } else {
        tangent = Vec3(0.0, -n.z, n.y).normalized();
    }
    bitangent = n.cross(tangent);
}

// Refract function
bool Material::refract(const Vec3& incident, const Vec3& normal, Vec3& refracted, double eta_t) const {
    double eta_i = 1.0; // Refractive index of air
    double cosi = -std::max(-1.0, std::min(1.0, incident.dot(normal)));
    double eta = eta_i / eta_t;
    if (cosi < 0) {
        // Ray is inside the medium
        cosi = -cosi;
        eta = eta_t / eta_i;
    }
    double k = 1 - eta * eta * (1 - cosi * cosi);
    if (k < 0) {
        return false; // Total internal reflection
    } else {
        refracted = incident * eta + normal * (eta * cosi - sqrt(k));
        return true;
    }
}

// Reflect function
Vec3 Material::reflect(const Vec3& incident, const Vec3& normal) const {
    return incident - 2 * incident.dot(normal) * normal;
}


Vec3 Material::sampleDiffuseTransmission(const Vec3& normal, const Vec3& viewDir,
                                         Vec3& outDir, double& pdf) const {
    // Cosine-weighted hemisphere sampling for diffuse transmission
    double u1 = Random::uniformDouble();
    double u2 = Random::uniformDouble();

    double r = sqrt(u1);
    double theta = 2.0 * M_PI * u2;

    double x = r * cos(theta);
    double y = r * sin(theta);
    double z = sqrt(std::max(0.0, 1.0 - u1));

    // Orthonormal basis (u, v, w)
    Vec3 w = normal;
    Vec3 helper = (fabs(w.x) > 0.1) ? Vec3(0.0, 1.0, 0.0) : Vec3(1.0, 0.0, 0.0);
    Vec3 u = w.cross(helper).normalized();
    Vec3 v = w.cross(u);

    // Convert to world coordinates
    Vec3 direction = u * x + v * y + w * z;
    outDir = direction.normalized();
    pdf = std::max(0.0, outDir.dot(normal)) / M_PI;

    // BRDF for diffuse transmission
    return transmissionFilter * transparency / M_PI;
}


Vec3 Material::sampleBRDF(const Vec3& normal, const Vec3& viewDir,
                          Vec3& outDir, double& pdf, const Vec3& diffuseColor) const {
    // Adjust kd and ks based on metallic property
    double totalReflectance = kd + ks;
    double adjustedKd = kd;
    double adjustedKs = ks;
    Vec3 adjustedSpecularColor = specularColor;

    if (metallic > 0.0) {
        adjustedKd = kd * (1.0 - metallic);
        adjustedKs = ks * metallic;
        adjustedSpecularColor = diffuseColor; // For metallic materials, specular color is diffuse color
    }

    totalReflectance = adjustedKd + adjustedKs;

    if (totalReflectance == 0.0) {
        pdf = 0.0;
        return Vec3(0.0, 0.0, 0.0);
    }

    double randomSample = Random::uniformDouble();

    if (randomSample < adjustedKd / totalReflectance) {
        // Diffuse reflection
        outDir = sampleCosineHemisphere(normal);
        pdf = std::max(0.0, outDir.dot(normal)) / M_PI;

        Vec3 brdf = diffuseColor * adjustedKd / M_PI;
        return brdf;
    } else {
        // Specular reflection using GGX microfacet model
        outDir = sampleGGX(normal, viewDir, roughness);
        Vec3 halfVector = (viewDir + outDir).normalized();
        double D = ggxNormalDistribution(halfVector, normal, roughness);
        double G = ggxGeometryFunction(viewDir, outDir, normal, roughness);
        Vec3 F = fresnelSchlick(adjustedSpecularColor, std::max(0.0, halfVector.dot(viewDir)));

        double NdotV = std::max(0.0, normal.dot(viewDir));
        double NdotL = std::max(0.0, normal.dot(outDir));
        double denominator = 4.0 * NdotV * NdotL + 1e-6;
        Vec3 brdf = (D * G * F) / denominator;

        // PDF calculation
        double pdfDenominator = 4.0 * fabs(viewDir.dot(halfVector));
        pdf = D * normal.dot(halfVector) / std::max(pdfDenominator, 1e-6);

        return brdf * adjustedKs;
    }
}
// Cosine-weighted hemisphere sampling for diffuse reflection
Vec3 Material::sampleCosineHemisphere(const Vec3& normal) const {
    double u1 = Random::uniformDouble();
    double u2 = Random::uniformDouble();

    double r = sqrt(u1);
    double theta = 2 * M_PI * u2;

    double x = r * cos(theta);
    double y = r * sin(theta);
    double z = sqrt(std::max(0.0, 1.0 - u1));

    // Orthonormal basis (u, v, w)
    Vec3 w = normal;
    Vec3 helper = (fabs(w.x) > 0.1) ? Vec3(0.0, 1.0, 0.0) : Vec3(1.0, 0.0, 0.0);
    Vec3 u = w.cross(helper).normalized();
    Vec3 v = w.cross(u);

    // Convert to world coordinates
    Vec3 direction = u * x + v * y + w * z;
    return direction.normalized();
}

// Phong lobe sampling for specular reflection
Vec3 Material::samplePhongLobe(const Vec3& reflectionDir, double exponent) const {
    double u1 = Random::uniformDouble();
    double u2 = Random::uniformDouble();

    double cosTheta = pow(u1, 1.0 / (exponent + 1.0));
    double sinTheta = sqrt(1.0 - cosTheta * cosTheta);
    double phi = 2 * M_PI * u2;

    double x = sinTheta * cos(phi);
    double y = sinTheta * sin(phi);
    double z = cosTheta;

    // Orthonormal basis around reflection direction
    Vec3 w = reflectionDir;
    Vec3 helper = (fabs(w.x) > 0.1) ? Vec3(0.0, 1.0, 0.0) : Vec3(1.0, 0.0, 0.0);
    Vec3 u = w.cross(helper).normalized();
    Vec3 v = w.cross(u);

    // Convert to world coordinates
    Vec3 direction = u * x + v * y + w * z;
    return direction.normalized();
}

// Helper functions (optional if you prefer to keep them)
void Material::sampleDiffuse(const Vec3& normal, Vec3& outDir, double& pdf) const {
    outDir = sampleCosineHemisphere(normal);
    pdf = std::max(0.0, outDir.dot(normal)) / M_PI;
}

void Material::sampleSpecular(const Vec3& normal, const Vec3& viewDir,
                              Vec3& outDir, double& pdf) const {
    Vec3 reflectionDir = (viewDir - 2 * viewDir.dot(normal) * normal).normalized();
    outDir = samplePhongLobe(reflectionDir, specularExponent);
    double cosAlpha = std::max(0.0, outDir.dot(reflectionDir));
    pdf = (specularExponent + 1) * pow(cosAlpha, specularExponent) / (2 * M_PI);
}

Vec3 Material::sampleGGX(const Vec3& normal, const Vec3& viewDir, double roughness) const {
    // Generate random numbers
    double u1 = Random::uniformDouble();
    double u2 = Random::uniformDouble();

    // Compute theta and phi
    double alpha = roughness * roughness;
    double phi = 2.0 * M_PI * u1;
    double cosTheta = sqrt((1.0 - u2) / (1.0 + (alpha * alpha - 1.0) * u2));
    double sinTheta = sqrt(1.0 - cosTheta * cosTheta);

    // Compute half vector in tangent space
    double x = sinTheta * cos(phi);
    double y = sinTheta * sin(phi);
    double z = cosTheta;

    // Orthonormal basis
    Vec3 tangent, bitangent;
    orthonormalBasis(normal, tangent, bitangent);

    // Transform h to world space
    Vec3 h = (tangent * x + bitangent * y + normal * z).normalized();

    // Compute outgoing direction
    Vec3 outDir = (2.0 * viewDir.dot(h) * h - viewDir).normalized();

    // Ensure the outgoing direction is in the same hemisphere as the normal
    if (outDir.dot(normal) <= 0.0) {
        outDir = -outDir;
    }

    return outDir;
}

double Material::ggxNormalDistribution(const Vec3& h, const Vec3& n, double roughness) const {
    double alpha = roughness * roughness;
    double alpha2 = alpha * alpha;
    double NdotH = std::max(n.dot(h), 0.0);
    double NdotH2 = NdotH * NdotH;

    double denom = NdotH2 * (alpha2 - 1.0) + 1.0;
    denom = M_PI * denom * denom;

    return alpha2 / denom;
}

double Material::ggxGeometryFunction(const Vec3& v, const Vec3& l, const Vec3& n, double roughness) const {
    double NdotV = std::max(n.dot(v), 0.0);
    double NdotL = std::max(n.dot(l), 0.0);
    double r = roughness + 1.0;
    double k = (r * r) / 8.0;

    double G1V = NdotV / (NdotV * (1.0 - k) + k);
    double G1L = NdotL / (NdotL * (1.0 - k) + k);

    return G1V * G1L;
}

Vec3 Material::evaluateBRDF(const Vec3& normal, const Vec3& viewDir, const Vec3& lightDir, const Vec3& diffuseColor) const {
    Vec3 halfVector = (viewDir + lightDir).normalized();
    double NdotL = std::max(0.0, normal.dot(lightDir));
    double NdotV = std::max(0.0, normal.dot(viewDir));
    double NdotH = std::max(0.0, normal.dot(halfVector));
    double VdotH = std::max(0.0, viewDir.dot(halfVector));

    // Diffuse component
    Vec3 diffuse = diffuseColor * kd / M_PI;

    // Specular component (using GGX)
    double D = ggxNormalDistribution(halfVector, normal, roughness);
    double G = ggxGeometryFunction(viewDir, lightDir, normal, roughness);
    Vec3 F = fresnelSchlick(specularColor, VdotH);
    Vec3 specular = (D * G * F) / (4.0 * NdotL * NdotV + 1e-6);

    return diffuse + specular;
}

