#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <string>

// Forward declarations of classes
class Camera;
class Vec3;
class Ray;

// Function to write a PPM image file
void writePPM(const std::string& filename, int width, int height, const std::string& pixelData);

#endif // RAYTRACER_H
