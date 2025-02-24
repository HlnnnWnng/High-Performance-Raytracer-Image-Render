// Camera.h
#ifndef CAMERA_H
#define CAMERA_H

#include "Vec3.h"
#include "Ray.h"

class Camera {
public:
    // Camera parameters
    Vec3 position;     // Camera position in world coordinates
    Vec3 target;       // Look-at point
    Vec3 up;           // Up direction

    double fov;        // Field of view in degrees
    double aspectRatio; // Aspect ratio of the image (width / height)

    // Image plane vectors
    Vec3 horizontal;       // Horizontal axis of the camera
    Vec3 vertical;         // Vertical axis of the camera
    Vec3 lowerLeftCorner;  // Lower-left corner of the image plane

    double exposure;

    // Image dimensions
    int imageWidth;   // Width of the image in pixels
    int imageHeight;  // Height of the image in pixels

    // **New parameters for depth of field**
    double lensRadius;    // Radius of the lens aperture
    double focusDist;     // Focus distance

    // Constructors
    Camera();
    Camera(const Vec3& position_, const Vec3& target_, const Vec3& up_,
               double fov_, double aspectRatio_,
               double aperture, double focusDist_);

    // Set image dimensions
    void setImageDimensions(int width, int height);

    // Generate a ray going through pixel (u, v) on the image plane
    Ray getRay(double u, double v) const;

private:
    // Helper method to initialize camera parameters
    void initialize();
    // Camera coordinate system basis vectors
    Vec3 u, v, w;
    // Random point in unit disk for lens sampling
    Vec3 randomInUnitDisk() const;
};

#endif // CAMERA_H
