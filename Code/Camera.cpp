// Camera.cpp
#include "Camera.h"
#include <cmath>
#include <random>

// Default constructor
Camera::Camera()
    : position(0.0, 0.0, 0.0),
      target(0.0, 0.0, -1.0),
      up(0.0, 1.0, 0.0),
      fov(90.0),
      aspectRatio(16.0 / 9.0),
      imageWidth(400),
      imageHeight(225) {
    initialize();
}

// Parameterized constructor
Camera::Camera(const Vec3& position_, const Vec3& target_, const Vec3& up_,
               double fov_, double aspectRatio_,
               double aperture, double focusDist_)
    : position(position_), target(target_), up(up_), fov(fov_), aspectRatio(aspectRatio_), lensRadius(aperture / 2.0), focusDist(focusDist_) {
    initialize();
}

// Set image dimensions
void Camera::setImageDimensions(int width, int height) {
    imageWidth = width;
    imageHeight = height;
    aspectRatio = static_cast<double>(imageWidth) / imageHeight;
    initialize(); // Recalculate camera parameters based on new aspect ratio
}

void Camera::initialize() {
    // Calculate the camera coordinate system
    w = (position - target).normalized();          // Camera backward vector
    u = up.cross(w).normalized();                  // Camera right vector
    v = w.cross(u);                                // Camera up vector

    // Swap the cross product operands for 'u'
    u = w.cross(up).normalized();                       // Corrected Camera right vector

    // Convert field of view from degrees to radians
    double theta = fov * M_PI / 180.0;
    double viewportHeight = 2.0 * std::tan(theta / 2);
    double viewportWidth = aspectRatio * viewportHeight;

    // Calculate the horizontal and vertical vectors scaled by the viewport size
    horizontal = u * viewportWidth * focusDist;
    vertical = v * viewportHeight * focusDist;

    // Calculate the lower-left corner of the image plane
    lowerLeftCorner = position - horizontal / 2.0 - vertical / 2.0 - w * focusDist;

}


// Generate a ray through the pixel at normalized coordinates (u, v)
Ray Camera::getRay(double s, double t) const {
    // Random point in lens aperture
    Vec3 rd = lensRadius * randomInUnitDisk();
    Vec3 offset = u * rd.x + v * rd.y;

    Vec3 direction = lowerLeftCorner + horizontal * s + vertical * t - position - offset;

    return Ray(position + offset, direction.normalized());
}

Vec3 Camera::randomInUnitDisk() const {
    static std::mt19937 generator(std::random_device{}());
    static std::uniform_real_distribution<double> distribution(-1.0, 1.0);
    Vec3 p;
    do {
        p = Vec3(distribution(generator), distribution(generator), 0);
    } while (p.lengthSquared() >= 1.0);
    return p;
}

