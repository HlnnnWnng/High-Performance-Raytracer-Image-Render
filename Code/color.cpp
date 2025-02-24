#include "color.h"

// Write the color to an output stream as three integers in [0,255]
void writeColor(std::ostream& out, const Vec3& pixelColor){
    // Clamp color components
    double r = std::min(std::max(pixelColor.x, 0.0), 1.0);
    double g = std::min(std::max(pixelColor.y, 0.0), 1.0);
    double b = std::min(std::max(pixelColor.z, 0.0), 1.0);

    // Convert to [0,255] range
    int ir = static_cast<int>(255.999 * r);
    int ig = static_cast<int>(255.999 * g);
    int ib = static_cast<int>(255.999 * b);

    out << ir << ' ' << ig << ' ' << ib << '\n';
}
