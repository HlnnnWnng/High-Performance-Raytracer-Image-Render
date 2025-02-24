#ifndef COLOR_H
#define COLOR_H

#include "Vec3.h"
#include <iostream>

// Write the color to an output stream as three integers in [0,255]
void writeColor(std::ostream& out, const Vec3& pixelColor);

#endif // COLOR_H
