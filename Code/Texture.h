#ifndef TEXTURE_H
#define TEXTURE_H

#include "Vec3.h"
#include <string>
#include <vector>

class Texture {
public:
    Texture(const std::string& filename);

    Vec3 getColor(double u, double v) const;

private:
    int width, height;
    std::vector<Vec3> data;
};

#endif // TEXTURE_H
