#include "Texture.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>

Texture::Texture(const std::string& filename) {
    std::ifstream ppmFile(filename, std::ios::binary);
    if (!ppmFile.is_open()) {
        std::cerr << "Failed to open texture file: " << filename << std::endl;
        width = height = 0;
        return;
    }

    std::string header;
    ppmFile >> header;
    if (header != "P6") {
        std::cerr << "Invalid PPM file (must be binary P6): " << filename << std::endl;
        width = height = 0;
        ppmFile.close();
        return;
    }

    // Read width, height, and max color value
    int maxColorValue;
    // Skip comments and read width and height
    while (true) {
        char c = ppmFile.peek();
        if (c == '#') {
            ppmFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        } else if (std::isspace(c)) {
            ppmFile.get();
        } else {
            break;
        }
    }

    ppmFile >> width >> height;

    ppmFile >> maxColorValue;
    ppmFile.get(); // Consume the newline character after max color value

    if (width <= 0 || height <= 0 || maxColorValue != 255) {
        std::cerr << "Invalid PPM file parameters: " << filename << std::endl;
        width = height = 0;
        ppmFile.close();
        return;
    }

    data.resize(width * height);

    // Read pixel data
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            unsigned char rgb[3];
            ppmFile.read(reinterpret_cast<char*>(rgb), 3);
            if (ppmFile.eof()) {
                std::cerr << "Unexpected end of file in PPM file: " << filename << std::endl;
                width = height = 0;
                ppmFile.close();
                return;
            }
            int index = j * width + i;
            data[index] = Vec3(rgb[0], rgb[1], rgb[2]) / 255.0;
        }
    }

    ppmFile.close();
}

Vec3 Texture::getColor(double u, double v) const {
    if (width == 0 || height == 0) {
        return Vec3(0.0, 0.0, 0.0); // Return black if texture is invalid
    }

    // Wrap texture coordinates
    u = u - std::floor(u);
    v = v - std::floor(v);

    int i = static_cast<int>(u * width);
    int j = static_cast<int>((1 - v) * height); // Flip v to match image coordinates

    // Clamp indices
    i = std::min(std::max(i, 0), width - 1);
    j = std::min(std::max(j, 0), height - 1);

    return data[j * width + i];
}
