# High Performance Raytracer Image Render

This repository contains the source code and supporting files for a raytracer project. The project implements a full-featured raytracer that supports a wide range of rendering techniques—from basic image generation and intersection tests to advanced features like BRDF sampling, depth-of-field effects, and subsurface scattering.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Installation and Compilation](#installation-and-compilation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Dependencies](#dependencies)
- [References](#references)

## Overview
This raytracer project is designed to render 3D scenes by simulating the interaction of light with objects. It reads scene configurations from JSON files and supports various geometric primitives (spheres, triangles, cylinders, and meshes). The renderer uses techniques including:
- **Ray-Object Intersection Tests** (with support for spheres, triangles using the Möller-Trumbore algorithm, and cylinders)
- **Shading Models** (Blinn-Phong, reflection, refraction, and advanced BRDF sampling)
- **Texture Mapping** (with UV coordinate generation for different primitives)
- **Acceleration Structures** (Bounding Volume Hierarchy - BVH for efficient intersection testing)
- **Advanced Sampling Techniques** (pixel sampling, lens sampling for depth-of-field, and light sampling for soft shadows)
- **Additional Effects** such as tone mapping, transmission filters, and BSSRDF-based subsurface scattering.

## Features
- **Image Writing:** Generates images in PPM format (currently in ASCII P3, with potential for binary P6 for performance improvements).
- **Camera and Ray Generation:** Implements a virtual pinhole camera using the Look-At method. Generates rays for each pixel based on the camera’s orientation, field of view, and image plane geometry.
- **Intersection Tests:** Provides intersection routines for:
  - **Spheres:** Solving quadratic equations for ray-sphere intersections.
  - **Triangles:** Using the Möller-Trumbore algorithm.
  - **Cylinders:** Including tests for lateral surfaces and caps.
- **Scene Management and JSON Parsing:** Uses the [nlohmann/json](https://github.com/nlohmann/json) library to load scene configurations, camera parameters, materials, and lighting information.
- **Shading and Lighting:** Implements the Blinn-Phong shading model (with ambient, diffuse, and specular components) along with shadow casting. Supports reflection and refraction via recursive ray tracing and uses Fresnel equations for blending.
- **Texture Mapping:** Maps textures onto objects by computing UV coordinates for spheres, cylinders, and triangles. Also supports loading meshes from `.obj` files.
- **Acceleration Structure (BVH):** Integrates a Bounding Volume Hierarchy to reduce the number of ray-object intersection tests, significantly boosting rendering performance for complex scenes.
- **Advanced Sampling Techniques:**
  - **Pixel Sampling:** Stratified sampling (e.g., 4x4 subpixels) to reduce aliasing.
  - **Lens Sampling:** Simulates depth-of-field effects by generating rays over a finite-sized aperture.
  - **BRDF Sampling:** Implements both diffuse (cosine-weighted hemisphere sampling) and specular (GGX microfacet model) sampling for physically based path tracing.
  - **Light Sampling:** Provides soft shadow effects by sampling area lights.
- **Additional Effects:**
  - **Tone Mapping:** Supports Reinhard and Exponential tone mapping (with Exponential tone mapping currently adopted).
  - **Transmission Filter and BSSRDF:** Models light transmission through transparent materials and simulates subsurface scattering for translucent surfaces.

## Installation and Compilation
1. **Clone the Repository:**
   ```bash
   git clone https://github.com/HlnnnWnng/High-Performance-Raytracer-Image-Render.git
   cd High-Performance-Raytracer-Image-Render
   ```

2. **Compile the Project:** A Makefile is provided. Simply run:
   ```bash
   make
   ```
   This will compile the source code and generate the executable (e.g., `raytracer`).

## Usage
Run the raytracer from the command line using the following format:
```bash
./raytracer <path-to-json> <pixel-samples> <output-file>
```
For example:
```bash
./raytracer Code/jsons/test.json 64 output.ppm
```
- **JSON File:** Specifies the scene configuration (camera, objects, lights, materials, textures).
- **Pixel Samples:** Number of samples per pixel for anti-aliasing.
- **Output File:** Name of the generated image in PPM format (which can be converted to PNG using ImageMagick).

## Project Structure
```
├── Code
│   ├── jsons         # JSON scene configuration files
│   ├── models        # .obj and .mtl files for mesh testing
│   └── [source files]# C++ source code for the raytracer
├── TestSuite         # Test cases and sample scenes with expected outputs
├── Makefile          # Build script to compile the project
└── README.md         # This file
```

## Dependencies
- **C++ Compiler:** C++11 or later (e.g., g++, clang++)
- **nlohmann/json:** For parsing JSON files ([GitHub Repository](https://github.com/nlohmann/json))
- Standard libraries such as `<iostream>`, `<fstream>`, `<cmath>`, etc.

## References
- [Möller-Trumbore Algorithm (MT97)](https://www.graphics.cornell.edu/pubs/1997/MT97.pdf)
- [Reinhard Tone Mapping (RSSF23)](https://www.cs.cmu.edu/~ph/865/2012/readings/Reinhard2002HDR.pdf)
- Additional details and theory can be found in the project report.
```
