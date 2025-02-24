# Makefile
# Compiler (Use Homebrew's g++-14)
CXX = g++-14

# Compiler Flags
CXXFLAGS = -Wall -Wextra -std=c++17 -O3


# Target executable
TARGET = raytracer

# Source files
SRCS = raytracer.cpp Camera.cpp Vec3.cpp color.cpp Sphere.cpp Triangle.cpp Cylinder.cpp Scene.cpp Material.cpp PointLight.cpp Texture.cpp Mesh.cpp AABB.cpp BVHNode.cpp MeshBVHNode.cpp RectangularAreaLight.cpp TriangleAreaLight.cpp
# Object files
OBJS = $(SRCS:.cpp=.o)

# Header files
HEADERS = raytracer.h Camera.h Vec3.h Ray.h color.h Shape.h Sphere.h Triangle.h Cylinder.h Scene.h Material.h PointLight.h Texture.h Mesh.h AABB.h BVHNode.h MeshBVHNode.h json.hpp TriangleAreaLight.h

# Default rule
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Compile source files into object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)
