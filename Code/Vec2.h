// Vec2.h
#ifndef VEC2_H
#define VEC2_H

class Vec2 {
public:
    double x, y;

    // Constructors
    Vec2() : x(0), y(0) {}
    Vec2(double x_, double y_) : x(x_), y(y_) {}

    // Operator overloading
    Vec2 operator+(const Vec2& v) const { return Vec2(x + v.x, y + v.y); }
    Vec2 operator-(const Vec2& v) const { return Vec2(x - v.x, y - v.y); }
    Vec2 operator*(double t) const { return Vec2(x * t, y * t); }
    Vec2 operator/(double t) const { return Vec2(x / t, y / t); }

    Vec2& operator+=(const Vec2& v) { x += v.x; y += v.y; return *this; }
    Vec2& operator-=(const Vec2& v) { x -= v.x; y -= v.y; return *this; }
    Vec2& operator*=(double t) { x *= t; y *= t; return *this; }
    Vec2& operator/=(double t) { x /= t; y /= t; return *this; }

    // Dot product
    double dot(const Vec2& v) const { return x * v.x + y * v.y; }

    // Length and normalization
    double length() const { return std::sqrt(x * x + y * y); }
    Vec2 normalized() const { return *this / length(); }

    // Scalar multiplication from the left
    friend Vec2 operator*(double t, const Vec2& v) { return Vec2(t * v.x, t * v.y); }
};

#endif // VEC2_H
