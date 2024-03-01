#include <iomanip>
#include "Color.h"

Color::Color() {
    this->r = 0;
    this->g = 0;
    this->b = 0;
}

Color::Color(double r, double g, double b)
{
    this->r = r;
    this->g = g;
    this->b = b;
}

Color::Color(const Color &other)
{
    this->r = other.r;
    this->g = other.g;
    this->b = other.b;
}

std::ostream &operator<<(std::ostream &os, const Color &c)
{
    os << std::fixed << std::setprecision(0) << "rgb(" << c.r << ", " << c.g << ", " << c.b << ")";
    return os;
}

Color operator+(const Color &c1, const Color &c2)
{
    Color result;
    result.r = c1.r + c2.r;
    result.g = c1.g + c2.g;
    result.b = c1.b + c2.b;
    return result;
}

Color operator-(const Color &c1, const Color &c2)
{
    Color result;
    result.r = c1.r - c2.r;
    result.g = c1.g - c2.g;
    result.b = c1.b - c2.b;
    return result;
}

Color operator/(const Color &c1, int n)
{
    Color result;
    result.r = c1.r / n;
    result.g = c1.g / n;
    result.b = c1.b / n;
    return result;
}

Color operator*(double n, const Color &c1)
{
    Color result;
    result.r = c1.r * n;
    result.g = c1.g * n;
    result.b = c1.b * n;
    return result;
}

