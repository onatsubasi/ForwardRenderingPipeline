#ifndef __COLOR_H__
#define __COLOR_H__

class Color
{
public:
    double r, g, b;

    Color();
    Color(double r, double g, double b);
    Color(const Color &other);
    friend std::ostream &operator<<(std::ostream &os, const Color &c);
    friend Color operator+(const Color &c1, const Color &c2);
    friend Color operator-(const Color &c1, const Color &c2);
    friend Color operator/(const Color &c1, int n);
    friend Color operator*(double n, const Color &c1);
};

#endif