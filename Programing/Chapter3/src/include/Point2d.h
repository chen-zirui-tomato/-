#include <cmath>

class Point2d {
public:
    double x, y;
    Point2d(double x=0, double y=0) : x(x), y(y) {}
    double distance(const Point2d& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        return sqrt(dx*dx + dy*dy);
    }
    Point2d operator+(const Point2d& other) const {
        return Point2d(x + other.x, y + other.y);
    }
};