#ifndef FUNCTION
#define FUNCTION

enum class SplineType{
    Natural,
    Complete,
    Periodic
};

class Function{
public:
    virtual double operator()(double x) const = 0;
    virtual double derivative(double x) const {
        double h = 1e-6;
        return (operator()(x+h) - operator()(x-h))/(2*h);
    }
    virtual ~Function(){}

    double bisection_method(double a, double b, double tolerance = 1e-6, int max_iterations = 100) {
        if (operator()(a) * operator()(b) >= 0) {
            std::cerr << "函数在区间端点的值有相同符号，无法应用二分法。" << std::endl;
            return -1;
        }

        double c;
        for (int i = 0; i < max_iterations; ++i) {
            // 计算中点
            c = (a + b) / 2;
            if (std::abs(operator()(c)) < tolerance) {
                return c; // 找到根
            }
            // 更新区间
            if (operator()(c) * operator()(a) < 0) {
                b = c; // 根在 [a, c] 区间
            } else {
                a = c; // 根在 [c, b] 区间
            }
        }

        return c; // 如果未找到根，返回最后的估计
    }

};

#endif