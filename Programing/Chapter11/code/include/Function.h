#ifndef DYNAMIC_FUNCTION_7_H
#define DYNAMIC_FUNCTION_7_H

#include <exprtk.hpp>

class Function{
public:
    virtual double operator()(double x_val) const = 0;
    virtual ~Function() {}
};

class phi : public Function {
    public:
        double operator()(double x) const override {
            if(9.0/20 <= x && x < 1.0/2)
                return 20*(x - 9.0/20);
            else if(0.5 <= x && x <= 11.0/20)
                return -20*(x - 11.0/20);
            else 
                return 0;
        }
    };

#endif // DYNAMIC_FUNCTION_7_H
