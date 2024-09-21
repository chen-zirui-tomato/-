#ifndef FUNCTION
#define FUNCTION

class Function{
public:
    virtual double operator()(double x) const = 0;
    virtual double derivative(double x) const {
        double h = 1e-6;
        return (operator()(x+h) - operator()(x-h))/(2*h);
    }
};

#endif