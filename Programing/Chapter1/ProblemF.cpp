#include "Function.hpp"
#include "EquationSolver.hpp"
#include "ProblemA.hpp"
#include <iostream>
#include <cmath>


class F1 : public Function {
public:
    double operator()(double x) const {
        int l = 89, h=49, D=55;
        double beta1=11.5*Pi/180.;
        x = l*sin(beta1)*sin(x)*cos(x)+l*cos(beta1)*pow(sin(x),2)
            -cos(x)*(sin(beta1)*(h+0.5*D)-0.5*D*tan(beta1)) 
            -sin(x)*(cos(beta1)*(h+0.5*D)-0.5*D);
        return x;
    }
};

class F2 : public Function {
public:
    double operator()(double x) const {
        int l = 89, h=49, D=30;
        double beta1=11.5*Pi/180.;
        x = l*sin(beta1)*sin(x)*cos(x)+l*cos(beta1)*pow(sin(x),2)
            -cos(x)*(sin(beta1)*(h+0.5*D)-0.5*D*tan(beta1)) 
            -sin(x)*(cos(beta1)*(h+0.5*D)-0.5*D);
        return x;
    }
};

int main() {
    std::cout << "Attention! The output value of this method is the angle in radians " << std::endl;
    std::cout << "Please choose the problem to solve:" << std::endl;
    int n;

    Function *selected_function = nullptr;
    

    std::cin >> n;
    switch (n) {
        case 1:
            selected_function = new F1();
            solve_f4(selected_function);
            break;
        case 2:
            selected_function = new F2();
            solve_f4(selected_function);
            break;
        case 3:
            selected_function = new F1();
            solve_f3(selected_function);
            break;
        default:
            std::cout << "Invalid choice. Exiting." << std::endl;
            break;
    }
    delete selected_function;
    return 0;
}