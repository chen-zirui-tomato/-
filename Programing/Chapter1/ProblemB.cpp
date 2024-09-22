#include "Function.hpp"
#include "EquationSolver.hpp"
#include "ProblemA.hpp"
#include <iostream>
#include <cmath>

using namespace std;

class F1 : public Function {
public:
    double operator()(double x) const {
        return 1/x-tan(x);
    }
};

class F2 : public Function {
public:
    double operator()(double x) const {
        return 1/x-std::pow(2,x);
    }
};

class F3 : public Function {
public:
    double operator()(double x) const {
        return std::pow(2,-x)+std::pow(e,x)+2*cos(x)-6;
    }
};

class F4 : public Function {
public:
    double operator()(double x) const {
        double r1 = (pow(x,3)+4*pow(x,2)+3*x+5);
        double r2 = (2*pow(x,3)-9*pow(x,2)+18*x-2);
        return r1/r2;
    }
};


int main() {

    int n;
    std::cout << "Select the number of test: ";
    std::cin >> n;
    
    //选择函数
    Function *selected_function = nullptr;
    switch (n) {
        case 1:
            selected_function = new F1();
            break;
        case 2:
            selected_function = new F2();
            break;
        case 3:
            selected_function = new F3();
            break;
        case 4:
            selected_function = new F4();
            break;
        default:
            std::cout << "Invalid choice. Exiting." << std::endl;
            break;
    }
    
    solve_f1(selected_function);
    delete selected_function;

    return 0;
}
