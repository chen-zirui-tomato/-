#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>

using namespace std;
const double Pi = acos(-1.);
const double e = std::exp(1);

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

void solve_f1(double eps = 1e-6,double delta = 1e-6, int max_iter = 50) {
    double a, b;
    int n;

    std::cout << "Solving with Bisection Method" << std::endl;
    std::cout << "Select the number of test: ";
    std::cin >> n;

    std::cout<<"Enter the lower bound of the interval: "<< std::endl;
    std::cin>>a;
    std::cout<<"Enter the upper bound of the interval: "<< std::endl;
    std::cin>>b;
    std::cout<<"Enter the default value of eps,delta,max_iter ";
    std::cout<<"(or press enter to use default values)"<< std::endl;
    std::cin.ignore();
    char c = std::cin.get(); 

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
            return;
    }

    if (c == '\n'){
        Bisection_Method solver_f1(*selected_function, a, b, eps, delta, max_iter);
        double x = solver_f1.solve();
        std::cout << "A root with default values is: " << x << std::endl;
    }
    else{
        std::cin>>eps>>delta>>max_iter;
        Bisection_Method solver_f1(*selected_function, a, b, eps, delta, max_iter);
        double x = solver_f1.solve();
        std::cout << "A root with custom values is: " << x << std::endl;
    }
    return;
}

int main() {

    solve_f1();

    return 0;
}
