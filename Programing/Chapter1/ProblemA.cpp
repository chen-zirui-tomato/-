#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>

const double Pi = acos(-1.);

class F1 : public Function {
public:
    double operator()(double x) const {
        return x * x * x - 3 * x * x + 2 * x - 5;
    }
};

void solve_f1(double eps = 1e-6,double delta = 1e-6, int max_iter = 50) {
    double a, b;
    std::cout << "Solving with Bisection Method" << std::endl;
    std::cout<<"Enter the lower bound of the interval: "<< std::endl;
    std::cin>>a;
    std::cout<<"Enter the upper bound of the interval: "<< std::endl;
    std::cin>>b;
    std::cout<<"Enter the default value of eps,delta,max_iter ";
    std::cout<<"(or press enter to use default values)"<< std::endl;
    std::cin.ignore();
    char c = std::cin.get(); 

    if (c == '\n'){
        Bisection_Method solver_f1(F1(), a, b, eps, delta, max_iter);
        double x = solver_f1.solve();
        std::cout << "A root with default values is: " << x << std::endl;
    }
    else{
        std::cin>>eps>>delta>>max_iter;
        Bisection_Method solver_f1(F1(), a, b, eps, delta, max_iter);
        double x = solver_f1.solve();
        std::cout << "A root with custom values is: " << x << std::endl;
    }
    return;
}

void solve_f2(double eps = 1e-7, int max_iter = 50) {
    double x0;
    std::cout << "Solving with Newton Method" << std::endl;
    std::cout<<"Enter the initial value of x0: "<< std::endl;
    std::cin>>x0;
    std::cout<<"Enter the default value of eps,max_iter ";
    std::cout<<"(or press enter to use default values)"<< std::endl;
    std::cin.ignore();
    char c = std::cin.get(); 

    if (c == '\n'){
        Newton_Method solver_f1(F1(), x0, eps, max_iter);
        double x = solver_f1.solve();
        std::cout << "A root with default values is: " << x << std::endl;
    }
    else{
        std::cin>>eps>>max_iter;
        Newton_Method solver_f1(F1(), x0, eps, max_iter);
        double x = solver_f1.solve();
        std::cout << "A root with custom values is: " << x << std::endl;
    }
    return ;
}

void solve_f3(double eps = 1e-7, double delta = 1e-6, int max_iter = 50) {
    double a, b;
    std::cout << "Solving with Secant Method" << std::endl;
    std::cout<<"Enter the value of x0, x1:"<< std::endl;
    std::cin>>a>>b;
    std::cout<<"Enter the default value of eps,delta,max_iter ";
    std::cout<<"(or press enter to use default values)"<< std::endl;
    std::cin.ignore();
    char c = std::cin.get(); 

    if (c == '\n'){
        Secant_Method solver_f1(F1(), a, b, eps, delta, max_iter);
        double x = solver_f1.solve();
        std::cout << "A root with default values is: " << x << std::endl;
    }
    else{
        std::cin>>eps>>delta>>max_iter;
        Secant_Method solver_f1(F1(), a, b, eps, delta, max_iter);
        double x = solver_f1.solve();
        std::cout << "A root with custom values is: " << x << std::endl;
    }
    return ;
}

int main() {

    // solve_f1();
    //solve_f2();
    solve_f3();

    return 0;
}