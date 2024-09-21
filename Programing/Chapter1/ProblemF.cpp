#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>

const double Pi = acos(-1.);

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

void solve_f2(int n,double eps = 1e-7, int max_iter = 50) {
    double x0;
    std::cout << "Solving with Newton Method" << std::endl;
    std::cout<<"Enter the initial value of x0: "<< std::endl;
    std::cin>>x0;
    x0 = x0*Pi/180.;
    std::cout<<"Enter the default value of eps,max_iter ";
    std::cout<<"(or press enter to use default values)"<< std::endl;
    std::cin.ignore();
    char c = std::cin.get(); 

    Function *selected_function = nullptr;
    switch (n) {
        case 1:
            selected_function = new F1();
            break;
        case 2:
            selected_function = new F2();
            break;
        default:
            std::cout << "Invalid choice. Exiting." << std::endl;
            return;
    }

    if (c == '\n'){
        Newton_Method solver_f1(*selected_function, x0, eps, max_iter);
        double x = solver_f1.solve();
        std::cout << "A root with default values is: " << x/Pi*180. << std::endl;
    }
    else{
        std::cin>>eps>>max_iter;
        Newton_Method solver_f1(*selected_function, x0, eps, max_iter);
        double x = solver_f1.solve();
        std::cout << "A root with custom values is: " << x/Pi*180. << std::endl;
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
        std::cout << "A root with default values is: " << x/Pi*180. << std::endl;
    }
    else{
        std::cin>>eps>>delta>>max_iter;
        Secant_Method solver_f1(F1(), a, b, eps, delta, max_iter);
        double x = solver_f1.solve();
        std::cout << "A root with custom values is: " << x/Pi*180. << std::endl;
    }
    return ;
}

int main() {
    std::cout << "Please choose the problem to solve:" << std::endl;
    int n;
    std::cin >> n;
    switch (n) {
        case 1:
            solve_f2(n);
            break;
        case 2:
            solve_f2(n);
            break;
        case 3:
            solve_f3();
            break;
        default:
            std::cout << "Invalid choice. Exiting." << std::endl;
    }
    return 0;
}