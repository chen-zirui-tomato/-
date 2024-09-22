#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>

const double Pi = acos(-1.);
const double e = std::exp(1);

void solve_f1(Function *func ,double eps = 1e-7,double delta = 1e-6, int max_iter = 50) {
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

    //选择函数

    if (c == '\n'){
        Bisection_Method bisectionsolver(*func, a, b, eps, delta, max_iter);
        double x = bisectionsolver.solve();
        std::cout << "A root with default values is: " << x << std::endl;
    }
    else{
        std::cin>>eps>>delta>>max_iter;
        Bisection_Method bisectionsolver(*func, a, b, eps, delta, max_iter);
        double x = bisectionsolver.solve();
        std::cout << "A root with custom values is: " << x << std::endl;
    }
    return;
}

void solve_f2(Function *func ,double eps = 1e-7, int max_iter = 50) {
    double x0;
    std::cout << "Solving with Newton Method" << std::endl;
    std::cout<<"Enter the initial value of x0: "<< std::endl;
    std::cin>>x0;
    std::cout<<"Enter the default value of eps,max_iter ";
    std::cout<<"(or press enter to use default values)"<< std::endl;
    std::cin.ignore();
    char c = std::cin.get(); 

    if (c == '\n'){
        Newton_Method newtonsolver(*func, x0, eps, max_iter);
        double x = newtonsolver.solve();
        std::cout << "A root with default values is: " << x << std::endl;
    }
    else{
        std::cin>>eps>>max_iter;
        Newton_Method newtonsolver(*func, x0, eps, max_iter);
        double x = newtonsolver.solve();
        std::cout << "A root with custom values is: " << x << std::endl;
    }
    return ;
}

void solve_f3(Function *func ,double eps = 1e-7, double delta = 1e-6, int max_iter = 50) {
    double a, b;
    std::cout << "Solving with Secant Method" << std::endl;
    std::cout<<"Enter the value of x0, x1:"<< std::endl;
    std::cin>>a>>b;
    std::cout<<"Enter the default value of eps,delta,max_iter ";
    std::cout<<"(or press enter to use default values)"<< std::endl;
    std::cin.ignore();
    char c = std::cin.get(); 

    if (c == '\n'){
        Secant_Method secantsolver(*func, a, b, eps, delta, max_iter);
        double x = secantsolver.solve();
        std::cout << "A root with default values is: " << x << std::endl;
    }
    else{
        std::cin>>eps>>delta>>max_iter;
        Secant_Method secantsolver(*func, a, b, eps, delta, max_iter);
        double x = secantsolver.solve();
        std::cout << "A root with custom values is: " << x << std::endl;
    }
    return ;
}

//计算和角度
void solve_f4(Function *func ,double eps = 1e-7, int max_iter = 50) {
    double x0;
    std::cout << "Solving with Newton Method" << std::endl;
    std::cout<<"Enter the initial value of x0: "<< std::endl;
    std::cin>>x0;
    x0 = x0*Pi/180.;
    std::cout<<"Enter the default value of eps,max_iter ";
    std::cout<<"(or press enter to use default values)"<< std::endl;
    std::cin.ignore();
    char c = std::cin.get(); 

    if (c == '\n'){
        Newton_Method newtonsolver(*func, x0, eps, max_iter);
        double x = newtonsolver.solve();
        std::cout << "A root with default values is: " << x << std::endl;
    }
    else{
        std::cin>>eps>>max_iter;
        Newton_Method newtonsolver(*func, x0, eps, max_iter);
        double x = newtonsolver.solve();
        std::cout << "A root with custom values is: " << x << std::endl;
    }
    return ;
}
