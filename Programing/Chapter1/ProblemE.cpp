#include "Function.hpp"
#include "EquationSolver.hpp"
#include "ProblemA.hpp"
#include <iostream>
#include <cmath>


class F1 : public Function {
public:
    double operator()(double x) const {
        return -1.24+0.5*Pi-asin(x)-x*pow(1-x*x,0.5);
    }
};

int main() {

    std::cout<<"Enter the number of the funcion to use (1,2,3): "<< std::endl;
    int n;

    std::cin>>n;

    Function *selected_function = nullptr;
    selected_function = new F1();

    switch(n){
        case 1:
            solve_f1(selected_function ,1e-2);
            break;
        case 2:
            solve_f2(selected_function ,1e-2);
            break;
        case 3:
            solve_f3(selected_function ,1e-2);
            break;
        default:
            std::cout<<"Invalid input"<<std::endl;
            break;
    } 
    delete selected_function;
    return 0;
}