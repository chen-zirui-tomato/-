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
        return x-tan(x);
    }
};

void solve_f2(double eps = 1e-6, int max_iter = 50) {
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

int main() {

    solve_f2();

    return 0;
}