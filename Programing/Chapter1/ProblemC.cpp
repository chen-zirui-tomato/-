#include "Function.hpp"
#include "EquationSolver.hpp"
#include "ProblemA.hpp"
#include <iostream>
#include <cmath>

using namespace std;

class F1 : public Function {
public:
    double operator()(double x) const {
        return x-tan(x);
    }
};

int main() {
    Function *selected_function = nullptr;
    selected_function = new F1();

    solve_f2(selected_function);
    delete selected_function;
    
    return 0;
}