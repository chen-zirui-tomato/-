#include "Function.hpp"
#include "EquationSolver.hpp"
#include "ProblemA.hpp"
#include <iostream>
#include <cmath>

using namespace std;

class F1 : public Function {
public:
    double operator()(double x) const {
        return sin(x/2)-1;
    }
};

class F2 : public Function {
public:
    double operator()(double x) const {
        return pow(e,2)-tan(x);
    }
};

class F3 : public Function {
public:
    double operator()(double x) const {
        return pow(x,3)-12*pow(x,2)+3*x+1;
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
        default:
            std::cout << "Invalid choice. Exiting." << std::endl;
            break;;
    }
    
    solve_f3(selected_function);
    delete selected_function;

    return 0;
}
