#ifndef EQUATIONSOLVER
#define EQUATIONSOLVER

#include <iostream>
#include "Function.hpp"
#include <cmath>
#include <limits>

class EquationSolver{
protected:
    const Function& F;
public:
    EquationSolver(const Function& F) : F(F) {}
    virtual double solve() = 0;
};

class Bisection_Method : public EquationSolver{
private:
    double a,b;
    double eps,delta;
    int max_iter;
public:
    Bisection_Method(const Function& F, double a, double b, 
    double eps, double delta ,int max_iter) : 
    EquationSolver(F), a(a), b(b), eps(eps), delta(delta), max_iter(max_iter) {}

    double solve() override{
        std::cout << "使用二分法求解方程..." << std::endl;
        double fa = F(a);
        double fb = F(b);
        double h = b-a;

        //检查输入区间是否合法
        if(fa * fb > 0){
            std::cout << "错误!f(a)与f(b)应有不同的符号" << std::endl;
            return std::numeric_limits<double>::quiet_NaN();
        }
        for(int i=0; i < max_iter; i++){
            h = h/2;
            double c = a + h;
            double fc = F(c);

            //检查精度
            if(fabs(fc) < eps || h < delta ){
                std::cout << "方程的解:x = " << c << std::endl;
                return c ;
            }

            //更新区间
            else if(fa * fc > 0){
                a = c;
            }    
        }
    std::cout << "达到最大迭代次数，未找到精确解。"<<std::endl;
    std::cout << "可尝试使用更精确的初值"<< std::endl;
    return std::numeric_limits<double>::quiet_NaN();
    }
};

class Newton_Method : public EquationSolver{
private:
    double x0;
    double eps;
    int max_iter;
public:
    Newton_Method(const Function &F, double x0,
        double eps, double max_iter) :
        EquationSolver(F), x0(x0), eps(eps), max_iter(max_iter) {}

    double solve() override{
        std::cout << "使用牛顿法求解方程..." << std::endl;
        double dfx = F.derivative(x0);
        if(fabs(dfx) < eps){
            std::cout << "导数接近0,无法找到精确解。" << std::endl;
        }
        for(int i = 0; i < max_iter; i++){
            double fx = F(x0);
            if(fabs(fx) < eps){
                std::cout << "方程的解: x = " << x0 << std::endl;
                return x0;
            }
            else{
                x0 = x0 - fx/dfx;
                dfx = F.derivative(x0);
            }
        }
        std::cout << "达到最大迭代次数，未找到精确解。"<<std::endl;
        std::cout << "可尝试使用更精确的初值"<< std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
};

class Secant_Method : public EquationSolver{
private:
    double a,b;
    double eps,delta;
    int max_iter;
public:
    Secant_Method(const Function& F, double a, double b, 
    double eps , double delta , int max_iter) : 
    EquationSolver(F), a(a), b(b), eps(eps), delta(delta), max_iter(max_iter) {}

    double solve() override{
        std::cout << "使用割线法求解方程..." << std::endl;
        double fa = F(a);
        double fb = F(b);
        for(int i = 0;i < max_iter;i++){
                double s=(b-a)/(fb-fa);
                a = b,fa = fb;
                b = b-s*fb,fb = F(b);
                if(fabs(b-a) < delta || fabs(fb) < eps){
                    std::cout << "找到的解: x = " << b << std::endl;
                    return b;
                }
            }
            std::cout << "达到最大迭代次数，未找到精确解" << std::endl;
            std::cout << "可尝试使用更精确的初值"<< std::endl;
            return std::numeric_limits<double>::quiet_NaN();
        }
};

#endif