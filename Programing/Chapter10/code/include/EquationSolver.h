#pragma once
#include"Function.h"
#include <iostream>
#include <cmath>
#include <exprtk.hpp>
#include <memory>

class EquationSolver {
public:
    virtual ~EquationSolver() = default;
    virtual void solve(const DynamicFunction& Func, double t_end, double h) = 0;

    virtual std::string name() const = 0;

    virtual void print_result(const std::vector<double>& u) const{
        std::cout << "Result from " << name() << " solver:\n";
        for(size_t i = 0; i < u.size(); ++i){
            std::cout << "u(" << i << ") = " << u[i] << std::endl;
        }
    };
protected:
    //四阶RK方法，用于求解初始解
    double RK4(const DynamicFunction& Func, double x, double y, double h){
        double k1 = Func(x, y);
        double k2 = Func(x + h/2, y + h*k1/2);
        double k3 = Func(x + h/2, y + h*k2/2);
        double k4 = Func(x + h, y + h*k3);
        return h/6*(k1 + 2*k2 + 2*k3 + k4);
    }
    //不动点迭代法，用于求解隐格式
    double fixedPointIteration(const std::function<double(double)>& G, double h,
                               double initial_guess, double tolerance = 1e-6, int max_iter = 100){
        double u_k = initial_guess;
        for(int i = 0; i < max_iter; ++i){
            double u_next = G(u_k);
            if(std::abs(u_next - u_k) < tolerance){
                return u_next;
            }
            u_k = u_next;
        }
        std::cerr << "fixedPointIteration: Failed to converge " << std::endl;
        return -1;
    }
};


class AdamsBashforthSolver : public EquationSolver {
private:
    int order_;
    //储存解向量
    std::vector<double> u_;
    //储存系数
    std::vector<double> coefficients_;
public:
    AdamsBashforthSolver(int order) : order_(order) {
        loadCoefficients(order);
    }
    void solve(const DynamicFunction& Func, double t_end, double h)override{
        //
        switch (order_) {
            case 1: integrateAB1(Func, t_end, h); break;
            case 2: integrateAB2(Func, t_end, h); break;
            case 3: integrateAB3(Func, t_end, h); break;
            case 4: integrateAB4(Func, t_end, h); break;
            default: throw std::invalid_argument("Invalid order for Adams-Bashforth");
        }
    }
    std::string name() const override{
        return "Adams-Bashforth";
    }
private:
    void loadCoefficients(int order) {
        // 预定义Adams-Bashforth系数表（阶数 → 系数）
        static const std::map<int, std::vector<double>> ab_coeffs = {
            {1, {1.0}},
            {2, {3.0/2, -1.0/2}},
            {3, {23.0/12, -16.0/12, 5.0/12}},
            {4, {55.0/24, -59.0/24, 37.0/24, -9.0/24}}
        };
        if (ab_coeffs.find(order) == ab_coeffs.end()) {
            throw std::invalid_argument("Invalid order for Adams-Bashforth");
        }
        coefficients_ = ab_coeffs.at(order);
    }

    void integrateAB1(const DynamicFunction& Func, double t_end, double h) {
        //
        int n = static_cast<int>(t_end / h);
        u_.resize(n + 1);
        u_[0] = Func(0, 0);
        for(int i = 1; i <= n; ++i){
            u_[i] = u_[i-1] + h*(coefficients_[0]*Func(u_[i-1], (i-1)*h));
        }
    }
    void integrateAB2(const DynamicFunction& Func, double t_end, double h) {
        //
        int n = static_cast<int>(t_end / h);
        u_.resize(n + 1);
        u_[0] = Func(0, 0);
        u_[1] = RK4(Func, 0, u_[0], h);
        for(int i = 2; i <= n; ++i){
            u_[i] = u_[i-1] + h*(coefficients_[0]*Func(u_[i-1], (i-1)*h) 
                            + coefficients_[1]*Func(u_[i-2], (i-2)*h));
        }
    }
    void integrateAB3(const DynamicFunction& Func, double t_end, double h) {
        //
        int n = static_cast<int>(t_end / h);
        u_.resize(n + 1);
        u_[0] = Func(0, 0);
        u_[1] = RK4(Func, 0, u_[0], h);
        u_[2] = RK4(Func, h, u_[1], h);
        for(int i = 3; i <= n; ++i){
            u_[i] = u_[i-1] + h*(coefficients_[0]*Func(u_[i-1], (i-1)*h) 
                            + coefficients_[1]*Func(u_[i-2], (i-2)*h)
                            + coefficients_[2]*Func(u_[i-3], (i-3)*h));      
        }
    }
    void integrateAB4(const DynamicFunction& Func, double t_end, double h) {
        //
        int n = static_cast<int>(t_end / h);
        u_.resize(n + 1);
        //u_[0]是给定的初始条件
        u_[0] = Func(0, 0);
        u_[1] = RK4(Func, 0, u_[0], h);
        u_[2] = RK4(Func, h, u_[1], h);
        u_[3] = RK4(Func, h+h, u_[2], h);
        for(int i = 4; i <= n; ++i){
            u_[i] = u_[i-1] + h*(coefficients_[0]*Func(u_[i-1], (i-1)*h) 
                            + coefficients_[1]*Func(u_[i-2], (i-2)*h)
                            + coefficients_[2]*Func(u_[i-3], (i-3)*h)
                            + coefficients_[3]*Func(u_[i-4], (i-4)*h));      
        }
    }
};


//=================================================================

class AdamsMoultonSolver : public EquationSolver {
private:
    int order_;
    //储存解向量
    std::vector<double> u_;
    //储存系数
    std::vector<double> coefficients_;
public:
    AdamsMoultonSolver(int order) : order_(order) {}
    void solve(const DynamicFunction& Func, double t_end, double h)override{
        //
    }
    std::string name() const override{
        return "Adams-Moulton";
    }
private: 
    void loadCoefficients(int order) {
        static const std::map<int, std::vector<double>> am_coeffs = {
            {1, {1.0}},
            {2, {1.0/2, 1.0/2}},
            {3, {5.0/12, 8.0/12, -1.0/12}},
            {4, {9.0/24, 19.0/24, -5.0/24, -1.0/24}},
            {5, {251.0/720, 646.0/720, -264.0/720, -106.0/720, 19.0/720}}
        };
        if(am_coeffs.find(order) == am_coeffs.end()){
            throw std::invalid_argument("Invalid order for Adams-Moulton");
        }
        coefficients_ = am_coeffs.at(order);
    }

    void integrateAM1(const DynamicFunction& Func, double t_end, double h) {
        //
        int n = static_cast<int>(t_end / h);
        u_.resize(n + 1);
        u_[0] = Func(0, 0);
        for(int i = 1; i <= n; ++i){
            u_[i] = u_[i-1] + h*(coefficients_[0]*Func(u_[i], i*h));
        }
    }
    void integrateAM2(const DynamicFunction& Func, double t_end, double h) {
        //
        int n = static_cast<int>(t_end / h);
        u_.resize(n + 1);
        u_[0] = Func(0, 0);
        for(int i = 1; i <= n; ++i){
            u_[i] = u_[i-1] + h*(coefficients_[0]*Func(u_[i], i*h) 
                            + coefficients_[1]*Func(u_[i-1], (i-1)*h));
        }
    }
    void integrateAM3(const DynamicFunction& Func, double t_end, double h) {
        //
        int n = static_cast<int>(t_end / h);
        u_.resize(n + 1);
        u_[0] = Func(0, 0);
        u_[1] = RK4(Func, 0, u_[0], h);
        for(int i = 2; i <= n; ++i){
            u_[i] = u_[i-1] + h*(coefficients_[0]*Func(u_[i], i*h) 
                            + coefficients_[1]*Func(u_[i-1], (i-1)*h)
                            + coefficients_[2]*Func(u_[i-2], (i-2)*h));
        }
    }
    void integrateAM4(const DynamicFunction& Func, double t_end, double h) {
        //
        int n = static_cast<int>(t_end / h);
        u_.resize(n + 1);
        u_[0] = Func(0, 0);
        u_[1] = RK4(Func, 0, u_[0], h);
        u_[2] = RK4(Func, h, u_[1], h);
        for(int i = 3; i <= n; ++i){
            u_[i] = u_[i-1] + h*(coefficients_[0]*Func(u_[i], i*h)  
                            + coefficients_[1]*Func(u_[i-1], (i-1)*h)   
                            + coefficients_[2]*Func(u_[i-2], (i-2)*h)
                            + coefficients_[3]*Func(u_[i-3], (i-3)*h));
        }
    }
    void integrateAM5(const DynamicFunction& Func, double t_end, double h) {
        //
        int n = static_cast<int>(t_end / h);
        u_.resize(n + 1);
        u_[0] = Func(0, 0);
        u_[1] = RK4(Func, 0, u_[0], h);
        u_[2] = RK4(Func, h, u_[1], h);
        u_[3] = RK4(Func, h+h, u_[2], h);
        for(int i = 4; i <= n; ++i){
            u_[i] = u_[i-1] + h*(coefficients_[0]*Func(u_[i], i*h) 
                            + coefficients_[1]*Func(u_[i-1], (i-1)*h)
                            + coefficients_[2]*Func(u_[i-2], (i-2)*h) 
                            + coefficients_[3]*Func(u_[i-3], (i-3)*h)   
                            + coefficients_[4]*Func(u_[i-4], (i-4)*h));      
        }
    }
};


//==================================================================

class BDFSolver : public EquationSolver {
private:
    int order_;
    //储存解向量
    std::vector<double> u_;
    //存储系数
    std::vector<double> coefficients_;
public:
    BDFSolver(int order) : order_(order) {}
    void solve(const DynamicFunction& Func, double t_end, double h)override{
        //
    }
    std::string name() const override{
        return "BDF";
    }
private:
    void loadCoefficients(int order) {
        static const std::map<int, std::vector<double>> bdf_coeffs = {
            {1, {1.0, -1.0, 1.0}},
            {2, {1.0, -4.0/3, 1.0/3, 2.0/3}},
            {3, {1.0, -18.0/11, 9.0/11, -2.0/11, 6.0/11}},
            {4, {1.0, -48.0/25, 36.0/25, -16.0/25, 3.0/25, 12.0/25}}
        };
    }

    void integrateBDF1(const DynamicFunction& Func, double t_end, double h) {
        //
        int n = static_cast<int>(t_end / h);
        u_.resize(n + 1);
        u_[0] = Func(0, 0);
        for(int i = 1; i <= n; ++i){
            u_[i] = -coefficients_[1]*u_[i-1] + h*coefficients_[2]*Func(u_[i], i*h);
        }
    }
    void integrateBDF2(const DynamicFunction& Func, double t_end, double h) {
        //
        int n = static_cast<int>(t_end / h);
        u_.resize(n + 1);
        u_[0] = Func(0, 0);
        u_[1] = RK4(Func, 0, u_[0], h);
        for(int i = 2; i <= n; ++i){
            u_[i] = -coefficients_[1]*u_[i-1] 
                    -coefficients_[2]*u_[i-2] + h*coefficients_[3]*Func(u_[i], i*h);
        }
    }
    void integrateBDF3(const DynamicFunction& Func, double t_end, double h) {
        //
        int n = static_cast<int>(t_end / h);
        u_.resize(n + 1);
        u_[0] = Func(0, 0);
        u_[1] = RK4(Func, 0, u_[0], h);
        u_[2] = RK4(Func, h, u_[1], h);
        for(int i = 3; i <= n; ++i){
            u_[i] = -coefficients_[1]*u_[i-1] 
                    -coefficients_[2]*u_[i-2] 
                    -coefficients_[3]*u_[i-3] + h*coefficients_[4]*Func(u_[i], i*h);
        }
    }
    void integrateBDF4(const DynamicFunction& Func, double t_end, double h) {
        //
        int n = static_cast<int>(t_end / h);
        u_.resize(n + 1);
        u_[0] = Func(0, 0);
        u_[1] = RK4(Func, 0, u_[0], h);
        u_[2] = RK4(Func, h, u_[1], h);
        u_[3] = RK4(Func, h+h, u_[2], h);
        for(int i = 4; i <= n; ++i){
            u_[i] = -coefficients_[1]*u_[i-1] 
                    -coefficients_[2]*u_[i-2] 
                    -coefficients_[3]*u_[i-3] 
                    -coefficients_[4]*u_[i-4] + h*coefficients_[5]*Func(u_[i], i*h);
        }
    }
};


inline std::unique_ptr<EquationSolver> createAdamsBashforthSolver(int order){
    return std::unique_ptr<EquationSolver>(new AdamsBashforthSolver(order));
}

inline std::unique_ptr<EquationSolver> createAdamsMoultonSolver(int order){
    return std::unique_ptr<EquationSolver>(new AdamsMoultonSolver(order));
}

inline std::unique_ptr<EquationSolver> createBDFSolver(int order){
    return std::unique_ptr<EquationSolver>(new BDFSolver(order));
}











