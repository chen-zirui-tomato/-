#pragma once
#include"Function_with_7_param.h"
#include <iostream>
#include <cmath>
#include <exprtk.hpp>
#include <memory>
#include "vector_plus_and_muliple.h"


class EquationSolver {
public:
    virtual ~EquationSolver() = default;
    virtual void solve(const HardcodedFunction_7& Func, double t_end, int N, 
                       int condition_type, std::string filename, bool is_test = false) = 0;

    virtual std::string name() const = 0;
    virtual std::vector<double> get_u_start() const = 0;
    virtual std::vector<double> get_u_end() const = 0;
    // virtual void print_result(const std::vector<double>& u) const{
    //     std::cout << "Result from " << name() << " solver:\n";
    //     for(size_t i = 0; i < u.size(); ++i){
    //         std::cout << "u(" << i << ") = " << u[i] << std::endl;
    //     }
    // };
    //画图
    void plot(const std::vector<std::vector<double>>& u, const std::string& filename){
        std::ofstream ofs(filename);
        if (!ofs) {
            throw std::runtime_error("Failed to open file: " + filename);
        }
        for(const auto& vec : u){
            ofs << vec[0] << " " << vec[1] << " " << vec[2] << "\n";
        }
    }
    void plot_test(const std::vector<std::vector<double>>& u, const std::string& filename, double t_end, int N){
        std::ofstream ofs(filename);
        if (!ofs) {
            throw std::runtime_error("Failed to open file: " + filename);
        }
        int count = 0;
        double step = t_end / N;
        for(const auto& vec : u){
            ofs << vec[0] << " " << (count++)*step << "\n";
        }
    }
    //误差计算
    double error(const std::vector<double>& u_exact, const std::vector<double>& u_approx){
        double sum = 0.0;
        for(int i = 0; i < u_exact.size(); ++i){
            double diff = u_exact[i] - u_approx[i];
            sum += diff * diff;
        }
        return std::sqrt(sum);
    }
    double error_test(const std::vector<double>& u_exact, double t_end){
        double sum = 0.0;
        return std::sqrt(std::abs(u_exact[0] - exp(-t_end)));
    }
protected:
    //四阶RK方法，用于求解初始解
    std::vector<double> RK4(const HardcodedFunction_7& Func, double t, std::vector<double>& u, double h){
        std::vector<double> k1 = Func(u, t);
        
        std::vector<double> temp1 = k1 * h * 0.5;
        std::vector<double> u2 = u + temp1;
        std::vector<double> k2 = Func(u2, t + h/2);
        
        std::vector<double> temp2 = k2 * h * 0.5;
        std::vector<double> u3 = u + temp2;
        std::vector<double> k3 = Func(u3, t + h/2);
        
        std::vector<double> temp3 = k3 * h;
        std::vector<double> u4 = u + temp3;
        std::vector<double> k4 = Func(u4, t + h);
        
        std::vector<double> sum = k1 + (k2 * 2.0) + (k3 * 2.0) + k4;
        return sum * (h / 6.0) + u;
    }
    //不动点迭代法，用于求解隐格式
    std::vector<double> fixedPointIteration(const std::function<std::vector<double>(const std::vector<double>&)>& G, double h,
                               std::vector<double>& initial_guess, double tolerance = 1e-6, int max_iter = 100){
        std::vector<double> u_k = initial_guess;
        for(int i = 0; i < max_iter; ++i){
            std::vector<double> u_next = G(u_k);
            // if(std::abs(u_next - u_k) < tolerance){
            //     return u_next;
            // }
            double sum = 0.0;
            for (int i = 0; i < u_k.size(); ++i) {
                double diff = u_next[i] - u_k[i];
                sum += diff * diff;
            }
            // std::cerr << "sum: " << std::sqrt(sum) << std::endl;
            
            if(std::sqrt(sum) < tolerance)
                return u_next;
            u_k = u_next;
        }
        std::cerr << "fixedPointIteration: Failed to converge " << std::endl;
        return {};
    }
    //Newton迭代法，用于求解隐格式
    std::vector<double> newtonIteration(
        const std::function<std::vector<double>(const std::vector<double>&)>& G,
        const std::vector<double>& u_init,
        double tol = 1e-6, int max_iter = 20) 
    {
        std::vector<double> u = u_init;
        int n = u.size();
        for (int iter = 0; iter < max_iter; ++iter) {
            std::vector<double> F = G(u);  // G(u) = u - RHS(u)
            
            // 构造Jacobian（数值微分）
            std::vector<std::vector<double>> J(n, std::vector<double>(n));
            const double epsilon = 1e-8;
            for (int j = 0; j < n; ++j) {
                std::vector<double> u_perturb = u;
                u_perturb[j] += epsilon;
                std::vector<double> F_perturb = G(u_perturb);
                for (int i = 0; i < n; ++i) {
                    J[i][j] = (F_perturb[i] - F[i]) / epsilon;
                }
            }
    
            // 解 J * delta = -F
            std::vector<double> delta = gaussianElimination(J, F, true); // true 表示右边为 -F
    
            // 更新解
            for (int i = 0; i < n; ++i) {
                u[i] += delta[i];
            }
    
            // 检查收敛
            double err = 0.0;
            for (double x : delta) err += x * x;
            if (std::sqrt(err) < tol) return u;
        }
        throw std::runtime_error("Newton iteration failed to converge.");
    }
    std::vector<double> gaussianElimination(std::vector<std::vector<double>> A, std::vector<double> b, bool negate_b = false) {
        int n = A.size();
        if (negate_b) {
            for (int i = 0; i < n; ++i) b[i] = -b[i];
        }
        // Forward elimination
        for (int i = 0; i < n; ++i) {
            // Find pivot
            int max_row = i;
            for (int k = i + 1; k < n; ++k) {
                if (std::abs(A[k][i]) > std::abs(A[max_row][i]))
                    max_row = k;
            }
            std::swap(A[i], A[max_row]);
            std::swap(b[i], b[max_row]);
    
            // Eliminate
            for (int k = i + 1; k < n; ++k) {
                double factor = A[k][i] / A[i][i];
                for (int j = i; j < n; ++j) {
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }
    
        // Back substitution
        std::vector<double> x(n);
        for (int i = n - 1; i >= 0; --i) {
            x[i] = b[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= A[i][j] * x[j];
            }
            x[i] /= A[i][i];
        }
        return x;
    }
    int detect_period(const std::vector<std::vector<double>>& u, double tol = 1e-4) {
        std::vector<double> ref = u[0];
        double sum = 0;
        double sum_ = 1000; 
        int count = 1;
        for (int i = 100; i < u.size(); ++i) {
            for(int j = 0; j < u[0].size(); ++j){
                sum += (u[i][j] - ref[j]) * (u[i][j] - ref[j]);
                if (std::sqrt(sum) < tol) {
                    return i;
                }
                if(sum_ > std::sqrt(sum)){
                    sum_ = std::sqrt(sum);
                    count = i;
                }
            }
        }
        std::cerr << "周期未检测到（未返回初始值附近）" << std::endl;
        return count;
    }
};

class AdamsBashforthSolver : public EquationSolver {
    private:
        int order_;
        std::vector<std::vector<double>> u_;
        std::vector<double> coefficients_;
        std::vector<double> u_0;
        std::vector<double> u_1 = {0.994, 0, 0, 0, -2.0015851063790825224, 0}; 
        std::vector<double> u_2 = {0.879779227778, 0, 0, 0, -0.379677780949, 0};
        std::vector<double> u_3 = {1.0,0,0,0,0,0};
    
    public:
        AdamsBashforthSolver(int order) : order_(order) {
            loadCoefficients(order);
        }
    
        void solve(const HardcodedFunction_7& Func, double t_end, int N,
                   int condition_type, std::string filename, bool is_test = false) override {
            switch(condition_type){
                case 1: u_0 = u_1; break;
                case 2: u_0 = u_2; break;
                case 3: u_0 = u_3; break;
                default: throw std::invalid_argument("Invalid condition type");
            }
    
            integrateGeneralAB(Func, t_end, N);
    
            double err;
            if (is_test){
                plot_test(u_, filename, t_end, N);
                int N = u_.size() - 1;
                err = error_test(u_[N], t_end);
            }
            else{
                double h = t_end / N;
                plot(u_, filename);
                err = error(u_[0], u_[u_.size()-1]);
                // std::cerr << "周期 = " << detect_period(u_) << std::endl;
            }
            std::cout << "[AB" << order_ << "] 误差 = " << err <<", " << std::endl;
        }
    
        std::string name() const override {
            return "Adams-Bashforth";
        }

        std::vector<double> get_u_start() const override{
            return u_[0];
        }
        std::vector<double> get_u_end() const override{
            return u_[u_.size()-1];
        }
    
    private:
        void loadCoefficients(int order) {
            static const std::map<int, std::vector<double>> ab_coeffs = {
                {1, {1.0}},
                {2, {3.0/2, -1.0/2}},
                {3, {23.0/12, -16.0/12, 5.0/12}},
                {4, {55.0/24, -59.0/24, 37.0/24, -9.0/24}}
            };
            auto it = ab_coeffs.find(order);
            if (it == ab_coeffs.end()) throw std::invalid_argument("Unsupported AB order");
            coefficients_ = it->second;
        }
    
        void integrateGeneralAB(const HardcodedFunction_7& Func, double t_end, int N) {
            int n = N;
            double h = t_end / N;
            u_.resize(n + 1);
            u_[0] = u_0;
    
            // 使用 RK4 启动前 order_-1 步
            for (int i = 1; i < order_; ++i) {
                u_[i] = RK4(Func, (i-1)*h, u_[i-1], h);
            }
    
            // Adams-Bashforth 主循环
            for (int i = order_; i <= n; ++i) {
                std::vector<double> sum(6, 0.0);
                for (int j = 0; j < order_; ++j) {
                    std::vector<double> f = Func(u_[i - 1 - j], (i - 1 - j)*h);
                    for (int k = 0; k < 6; ++k)
                        sum[k] += coefficients_[j] * f[k];
                }
                u_[i] = u_[i-1];
                for (int k = 0; k < 6; ++k)
                    u_[i][k] += h * sum[k];
            }
        }
    };
    
//=================================================================

class AdamsMoultonSolver : public EquationSolver {
    private:
        int order_;
        std::vector<std::vector<double>> u_;
        std::vector<double> coefficients_;
        std::vector<double> u_0;
        std::vector<double> u_1 = {0.994, 0, 0, 0, -2.0015851063790825224, 0}; 
        std::vector<double> u_2 = {0.879779227778, 0, 0, 0, -0.379677780949, 0};
        std::vector<double> u_3 = {1.0, 0, 0, 0, 0, 0};
    
    public:
        AdamsMoultonSolver(int order) : order_(order) {
            loadCoefficients(order);
        }
    
        void solve(const HardcodedFunction_7& Func, double t_end, int N,
                   int condition_type, std::string filename, bool is_test = false) override {
            switch(condition_type){
                case 1: u_0 = u_1; break;
                case 2: u_0 = u_2; break;
                case 3: u_0 = u_3; break;
                default: throw std::invalid_argument("Invalid condition type");
            }
    
            integrateGeneralAM(Func, t_end, N);
    
            double err;
            if (is_test){
                plot_test(u_, filename, t_end, N);
                err = error_test(u_.back(), t_end);
            }
            else{
                double h = t_end / N;
                plot(u_, filename);
                err = error(u_[0], u_[u_.size()-1]);
            }
            std::cout << "[AM" << order_ << "] 误差 = " << err <<", " << std::endl;
        }
    
        std::string name() const override {
            return "Adams-Moulton";
        }

        std::vector<double> get_u_start() const override{
            return u_[0];
        }
        std::vector<double> get_u_end() const override{
            return u_[u_.size()-1];
        }
    
    private:
        void loadCoefficients(int order) {
            static const std::map<int, std::vector<double>> am_coeffs = {
                {1, {1.0}},
                {2, {0.5, 0.5}},
                {3, {5.0/12.0, 2.0/3.0, -1.0/12.0}},
                {4, {3.0/8.0, 19.0/24.0, -5.0/24.0, 1.0/24.0}},
                {5, {251.0/720.0, 646.0/720.0, -264.0/720.0, 106.0/720.0, -19.0/720.0}}
            };
            auto it = am_coeffs.find(order);
            if (it == am_coeffs.end()) throw std::invalid_argument("Unsupported AM order");
            coefficients_ = it->second;
        }
    
        void integrateGeneralAM(const HardcodedFunction_7& Func, double t_end, int N) {
            double h = t_end / N;
            int n = N;
            u_.resize(n + 1);
            u_[0] = u_0;
    
            // 前 order_-1 步用 RK4 初始化
            for (int i = 1; i < order_; ++i) {
                u_[i] = RK4(Func, (i - 1) * h, u_[i - 1], h);
            }
    
            for (int i = order_; i <= n; ++i) {
                double t_next = i * h;
                std::vector<double> u_prev = u_[i - 1];
    
                std::vector<double> history(u_0.size(), 0.0);
                for (int j = 0; j < order_-1; ++j) {
                    std::vector<double> fj = Func(u_[i - 1 - j], (i - 1 - j) * h);
                    for (size_t k = 0; k < fj.size(); ++k) {
                        history[k] += coefficients_[j+1] * fj[k];
                    }
                }
    
                std::vector<double> u_predict = RK4(Func, (i - 1) * h, u_prev, h);
                
                // auto G = [&](std::vector<double> u) -> std::vector<double>{
                //     return u_prev + h*(coefficients_[0]*Func(u, t_next) + history);
                // };
                // u_[i] = fixedPointIteration(G, h, u_predict, 1e-6, 100);
                auto G = [&](std::vector<double> u_trial) -> std::vector<double> {
                    std::vector<double> Fu = Func(u_trial, t_next);
                    std::vector<double> result(Fu.size());
                    for (size_t k = 0; k < Fu.size(); ++k) {
                        double rhs = u_prev[k] + h * (coefficients_[0] * Fu[k] + history[k]);
                        result[k] = u_trial[k] - rhs;
                    }
                    return result;
                };
                u_[i] = newtonIteration(G, u_predict);
                // u_[i] = u_prev + h*(coefficients_[0]*Func(u_prev, t_next) + history);
            }
        }
    };

//==================================================================
class BDFSolver : public EquationSolver {
    private:
        int order_;
        std::vector<std::vector<double>> u_;
        std::vector<double> coefficients_; // {a0, a1, a2, ..., as}
        std::vector<double> u_0;
        std::vector<double> u_1 = {0.994, 0, 0, 0, -2.0015851063790825224, 0}; 
        std::vector<double> u_2 = {0.879779227778, 0, 0, 0, -0.379677780949, 0};
        std::vector<double> u_3 = {1.0, 0, 0, 0, 0, 0};
    
    public:
        BDFSolver(int order) : order_(order) {
            loadCoefficients(order);
        }
    
        void solve(const HardcodedFunction_7& Func, double t_end, int N,
                   int condition_type, std::string filename, bool is_test = false) override {
            switch (condition_type) {
                case 1: u_0 = u_1; break;
                case 2: u_0 = u_2; break;
                case 3: u_0 = u_3; break;
                default: throw std::invalid_argument("Invalid condition type");
            }
    
            integrateGeneralBDF(Func, t_end, N);
    
            double h = t_end / N;
            if (is_test) {
                plot_test(u_, filename, t_end, N);
                std::cout << "[BDF" << order_ << "] test误差 = " << error_test(u_.back(), t_end) << std::endl;
            } else {
                plot(u_, filename);
                std::cout << "[BDF" << order_ << "] 误差 = " << error(u_[0], u_.back()) <<", " << std::endl;
    
                // int idx = detect_period(u_, 1e-4);
                // if (idx > 0) {
                //     std::cout << "[BDF" << order_ << "] 检测周期 ≈ " << h * idx << std::endl;
                //     std::cout << "[BDF" << order_ << "] 周期误差 = " << error(u_[0], u_[idx]) << std::endl;
                // }
            }
        }
    
        std::string name() const override {
            return "BDF";
        }

        std::vector<double> get_u_start() const override{
            return u_[0];
        }
        std::vector<double> get_u_end() const override{
            return u_[u_.size()-1];
        }
    
    private:
        void loadCoefficients(int order) {
            static const std::map<int, std::vector<double>> bdf_coeffs = {
                {1, {1.0, 1.0}},                               // u_{n+1} - u_n = h*f_{n+1}
                {2, {2.0/3, 4.0/3, -1.0/3}},                   // u_{n+1} - 4/3 u_n + 1/3 u_{n-1} = h*2/3 f_{n+1}
                {3, {6.0/11, 18.0/11, -9.0/11, 2.0/11}},       // ...
                {4, {12.0/25, 48.0/25, -36.0/25, 16.0/25, -3.0/25}}
            };
            auto it = bdf_coeffs.find(order);
            if (it == bdf_coeffs.end()) throw std::invalid_argument("Unsupported BDF order");
            coefficients_ = it->second;
        }
    
        void integrateGeneralBDF(const HardcodedFunction_7& Func, double t_end, int N) {
            double h = t_end / N;
            int n = N;
            u_.resize(n + 1);
            u_[0] = u_0;
    
            // 启动解：RK4 初始化前 order-1 步
            for (int i = 1; i < order_; ++i) {
                u_[i] = RK4(Func, (i - 1) * h, u_[i - 1], h);
            }
    
            for (int i = order_; i <= n; ++i) {
                double t_next = i * h;
                std::vector<double> history(u_0.size(), 0.0);
    
                // 累加历史项系数（u_n, u_{n-1}, ...）
                for (int j = 1; j <= order_; ++j) {
                    const auto& uj = u_[i - j];
                    for (size_t k = 0; k < uj.size(); ++k) {
                        history[k] += coefficients_[j] * uj[k];
                    }
                }
    
                std::vector<double> u_prev = u_[i - 1];
                std::vector<double> u_predict = RK4(Func, (i - 1) * h, u_prev, h);
    
                // 隐式 G(u) = u - (history + h*a0*f(u))
                auto G = [&](const std::vector<double>& u_trial) -> std::vector<double> {
                    std::vector<double> f = Func(u_trial, t_next);
                    std::vector<double> result(u_0.size());
                    for (size_t k = 0; k < result.size(); ++k) {
                        result[k] = u_trial[k] - (history[k] + h * coefficients_[0] * f[k]);
                    }
                    return result;
                };
    
                u_[i] = newtonIteration(G, u_predict);  // 推荐 Newton
            }
        }
    };

//==================================================================
class CRKSolver : public EquationSolver {
private:
    int order_;
    std::vector<std::vector<double>> u_;
    std::vector<double> coefficients_;
    std::vector<double> u_0;
    std::vector<double> u_1 = {0.994, 0, 0, 0, -2.0015851063790825224, 0}; 
    std::vector<double> u_2 = {0.879779227778, 0, 0, 0, -0.379677780949, 0};

public:
    CRKSolver(int order) : order_(order) {}

    void solve(const HardcodedFunction_7& Func, double t_end, int N,
               int condition_type, std::string filename, bool is_test = false) override {
        switch(condition_type){
            case 1: u_0 = u_1; break;
            case 2: u_0 = u_2; break;
            default: throw std::invalid_argument("Invalid condition type");
        }
        int n = N;
        double h = t_end / n;
        u_.resize(n + 1);
        u_[0] = u_0;

        for(int i = 1; i <= n; ++i){
            u_[i] = RK4(Func, (i-1)*h, u_[i-1], h);
        }
        
        plot(u_, filename);
        double err = error(u_[0], u_[u_.size()-1]);
        std::cout << "[CRK] 误差 = " << err <<", " << std::endl;
    }
            
    std::string name() const override {
        return "classical Runge-Kutta method";
    }
    std::vector<double> get_u_start() const override{
        return u_[0];
    }
    std::vector<double> get_u_end() const override{
        return u_[u_.size()-1];
    }

public:
};

//==================================================================
class ESDIRKSolver : public EquationSolver {
private:
    int order_;
    std::vector<std::vector<double>> u_;
    std::vector<double> coefficients_;
    std::vector<double> u_0;
    std::vector<double> u_1 = {0.994, 0, 0, 0, -2.0015851063790825224, 0}; 
    std::vector<double> u_2 = {0.879779227778, 0, 0, 0, -0.379677780949, 0};

public:
    ESDIRKSolver(int order) : order_(order) {}

    void solve(const HardcodedFunction_7& Func, double t_end, int N,
               int condition_type, std::string filename, bool is_test = false) override {
        switch(condition_type){
            case 1: u_0 = u_1; break;
            case 2: u_0 = u_2; break;
            default: throw std::invalid_argument("Invalid condition type");
        }
        int n = N;
        double h = t_end / N;
        u_.resize(n + 1);
        u_[0] = u_0;

        for(int i = 1; i <= n; ++i){
            u_[i] = ESDIRK4(Func, (i-1)*h, u_[i-1], h);
        }

        plot(u_, filename);
        double err = error(u_[0], u_[u_.size()-1]);
        std::cout << "[ESDIRK] 误差 = " << err <<", " << std::endl;
    }
            
    std::string name() const override {
        return "ESDIRK method";
    }
    std::vector<double> get_u_start() const override{
        return u_[0];
    }
    std::vector<double> get_u_end() const override{
        return u_[u_.size()-1];
    }

private:

    std::vector<double> ESDIRK4(const HardcodedFunction_7& Func, double t, const std::vector<double>& u, double h){
        std::vector<double> y1 = Func(u, t);

        // std::vector<double> y2 = Func(u + (y1 + y2)*h*0.25, t + h/2);
        auto G2 = [&](std::vector<double> u_trial) -> std::vector<double> {
            std::vector<double> result = Func(u + (y1 + u_trial)*h*0.25, t + h/2);
            return u_trial - result;
        };
        std::vector<double> y2_predict = RK4(Func, t, y1, h/2);
        std::vector<double> y2 = newtonIteration(G2, y2_predict);

        // std::vector<double> y3 = Func(u + h*(8611.0/62500*y1 - 1743.0/31250*y2 + 1.0/4*y3), t + 83.0/250*h);
        auto G3 = [&](std::vector<double> u_trail) -> std::vector<double>{
            std::vector<double> result = Func(u + h*(8611.0/62500*y1 - 1743.0/31250*y2 + 1.0/4*u_trail), t + 83.0/250*h);
            return u_trail - result;
        };
        std::vector<double> y3_predict = RK4(Func, t, y2, 83.0/250*h);
        std::vector<double> y3 = newtonIteration(G3, y3_predict);

        // std::vector<double> y4 = Func(u + h*(5012029.0/34652500*y1 - 654441.0/2922500*y2
        //                                    + 174375.0/388108*3*y3 + 1.0/4*y4), t + 31.0/50*h);
        auto G4 = [&](std::vector<double> u_trail) -> std::vector<double> {
            std::vector<double> result =  Func(u + h*(5012029.0/34652500*y1 - 654441.0/2922500*y2
                                                + 174375.0/388108*y3 + 1.0/4*u_trail), t + 31.0/50*h);
            return u_trail - result;
        };
        std::vector<double> y4_predict = RK4(Func, t, y3, 31.0/50*h);
        std::vector<double> y4 = newtonIteration(G4, y4_predict);

        // std::vector<double> y5 = Func(u + h*(15267082809.0/155376265600*y1 - 71443401.0/120774400*y2 
        //                                 + 730878875.0/902184768*y3 + 2285395.0/8070912*y4 + 1.0/4*y5), t + 17.0/20*h);
        auto G5 = [&](std::vector<double> u_trail) -> std::vector<double> {
            std::vector<double> result =  Func(u + h*(15267082809.0/155376265600*y1 - 71443401.0/120774400*y2 
                                        + 730878875.0/902184768*y3 + 2285395.0/8070912*y4 + 1.0/4*u_trail), t + 17.0/20*h);
            return u_trail - result;
        };
        std::vector<double> y5_predict = RK4(Func, t, y4, 17.0/20*h);
        std::vector<double> y5 = newtonIteration(G5, y5_predict);

        // std::vector<double> y6 = Func(u + h*(82889.0/524892*y1 + 15625.0/83664*y3 + 69875.0/102672*y4
        //                                 - 2260.0/8211*y5 + 1.0/4*y6), t + h);
        auto G6 = [&](std::vector<double> u_trail) -> std::vector<double> {
            std::vector<double> result = Func(u + h*(82889.0/524892*y1 + 15625.0/83664*y3 + 69875.0/102672*y4
                                        - 2260.0/8211*y5 + 1.0/4*u_trail), t + h);
            return u_trail - result;
        };
        std::vector<double> y6_predict = RK4(Func, t, y5, h);
        std::vector<double> y6 = newtonIteration(G6, y6_predict);

        std::vector<double> result = u + h*(82889.0/524892*y1 + 15625.0/83664*y3 + 69875.0/102672*y4
                                            - 2260.0/8211*y5 + 1.0/4*y6); 
        return result;
    }
};

//==================================================================
class GLRKSolver : public EquationSolver {
private:
    int order_;
    std::vector<std::vector<double>> u_;
    std::vector<std::vector<double>> coefficients_;
    std::vector<double> u_0;
    std::vector<double> u_1 = {0.994, 0, 0, 0, -2.0015851063790825224, 0}; 
    std::vector<double> u_2 = {0.879779227778, 0, 0, 0, -0.379677780949, 0};

public:
    GLRKSolver(int order) : order_(order) {}

    void solve(const HardcodedFunction_7& Func, double t_end, int N,
               int condition_type, std::string filename, bool is_test = false) override {
        loadCoefficients(order_);

        switch(condition_type){
            case 1: u_0 = u_1; break;
            case 2: u_0 = u_2; break;
            default: throw std::invalid_argument("Invalid condition type");
        }

        integrateGaussLegendre(Func, t_end, N);

        plot(u_, filename);
        double err = error(u_[0], u_[u_.size()-1]);
        std::cout << "[GLRK"<<order_<<"] 误差 = " << err << std::endl;
    }
            
    std::string name() const override {
        return "Gauss-Legendre Runge-Kutta method";
    }
    std::vector<double> get_u_start() const override{
        return u_[0];
    }
    std::vector<double> get_u_end() const override{
        return u_[u_.size()-1];
    }

private:
    void loadCoefficients(int order) {
        std::vector<std::vector<double>> coeffs_2, coeffs_3, coeffs_4, coeffs_5;
        switch (order) {
            case 2:
                coeffs_2 = {{(3-sqrt(3))/6, 1.0/4, (3-2*sqrt(3))/12},
                            {(3+sqrt(3))/6, (3+2*sqrt(3))/12, 1.0/4},
                            {0xffffff, 1.0/2, 1.0/2}}; break;
            case 3:
                coeffs_3 = {{(5-sqrt(15))/10, 5.0/36, (2.0/9-sqrt(15)/15), 5.0/36-sqrt(15)/30},
                            {1.0/2, 5.0/36+sqrt(15)/24, 2.0/9, 5.0/36-sqrt(15)/24},
                            {(5+sqrt(15))/10, 5.0/36-sqrt(15)/30, 2.0/9+sqrt(15)/15, 5.0/36},
                            {0xffffff, 15.0/18, 4.0/9, 5.0/18}}; break;
            case 4:
                coeffs_4 = {{0.06943184420297371, 0.08696371128436343, -0.026604180084998794, 0.012627462689404725, -0.003555149685795685},
                            {0.33000947820757187, 0.1881181174998681, 0.16303628871563644, -0.027880428602470822, 0.00673550059453814},
                            {0.6699905217924281, 0.16719192197418906, 0.35395300603374347, 0.16303628871563677, -0.014190694931141168},
                            {0.9305681557970262, 0.17748257225452302, 0.3134451147418678, 0.3526767575162716, 0.08696371128436374},
                            {0xffffff, 0.17392742256872684, 0.3260725774312731, 0.3260725774312731, 0.17392742256872684}
                        }; break;
            case 5:
                coeffs_5 = {{0.04691007703066802, 0.0592317212640473, -0.019570364359076053, 0.011254400818642964, -0.005593793660812192, 0.0015881129678660005},
                            {0.23076534494715845, 0.12815100567004536, 0.11965716762484149, -0.02459211461964212, 0.010318280670683319, -0.0027689943987695927},
                            {0.5, 0.11377628800422483, 0.26000465168064113, 0.14222222222222242, -0.02069031643095834, 0.004687154523869958},
                            {0.7692346550528415, 0.12123243692686453, 0.2289960545789994, 0.3090365590640866, 0.11965716762484174, -0.009687563141950701},
                            {0.9530899229693319, 0.11687532956022907, 0.24490812891049493, 0.2731900436258015, 0.25888469960875893, 0.05923172126404763},
                            {0xffffff, 0.11846344252809471, 0.2393143352496831, 0.2844444444444445, 0.2393143352496831, 0.11846344252809471}
                        }; break;
            default: throw std::invalid_argument("Unsupported order");
        }

        static const std::map<int, std::vector<std::vector<double>>> ab_coeffs = {
            {2, coeffs_2},
            {3, coeffs_3},
            {4, coeffs_4},
            {5, coeffs_5}
        };

        auto it = ab_coeffs.find(order);
        if (it == ab_coeffs.end()) throw std::invalid_argument("Unsupported AB order");
        coefficients_ = it->second;
    }

    void integrateGaussLegendre(const HardcodedFunction_7& Func, double t_end, int N) {
        int n = N;
        double h = t_end / n;
        u_.resize(n+1);
        u_[0] = u_0;
        /*Y外层的大小为order_,表示GL方法的步数
          Y内层的大小为6,即存储y1到y|order_的结果，每个yi是六阶向量*/
        std::vector<std::vector<double>> Y;

        for(int i = 1; i <= n; ++i){
            // double t = i*h;
            // //接下来计算求解Y
            // for(int j = 0; j < order_; ++j){
            //     std::vector<double> temp;
            //     for(int k = 1; k < order_; ++k){
            //         temp += coefficients_[j][k]*Y[k];
            //     }
            //     Y[j] = Func(u_[i-1] + temp, t + coefficients_[j][0]*h);
            // }
            // /*现在需要一个高阶的newton法来求解Y，但是newton法需要精确近似解，
            // 可以在第一次迭代用二分法？之后换牛顿法？*/
            std::vector<double> u_prev = u_[i-1];
            double t_prev = (i - 1) * h;
    
            // 初始猜测：每个 Y[j] = u_prev
            std::vector<std::vector<double>> Y(order_, u_prev);
            for (int j = 0; j < order_; ++j) {
                std::vector<double> f0 = Func(u_prev, t_prev + coefficients_[j][0] * h);
                for (int k = 0; k < 6; ++k)
                    Y[j][k] += h * coefficients_[j][1] * f0[k];
            }
    
            // 构造 G 函数用于 Newton
            auto G = [&](const std::vector<std::vector<double>>& Y_input) -> std::vector<std::vector<double>> {
                std::vector<std::vector<double>> result(order_, std::vector<double>(6, 0.0));
                for (int j = 0; j < order_; ++j) {
                    std::vector<double> sum(6, 0.0);
                    for (int k = 0; k < order_; ++k) {
                        auto fYk = Func(Y_input[k], t_prev + coefficients_[k][0] * h);  // c_k h
                        for (int l = 0; l < 6; ++l) {
                            sum[l] += coefficients_[j][k + 1] * fYk[l];  // A[j][k]
                        }
                    }
                    for (int l = 0; l < 6; ++l) {
                        result[j][l] = Y_input[j][l] - u_prev[l] - h * sum[l];
                    }
                }
                return result;
            };

            Y = highOrderNewton(G, Y, 1e-10, 50);

            // 计算 u_i
            std::vector<double> temp(6, 0.0);
            for(int j = 0; j < order_; ++j){
                temp += coefficients_[order_][j+1]*Y[j];
            }
            u_[i] = u_prev + h*temp;
        }
    };


    std::vector<std::vector<double>> highOrderNewton(
        const std::function<std::vector<std::vector<double>>(const std::vector<std::vector<double>>&)>& G,
        const std::vector<std::vector<double>>& Y_init,
        double tol = 1e-6, int max_iter = 20)
    {
        std::vector<std::vector<double>> Y = Y_init;  // shape: s × d
        int s = Y.size();        // number of stages
        int d = Y[0].size();     // dimension of each Y_j
        int n = s * d;           // total system dimension
        const double epsilon = 1e-8;
    
        for (int iter = 0; iter < max_iter; ++iter) {
            // 计算 G(Y)
            std::vector<std::vector<double>> G_Y = G(Y);  // size: s × d
    
            // 将 G_Y 拉平成一维向量 F
            std::vector<double> F(n);
            for (int i = 0; i < s; ++i)
                for (int j = 0; j < d; ++j)
                    F[i * d + j] = G_Y[i][j];
    
            // 数值构造雅可比矩阵 J (n × n)
            std::vector<std::vector<double>> J(n, std::vector<double>(n, 0.0));
            for (int k = 0; k < n; ++k) {
                int i_k = k / d, j_k = k % d;
                std::vector<std::vector<double>> Y_perturb = Y;
                Y_perturb[i_k][j_k] += epsilon;
    
                std::vector<std::vector<double>> G_perturb = G(Y_perturb);
                for (int l = 0; l < s; ++l)
                    for (int m = 0; m < d; ++m) {
                        int idx = l * d + m;
                        J[idx][k] = (G_perturb[l][m] - G_Y[l][m]) / epsilon;
                    }
            }
    
            // 求解线性系统 J * delta = -F
            std::vector<double> delta = gaussianElimination(J, F, true);  // true 表示 -F
    
            // 更新 Y
            for (int i = 0; i < s; ++i)
                for (int j = 0; j < d; ++j)
                    Y[i][j] += delta[i * d + j];
    
            // 检查收敛
            double err = 0.0;
            for (double x : delta) err += x * x;
            if (std::sqrt(err) < tol) return Y;
        }
    
        throw std::runtime_error("High-order Newton iteration failed to converge.");
    }
    
};

//==================================================================
class FSolver : public EquationSolver {
private:
    int order_;
    std::vector<std::vector<double>> u_;
    std::vector<double> coefficients_;
    std::vector<double> u_0;
    std::vector<double> u_1 = {0.994, 0, 0, 0, -2.0015851063790825224, 0}; 
    std::vector<double> u_2 = {0.879779227778, 0, 0, 0, -0.379677780949, 0};
public:
    FSolver(int order) : order_(order) {};

    void solve(const HardcodedFunction_7& Func, double t_end, int N,
               int condition_type, std::string filename, bool is_test = false) override {

        switch(condition_type){
            case 1: u_0 = u_1; break;
            case 2: u_0 = u_2; break;
            default: throw std::invalid_argument("Invalid condition type");
        }

        double E_abs = 1e-7;
        double E_rel = 1e-4;
        //安全因子
        double rho = 0.9;
        //误差阶数
        int q = 4;
        //缩放因子范围
        double rho_max = 1.5;
        double rho_min = 0.2;

        N = 100;
        int n = N;
        //N的作用只有定第一步的步长？
        double h = t_end / n;
        u_.clear();
        std::vector<double> u = u_0;
        u_.push_back(u);

        // for(int i = 1; i <= n; ++i){
        //     U1 = F(Func, (i-1)*h, u_[i-1], h, true);
        //     U2 = F(Func, (i-1)*h, u_[i-1], h, false);
        //     double eps = 0;
        //     h = h*std::min(rho_max, std::max(rho_min, rho*std::pow(eps, -1.0/(q+1))));
        //     for(size_t j = 0; j < u_[i-1].size(); ++j){
        //         double epsi = E_abs + E_rel*std::abs(u_[i-1][j]);
        //         eps += ((U1[j] - U2[j])/epsi)*((U1[j] - U2[j])/epsi);
        //     }
        //     eps = std::sqrt(eps/N);
        //     if(eps > 1.0){
        //         --i;
        //         ++n;
        //     }
        // }
        double t = 0;
        while(t < t_end){
            if(t+h > t_end)
                h = t_end - t;

            std::vector<double> U1 = F(Func, t, u, h, true);
            std::vector<double> U2 = F(Func, t, u, h, false);

            double eps = 0;
            for(size_t j = 0; j < u.size(); ++j){
                double epsi = E_abs + E_rel*std::abs(u[j]);
                double delta = (U1[j] - U2[j])/epsi;
                eps += delta*delta;
            }
            eps = std::sqrt(eps/u.size());

            if(eps <= 1.0){
                t += h;
                u = U2;
                u_.push_back(u);
            }

            double scale = std::pow(rho/eps, 1.0/(q+1));
            h = h*std::min(rho_max, std::max(rho_min, scale));
            if(h < 1e-10)
                throw std::runtime_error("Step size is too small.");
        }

        plot(u_, filename);
        double err = error(u_[0], u_[u_.size()-1]);
        std::cout << "[Fehlberg method] 误差 = " << err <<", " << std::endl;
    }

    std::string name() const override {
        return "Fehlberg method";
    }
    std::vector<double> get_u_start() const override{
        return u_[0];
    }
    std::vector<double> get_u_end() const override{
        return u_[u_.size()-1];
    }
private:
    std::vector<double> F(const HardcodedFunction_7& Func, double t, const std::vector<double>& u, double h, bool is_first) {
        auto k1 = Func(u, t);
        auto k2 = Func(u + h * 0.25 * k1, t + 0.25 * h);
        auto k3 = Func(u + h * (3.0/32 * k1 + 9.0/32 * k2), t + 3.0/8 * h);
        auto k4 = Func(u + h * (1932.0/2197 * k1 - 7200.0/2197 * k2 + 7296.0/2197 * k3), t + 12.0/13 * h);
        auto k5 = Func(u + h * (439.0/216 * k1 - 8.0 * k2 + 3680.0/513 * k3 - 845.0/4104 * k4), t + h);
        auto k6 = Func(u + h * (-8.0/27 * k1 + 2.0 * k2 - 3544.0/2565 * k3 + 1859.0/4104 * k4 - 11.0/40 * k5), t + 0.5 * h);

        if (is_first) {
            return u + h * (25.0/216 * k1 + 1408.0/2565 * k3 + 2197.0/4104 * k4 - 1.0/5 * k5); // 4阶
        } else {
            return u + h * (16.0/135 * k1 + 6656.0/12825 * k3 + 28561.0/56430 * k4 - 9.0/50 * k5 + 2.0/55 * k6); // 5阶
        }
    }
};

//====================================================================
class DPSolver : public EquationSolver {
private:
    int order_;
    std::vector<std::vector<double>> u_;
    std::vector<double> coefficients_;
    std::vector<double> u_0;
    std::vector<double> u_1 = {0.994, 0, 0, 0, -2.0015851063790825224, 0}; 
    std::vector<double> u_2 = {0.879779227778, 0, 0, 0, -0.379677780949, 0};
public:
    DPSolver(int order) : order_(order) {};

    void solve(const HardcodedFunction_7& Func, double t_end, int N,
                int condition_type, std::string filename, bool is_test = false) override {

        switch(condition_type){
            case 1: u_0 = u_1; break;
            case 2: u_0 = u_2; break;
            default: throw std::invalid_argument("Invalid condition type");
        }

        double E_abs = 1e-7;
        double E_rel = 1e-4;
        //安全因子
        double rho = 0.9;
        //误差阶数
        int q = 4;
        //缩放因子范围
        double rho_max = 1.5;
        double rho_min = 0.2;

        N = 100;
        int n = N;
        //N的作用只有定第一步的步长？
        double h = t_end / n;
        u_.clear();
        std::vector<double> u = u_0;
        u_.push_back(u);

        // for(int i = 1; i <= n; ++i){
        //     U1 = F(Func, (i-1)*h, u_[i-1], h, true);
        //     U2 = F(Func, (i-1)*h, u_[i-1], h, false);
        //     double eps = 0;
        //     h = h*std::min(rho_max, std::max(rho_min, rho*std::pow(eps, -1.0/(q+1))));
        //     for(size_t j = 0; j < u_[i-1].size(); ++j){
        //         double epsi = E_abs + E_rel*std::abs(u_[i-1][j]);
        //         eps += ((U1[j] - U2[j])/epsi)*((U1[j] - U2[j])/epsi);
        //     }
        //     eps = std::sqrt(eps/N);
        //     if(eps > 1.0){
        //         --i;
        //         ++n;
        //     }
        // }
        double t = 0;
        while(t < t_end){
            if(t+h > t_end)
                h = t_end - t;

            std::vector<double> U1 = F(Func, t, u, h, true);
            std::vector<double> U2 = F(Func, t, u, h, false);

            double eps = 0;
            for(size_t j = 0; j < u.size(); ++j){
                double epsi = E_abs + E_rel*std::abs(u[j]);
                double delta = (U1[j] - U2[j])/epsi;
                eps += delta*delta;
            }
            eps = std::sqrt(eps/u.size());

            if(eps <= 1.0){
                t += h;
                u = U2;
                u_.push_back(u);
            }

            double scale = std::pow(rho/eps, 1.0/(q+1));
            h = h*std::min(rho_max, std::max(rho_min, scale));
            if(h < 1e-10)
                throw std::runtime_error("Step size is too small.");
        }
        
        plot(u_, filename);
        double err = error(u_[0], u_[u_.size()-1]);
        std::cout << "[Dormand-Prince method] 误差 = " << err <<", " << std::endl;
    }

    std::string name() const override {
        return "Dormand-Prince method";
    }
    std::vector<double> get_u_start() const override{
        return u_[0];
    }
    std::vector<double> get_u_end() const override{
        return u_[u_.size()-1];
    }
private:
    std::vector<double> F(const HardcodedFunction_7& Func, double t, const std::vector<double>& u, double h, bool is_first){
        std::vector<double> y1 = Func(u, t);
        std::vector<double> y2 = Func(u + h/5*y1, t + h/5);
        std::vector<double> y3 = Func(u + h*(3.0/40*y1 + 9.0/40*y2), t + h*3.0/10);
        std::vector<double> y4 = Func(u + h*(44.0/45*y1 - 56.0/15*y2 + 32.0/9*y3), t + h*4.0/5);
        std::vector<double> y5 = Func(u + h*(19372.0/6561*y1 - 25360.0/2187*y2 + 64448.0/6561*y3 - 212.0/729*y4), t + h*8.0/9);
        std::vector<double> y6 = Func(u + h*(9017.0/3168*y1 - 355.0/33*y2 + 46732.0/5247*y3 + 49.0/176*y4 - 5103.0/18656*y5), t + h);
        std::vector<double> y7 = Func(u + h*(35.0/384*y1 + 500.0/1113*y3 + 125.0/192*y4 - 2187.0/6784*y5 + 11.0/84*y6), t + h*1.0);

        if(is_first){
            return u + h*(35.0/384*y1 + 500.0/1113*y3 + 125.0/192*y4 - 2187.0/6784*y5 + 11.0/84*y6);
        }
        else{
            return u + h*(5179.0/57600*y1 + 7571.0/16695*y3 + 393.0/640*y4 - 92097.0/339200*y5 + 187.0/2100*y6 + 1.0/40*y7);
        }
    }
};

//==================================================================
inline std::unique_ptr<EquationSolver> createAdamsBashforthSolver(int order){
    return std::unique_ptr<EquationSolver>(new AdamsBashforthSolver(order));
}

inline std::unique_ptr<EquationSolver> createAdamsMoultonSolver(int order){
    return std::unique_ptr<EquationSolver>(new AdamsMoultonSolver(order));
}

inline std::unique_ptr<EquationSolver> createBDFSolver(int order){
    return std::unique_ptr<EquationSolver>(new BDFSolver(order));
}

inline std::unique_ptr<EquationSolver> createCRKSolver(int order){
    return std::unique_ptr<EquationSolver>(new CRKSolver(order));
}

inline std::unique_ptr<EquationSolver> createESDIRKSolver(int order){
    return std::unique_ptr<EquationSolver>(new ESDIRKSolver(order));
}

inline std::unique_ptr<EquationSolver> createGLRKSolver(int order){
    return std::unique_ptr<EquationSolver>(new GLRKSolver(order));
}

inline std::unique_ptr<EquationSolver> createFSolver(int order){
    return std::unique_ptr<EquationSolver>(new FSolver(order));
}

inline std::unique_ptr<EquationSolver> createDPSolver(int order){
    return std::unique_ptr<EquationSolver>(new DPSolver(order));
}




