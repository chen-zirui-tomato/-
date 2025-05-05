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
            std::cout << "[AB" << order_ << "] 误差 = " << err << std::endl;
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
            std::cout << "[AM" << order_ << "] 误差 = " << err << std::endl;
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
                std::cout << "[BDF" << order_ << "] 误差 = " << error(u_[0], u_.back()) << std::endl;
    
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


inline std::unique_ptr<EquationSolver> createAdamsBashforthSolver(int order){
    return std::unique_ptr<EquationSolver>(new AdamsBashforthSolver(order));
}

inline std::unique_ptr<EquationSolver> createAdamsMoultonSolver(int order){
    return std::unique_ptr<EquationSolver>(new AdamsMoultonSolver(order));
}

inline std::unique_ptr<EquationSolver> createBDFSolver(int order){
    return std::unique_ptr<EquationSolver>(new BDFSolver(order));
}











