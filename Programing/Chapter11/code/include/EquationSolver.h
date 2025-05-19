#pragma once
#include "Function.h"
#include "vector_plus_and_multiple.h"
#include <Eigen/Dense>
#include <memory>
#include <iostream>

class Solver{
public:
    virtual ~Solver() = default;
    virtual void solve(double r) = 0;
    void construct_D(double h){
        m = 1.0/h - 1;
        D = Eigen::MatrixXd::Zero(m, m);
        double coeff = 1.0/h/h;
        
        for(int i=0; i<m; i++){
            D(i,i) = -2.0*coeff;
            if(i>0)
                D(i,i-1) = 1.0*coeff;
            if(i<m-1)
                D(i,i+1) = 1.0*coeff;
        }
    }
    void plot(const std::string& filename){
        std::ofstream ofs(filename);
        if(!ofs)
            throw std::runtime_error("Failed to open file: " + filename);
        int count = 0;
        std::cout << "u_ size = " << u_.size() << std::endl;
        for(const auto& p : u_){
            ofs << (++count)*h_ << " " << p << std::endl;
            std::cout << "count = " << count << " u = " << p << std::endl;
        }
        ofs.close();
    };
    double u_exact(double t, double x, int N_terms = 100) {
        double sum = 0.0;
        for (int k = 1; k <= N_terms; ++k) {
            double kpi = k * M_PI;
            double Ak = (40.0 / (k * k * M_PI * M_PI)) * (
                -sin(kpi * 9.0 / 20.0) + 
                2.0 * sin(kpi * 0.5) - 
                sin(kpi * 11.0 / 20.0)
            );
            sum += Ak * exp(-k * k * M_PI * M_PI * t) * sin(kpi * x);
        }
        return sum;
    }
    void compute_true_solution(double t, int N_terms = 100) {
        u_ = Eigen::VectorXd::Zero(19);
        for (int j = 1; j <= 19; ++j) {
            double x = j * h_;
            u_(j - 1) = u_exact(t, x, N_terms);
        }
    }
    
protected:
    int m;
    double h_ = 1.0/20;
    Eigen::MatrixXd D;
    
    // //数值解向量
    Eigen::VectorXd u_;
    // 初始条件
    phi I_C;
};

// 求解ODE ：du/dt = vDu
// ===========================================

class GLRKSolver : public Solver {
private:
    int order_ = 2;
    int step;
    std::vector<std::vector<double>> coefficients_= 
                            {{(3-sqrt(3))/6, 1.0/4, (3-2*sqrt(3))/12},      
                             {(3+sqrt(3))/6, (3+2*sqrt(3))/12, 1.0/4},
                             {0xffffff, 1.0/2, 1.0/2}};
public:
    GLRKSolver(int step) : step(step) {}

    void solve(double r) override {
        if(r == 0)
            r = 1.0 / (2 * h_);
        double k = r * h_ * h_;
        construct_D(h_);
        u_ = Eigen::VectorXd::Zero(m);
        for(int i = 1; i <= m; i++)
            u_(i-1) = I_C(i * h_); 

        integrateGaussLegendre(1, k);
    }

private:
    void integrateGaussLegendre(double t_end, double h) {
        auto F = [&](const std::vector<double>& u_vec, double t) -> std::vector<double> {
            Eigen::Map<const Eigen::VectorXd> u_map(u_vec.data(), u_vec.size());
            Eigen::VectorXd du = D * u_map;
            return std::vector<double>(du.data(), du.data() + du.size());
        };
        int n = t_end / h;
        // double h = t_end / n;
        std::vector<std::vector<double>> Y;

        for(int i = 1; i <= step; ++i){
            std::vector<double> u_prev;
            for(size_t j = 0; j < u_.size(); ++j)
                u_prev.push_back(u_(j));
            double t_prev = (i - 1) * h;
    
            // 构造 G 函数用于 Newton
            auto G = [&](const std::vector<std::vector<double>>& Y) -> std::vector<std::vector<double>> {
                std::vector<std::vector<double>> result(order_, std::vector<double>(u_prev.size(), 0.0));
                // j代表我们在计算Y[j]的零求解函数
                for (int j = 0; j < order_; ++j) {
                    std::vector<double> sum(u_prev.size(), 0.0);
                    for (int k = 0; k < order_; ++k) {
                        sum += coefficients_[j][k+1]*Y[k];
                    }
                    result[j] = Y[j] - F(u_prev + h*sum, t_prev + coefficients_[j][0] * h);
                }
                return result;
            };  

            std::vector<std::vector<double>> Y_predict(order_, u_prev);
            Y = highOrderNewton(G, Y_predict, 1e-10, 20);

            // 计算 u_i
            std::vector<double> temp(u_prev.size(), 0.0);
            for(int j = 0; j < order_; ++j){
                temp += coefficients_[order_][j+1]*Y[j];
            }
            u_prev = u_prev + h*temp;
            for(size_t j = 0; j < u_.size(); ++j)
                u_(j) = u_prev[j];
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
    
};


class Backward_Euler : public Solver{
private:
    int step;
public:
    Backward_Euler(int step): step(step){};
    void solve(double r) override{
        if(r == 0)
            r = 1.0/(2*h_);
        double k = r*h_*h_;
        construct_D(h_);
        // 初始化u_
        u_ = Eigen::VectorXd::Zero(m);
        for(int i=1; i<=m ; i++)
            u_(i-1) = I_C(i * h_);
        
        // 迭代求解
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, m);
        Eigen::MatrixXd A = I - k*D;
        for(int i=1; i<=step; i++){
            u_ = btc_step(A, u_);
        }
    };

    Eigen::VectorXd btc_step(const Eigen::MatrixXd& A, const Eigen::VectorXd& u_n){
        return A.fullPivLu().solve(u_n);
    }
};

//theta method that theta = 0.5
class Crank_Nicolson : public Solver {
private:
    int step;
public:
    Crank_Nicolson(int step): step(step){};
    void solve(double r) override {
        if(r == 0)
            r = 1.0 / (2 * h_);
        double k = r * h_ * h_;
        construct_D(h_);

        u_ = Eigen::VectorXd::Zero(m);
        for(int i = 1; i <= m; i++)
            u_(i-1) = I_C(i * h_); 

        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, m);
        Eigen::MatrixXd A = I - 0.5 * k * D;
        Eigen::MatrixXd B = I + 0.5 * k * D;

        for(int i = 1; i <= step; i++) {
            u_ = A.fullPivLu().solve(B * u_);
        }
    }
};

class Forward_Euler : public Solver{
private:
    int step;
public:
    Forward_Euler(int step): step(step){};
    void solve(double r) override{
        if(r == 0)
            r = 1.0/(2*h_);
        double k = r*h_*h_;
        construct_D(h_);
        // 初始化u_
        u_ = Eigen::VectorXd::Zero(m);
        for(int i=1; i<=m ; i++)
            u_(i-1) = I_C(i*h_);
        
        // 迭代求解
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, m);
        Eigen::MatrixXd A = I + k*D;
        for(int i=1; i<=step; i++){
            u_ = A*u_;
        }
    };
};

class one_stage_Runge_Kutta : public Solver {
private:
    int step;
public:
    one_stage_Runge_Kutta(int step): step(step){};
    void solve(double r) override {
        if(r == 0)
            r = 1.0 / (2 * h_);
        double k = r * h_ * h_;
        construct_D(h_);

        u_ = Eigen::VectorXd::Zero(m);
        for(int i = 1; i <= m; i++)
            u_(i-1) = I_C(i * h_);

        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, m);
        Eigen::MatrixXd A = I - 0.5 * k * D;
        Eigen::MatrixXd B = I + 0.5 * k * D;

        for(int i = 1; i <= step; i++) {
            u_ = A.fullPivLu().solve(B * u_);
        }
    }
};

inline std::unique_ptr<Solver> createBackwardEulerSolver(int step){
    return std::unique_ptr<Solver>(new Backward_Euler(step));
}

inline std::unique_ptr<Solver> createCrankNicolsonSolver(int step){
    return std::unique_ptr<Solver>(new Crank_Nicolson(step));
}

inline std::unique_ptr<Solver> createForwardEulerSolver(int step){
    return std::unique_ptr<Solver>(new Forward_Euler(step));
}

inline std::unique_ptr<Solver> createOneStageRungeKuttaSolver(int step){
    return std::unique_ptr<Solver>(new one_stage_Runge_Kutta(step));
}

inline std::unique_ptr<Solver> creatCollocationSolver(int step){
    return std::unique_ptr<Solver>(new GLRKSolver(step));
}