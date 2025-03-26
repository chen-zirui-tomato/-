#pragma once
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "Function.h"
#include <map>

class FD_Method{
protected:
    enum NormType{
        L2,
        L1,
        Linfty
    };
    int m;
    double h;
    std::vector<double> x_points;
    std::vector<double> y_points;
    //f是整个区域包括边界的函数
    std::vector<std::vector<double>> f_values;
    // std::vector<double> u_vector;
    // （u_11, u_21,..., u_m1, u_12, u_22,..., u_m2,..., u_1m, u_2m,..., u_mm）T
    DynamicFunction& f;
    DynamicFunction& b;
    DynamicFunction& gx;
    DynamicFunction& gy;
    //A只有区域内部
    Eigen::MatrixXd A;
    Eigen::VectorXd f_eigen;
    Eigen::VectorXd u_eigen;
    //g是边界函数
    std::vector<std::vector<double>> g_values_real;
    std::vector<std::vector<double>> g_values;
    // construct condition matrix as follow
    // | g_01    g_02   ... g_0m    |
    // | g_10    g_20   ... g_m0    |
    // | g_m+1,1 g_m+1,2... g_m+1,m |
    // | g_1,m+1 g_2,m+1... g_m,m+1 |

    //ghost vector G
    // std::vector<std::vector<double>> G;
    // construct ghost pixel
    // | G_-1,1  G_-1,2 ... G_-1,m |
    // | G_1,-1  G_2,-1 ... G_m,-1 |
    // | G_m+2,1 G_m+2,2... G_m+2,m |
    // | G_1,m+2 G_2,m+2... G_m,m+2 |
    ConditionType condition_type;

    void construct_A_matrix();

    void construct_RHS_vector();

    void construct_D_condition();

    void construct_M_condition();

    void construct_N_condition();

    double get_Length(double x[2], double y[2]){
        return sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]));
    };

    Eigen::VectorXd solve_zero_friendly(const Eigen::MatrixXd& A, const Eigen::VectorXd& f);

public:

    FD_Method(int m, DynamicFunction& f, DynamicFunction& b, DynamicFunction& gx, DynamicFunction& gy, 
            std::string condition_type);

    void construct_equation();

    void solver();

    void error_analysis(DynamicFunction& exact_solution);

    void error_L1(std::vector<double>& E);

    void error_L2(std::vector<double>& E);

    void error_Linfty(std::vector<double>& E);
};