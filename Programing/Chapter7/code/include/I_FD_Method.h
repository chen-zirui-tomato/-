#pragma once
#include "FD_Method.h"

class I_FD_Method : public FD_Method{
private:
    enum class Sign{
        up, down, left, right, null
    };
    double delta;
    double center[2];
    double radius;
    int count;
    DynamicFunction& centerFunc;
    Eigen::MatrixXd A_eigen_rebuild;
    Eigen::VectorXd f_eigen_rebuild;
    // 标记闭圆盘内部的点
    std::map<std::pair<int, int>, bool> mark_in;
    // 标记闭圆盘边界的点
    std::map<std::pair<int, int>, bool> mark_bound;
    

public:
    I_FD_Method(int m, DynamicFunction& f, DynamicFunction& b, DynamicFunction& gx, DynamicFunction& gy,
                DynamicFunction& centerFunc, 
                std::string condition_type, std::vector<double> center, double radius);

    void I_solver();

    void re_construct_equations();

    void re_construct_RHS(int i, int j, Sign sign1, Sign sign2);
    
    void judge(int i, int j, Sign sign);

    // 找到在sign方向上与u_ij相邻的ghost值并改变A矩阵和f_values
    void findBoundaryNeighbor(int i, int j, double k, int sign_, Sign sign);

    void build_Dirichlet_condition(int i, int j, double y, Sign sign);

    void build_Neumann_condition(int i, int j, double k, Sign sign);
};
