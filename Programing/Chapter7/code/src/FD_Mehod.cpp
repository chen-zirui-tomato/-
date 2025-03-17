#include "../include/FD_Method.h"

FD_Method::FD_Method(int m, DynamicFunction& f, DynamicFunction& gx, DynamicFunction& gy, 
    std::string condition_type): m(m), h(1.0/(m+1)), f(f), gx(gx), gy(gy)
    {   
        A.resize(m*m, m*m);
        A.setZero();
        f_values.resize(m+2, std::vector<double>(m+2, 0.0));
        g_values.resize(4, std::vector<double>(m, 0.0));
        for(int i=0; i < m+2; i++)
            x_points.push_back(i*h), y_points.push_back(i*h);
        for(int i=0; i < m+2; i++) for(int j=0; j < m+2; j++)
            f_values[i][j] = f(x_points[i], y_points[j])*h*h;
        for(int i = 1 ;i < m+1; i++){
            g_values[0][i] = -gy(x_points[0], y_points[i]);
            g_values[1][i] = -gx(x_points[i], y_points[0]);
            g_values[2][i] = gy(x_points[m+1], y_points[i]);
            g_values[3][i] = gx(x_points[i], y_points[m+1]);
        }
        if(condition_type == "Dirichlet")
            this->condition_type = ConditionType::Dirichlet;
        else if(condition_type == "Neumann")
            this->condition_type = ConditionType::Neumann;
    }

void FD_Method::construct_A_matrix(){
    // 1.内点
    for(int i = 2; i < m; i++) for(int j = 2; j < m; j++) //u索引
        A((j-1)*m+i-1, (j-1)*m+i-1) = 4,
        A((j-1)*m+i-1, (j-1)*m+i-2) = -1,
        A((j-1)*m+i-1, (j-1)*m+i) = -1, 
        A((j-1)*m+i-1, (j-2)*m+i-1) = -1,
        A((j-1)*m+i-1, j*m+i-1) = -1;
    // 2.边界
    for(int k = 2; k < m; k++) 
        A((k-1)*m, (k-1)*m) = 4,
        A((k-1)*m, (k-1)*m+1) = -1,
        A((k-1)*m, (k-2)*m) = -1,
        A((k-1)*m, k*m) = -1,

        A(k*m-1, k*m-1) = 4,
        A(k*m-1, k*m-2) = -1,
        A(k*m-1, (k-1)*m-1) = -1,
        A(k*m-1, (k+1)*m-1) = -1,

        A(k-1, k-1) = 4,
        A(k-1, k-2) = -1,
        A(k-1, k) = -1,
        A(k-1, m+k-1) = -1,

        A((m-1)*m+k-1, (m-1)*m+k-1) = 4,
        A((m-1)*m+k-1, (m-1)*m+k-2) = -1,
        A((m-1)*m+k-1, (m-2)*m+k-1) = -1,
        A((m-1)*m+k-1, (m-1)*m+k) = -1;

    // 3.角点
    A(0, 0) = 4, A(0, 1) = -1, A(0, m) = -1;
    A(m-1, m-1) = 4, A(m-1, m-2) = -1, A(m-1, 2*m-1) = -1;
    A((m-1)*m, (m-1)*m) = 4, A((m-1)*m, (m-1)*m+1) = -1, A((m-1)*m, (m-2)*m) = -1;
    A(m*m-1, m*m-1) = 4, A(m*m-1, m*m-2) = -1, A(m*m-1, (m-1)*m-1) = -1;

}

void FD_Method::construct_RHS_vector(){
    // 4.边界条件
    switch(condition_type){
        case ConditionType::Dirichlet:
            construct_D_condition();
            break;
        case ConditionType::Neumann:
            construct_N_condition();
            break;
        case ConditionType::Mixed:
            construct_M_condition();
            break;
    }
}

void FD_Method::construct_D_condition(){
    // 边界
    for(int k = 2; k < m; k++)
        f_values[1][k] += g_values[0][k-1],
        f_values[m][k] += g_values[2][k-1],
        f_values[k][1] += g_values[1][k-1],
        f_values[k][m] += g_values[3][k-1];
    // 角点
    f_values[1][1] += (g_values[0][0] + g_values[1][0]),
    f_values[m][1] += (g_values[2][0] + g_values[1][m-1]),
    f_values[1][m] += (g_values[0][m-1] + g_values[3][0]),
    f_values[m][m] += (g_values[2][m-1] + g_values[3][m-1]);
}

void FD_Method::construct_N_condition(){
    // 边界
    for(int k = 2; k < m; k++)
        A((k-1)*m, (k-1)*m) += -1, f_values[1][k] += -g_values[0][k-1]*h,
        A(k*m-1, k*m-1) += -1, f_values[m][k] += g_values[2][k-1]*h,
        A(k-1, k-1) += -1, f_values[k][1] += -g_values[1][k-1]*h,
        A((m-1)*m+k-1, (m-1)*m+k-1) += -1, f_values[k][m] += g_values[3][k-1]*h;

    // 角点
    A(0, 0) += -2, f_values[1][1] += -g_values[0][0]*h - g_values[1][0]*h,
    A(m-1, m-1) += -2, f_values[m][1] += g_values[2][0]*h - g_values[1][m-1]*h,
    A((m-1)*m, (m-1)*m) += -2, f_values[1][m] += -g_values[0][m-1]*h +g_values[3][0]*h,
    A(m*m-1, m*m-1) += -2, f_values[m][m] += g_values[2][m-1]*h + g_values[3][m-1]*h;
}

void FD_Method::construct_M_condition(){}

void FD_Method::construct_equation(){
    construct_A_matrix();
    construct_RHS_vector();
}

void FD_Method::solver(){
    f_eigen.resize(m*m);
    int index = -1;
    for(int j = 1; j < m+1; j++) for(int i = 1; i < m+1; i++)
        f_eigen(++index) = f_values[i][j];
    A *= 1/h/h;

    // for(int i = 0; i < m*m; i++){ for(int j = 0; j < m*m; j++){
    //     std::cout<<A(i,j);
    // }
    // std::cout<<std::endl; }
    // for(int i = 0; i < m*m; i++) std::cout<<f_eigen(i)<<std::endl;

    std::cerr<<"开始求解"<<std::endl;

    // Eigen::VectorXd u_eigen = solve_zero_friendly(A, f_eigen);
    // std::cout << "解向量:\n" << u_eigen << std::endl;
    u_eigen = A.fullPivLu().solve(f_eigen);

    std::cerr<<"求解完成"<<std::endl;
    std::ofstream outFile("../opt/data.txt");
    if (outFile.is_open()) {
        for (int i = 0; i < u_eigen.size(); ++i) {
            outFile << u_eigen(i);
            if ((i + 1) % m == 0 && i + 1 < u_eigen.size()) {
                outFile << "\n";
            } else if (i + 1 < u_eigen.size()) {
                outFile << " ";
            }
        }
        outFile.close();
    } else {
        std::cerr << "Unable to open file";
    }
}



Eigen::VectorXd FD_Method::solve_zero_friendly(const Eigen::MatrixXd& A, const Eigen::VectorXd& f) {
    // 使用全主元LU分解
    Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
    
    // 检查矩阵秩是否足够
    if (lu.rank() < A.cols()) {
        // 若有自由变量，强制自由变量为0
        Eigen::VectorXd x = Eigen::VectorXd::Zero(A.cols());
        x.head(lu.rank()) = lu.solve(f); // 只解非自由变量部分
        return x;
    } else {
        // 矩阵满秩，直接求解
        return lu.solve(f);
    }
}

void FD_Method::error_analysis(DynamicFunction& E){
    //真解E
    std::vector<double> E_values(m*m, 0.0);
    for(int j = 1; j < m+1; j++) for(int i = 1; i < m+1; i++)
        E_values.push_back(E(x_points[i], y_points[j]));
    //误差分析
    error_L2(E_values);
    error_L1(E_values);
    error_Linfty(E_values);

}

void FD_Method::error_L2(std::vector<double>& E_values){
    double error = 0.0;
    for(int i = 0; i < m*m; i++){
        error += pow(u_eigen(i) - E_values[i], 2);
    }
    error = sqrt(error*h);
    std::cout<<"L2误差为："<<error<<std::endl;
}

void FD_Method::error_L1(std::vector<double>& E_values){
    double error = 0.0;
    for(int i = 0; i < m*m; i++){
        error += fabs(u_eigen(i) - E_values[i]);
    }
    error = h*error;
    std::cout<<"L1误差为："<<error<<std::endl;
}

void FD_Method::error_Linfty(std::vector<double>& E_values){
    double error = 0.0;
    for(int i = 0; i < m*m; i++){
        error = (error > fabs(u_eigen(i) - E_values[i])) ? error : fabs(u_eigen(i) - E_values[i]);
    }
    std::cout<<"Linfty误差为："<<error<<std::endl;
}

















