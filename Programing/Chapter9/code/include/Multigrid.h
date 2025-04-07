#pragma once
#include "Function.h"
#include "Sparse.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <vector>

template<int Dim>
class Multigrid{
public:

    // h = 1/n ， A.size() = n-1
    double h;
    DynamicFunction& f;
    DynamicFunction& g;
    DynamicFunction& exact;
    //网格大小为1/int时对应的矩阵和右端项, 即int-1阶矩阵和右端项
    std::map<int, std::pair<sparseMatrix, sparseVector>> discretors;
    ConditionType conditionType;

    // sparseVector(*Restriction)(int , sparseVector& );
    // sparseVector(*Interpolation)(int , sparseVector& );
    // sparseVector(*Cycle)(int , sparseMatrix& A, sparseVector& b, sparseVector& v0, int nu1, int nu2);

    using IRestrictionFunc = std::function<sparseVector(int, const sparseVector&)>;
    using IInterpolationFunc = std::function<sparseVector(int, const sparseVector&)>;
    using ICycleFunc = std::function<sparseVector(int, sparseMatrix&, sparseVector&, sparseVector&, int nu1, int nu2,
                                              std::string, std::string)>;

    //真解
    sparseVector sol;

    Multigrid(DynamicFunction& f, DynamicFunction& g, std::string B, DynamicFunction& exact)
    :f(f), g(g), exact(exact){
        if(B == "Dirichlet") conditionType = ConditionType::Dirichlet;
        else if(B == "Neumann") conditionType = ConditionType::Neumann;
        else throw std::invalid_argument("Invalid boundary condition");
    };

    virtual void generateSparseMatrix(int level) = 0;

    sparseVector relaxation(int level, int nu1, sparseVector v0, sparseVector RHS);


    //level = n;
    //加权平均, v0.size() = n-1;
    sparseVector fullWeighting(int level, const sparseVector& v0){
        if(v0.size() != level - 1) throw std::invalid_argument("Invalid fullWeighting input size of v0");

        sparseVector v(level/2-1);
        // v.set_value(1, 1.0/4*(f(0,0) + 2*v0(1) + v0(2)));
        // v.set_value(level/2-1, 1.0/4*(v0(level-2) + 2*v0(level-1) + f(1,1)));
        for(int j = 1; j <= level/2-1; j++) 
            v.set_value(j, 1.0/4*(v0(2*j-1) + 2*v0(2*j) + v0(2*j+1)));
        return v;
    };
    //取偶
    sparseVector injection(int level, const sparseVector& v0){
        if(v0.size() != level-1) throw std::invalid_argument("Invalid injection input size of v0");

        sparseVector v(level/2-1);
        for(int j = 1; j <= level/2-1; j++)
            v.set_value(j, v0(2*j));
        return v;
    };

    //level = n/2;
    //线性插值
    sparseVector linear(int level, const sparseVector& v0){
        if(v0.size() != level-1) throw std::invalid_argument("Invalid linear input size of v0");

        //v0.size = level - 1;
        sparseVector v(level*2 - 1);
        for(int j = 1; j <= level*2 - 1; j++)
            if(j % 2 == 0) v.set_value(j, v0(j/2));
            else if(j == 1) v.set_value(j, 3.0/2*v0(1) - 1.0/2*v0(2));
            else if(j == level*2 - 1) v.set_value(j, -1.0/2*v0(level-2) + 3.0/2*v0(level-1));
            else v.set_value(j, 1.0/2*(v0((j-1)/2) + v0((j+1)/2)));
        return v;
    };
    //二次插值
    sparseVector quadratic(int level, const sparseVector& v0){
        if(v0.size() != level-1) throw std::invalid_argument("Invalid quadratic input size of v0");

        sparseVector v(level*2 - 1);
        for(int j = 1; j <= level*2 - 1; j++)
            if(j % 2 == 0) v.set_value(j, v0(j/2));
            else if(j == 1) v.set_value(j, 1.0/8*(3*v0(1) + 8*v0(2) - v0(3)));
            else if(j == level*2 - 1) v.set_value(j, 1.0/8*(-1*v0(level - 3) + 6*v0(level - 2) + 3*v0(level - 1)));
            else if(j == 3) v.set_value(3, 1.0/16*(9*(v0((j-1)/2) + v0((j+1)/2)) - (-1.0/2*v0(1) + 3.0/2*v0(2))));
            else if(j == level*2 - 3) v.set_value(level*2 - 3, 1.0/16*(9*(v0((j-1)/2) + v0((j+1)/2)) - (3.0/2*v0(level-2) - 1.0/2*v0(level-1))));
            else v.set_value(j, 1.0/16*(9*(v0((j-1)/2) + v0((j+1)/2)) - (v0((j-1)/2-1) + v0((j+1)/2+1))));
        return v;
    };


    //第一次V-Cycle，level = n = finestGridSize
    sparseVector VCycle(int level, sparseVector vh, sparseVector &RHS_level, int nu1, int nu2, IRestrictionFunc IRes, IInterpolationFunc IInter){
        //VC-1的实现
        generateSparseMatrix(level);
        vh = relaxation(level, nu1, vh, RHS_level);
        //接下来是VC-2的实现
        //f^2hu
        if(level != 4){
            sparseVector RHS(level/2-1);
            sparseVector residual = (RHS_level - (discretors[level].first*vh));
            std::cerr << "Residual norm at level " << level << ": " << residual.norm() << std::endl;
            RHS = IRes(level, residual);
            sparseVector v_2h(level/2-1);
            v_2h = VCycle(level/2, v_2h, RHS, nu1, nu2, IRes, IInter);
            vh = vh + IInter(level/2, v_2h);
            vh = relaxation(level, nu2, vh, RHS_level);
        }
        
        else {
            // Direct solver for coarsest grid (level==4)
            try {
                Eigen::MatrixXd A = discretors[level].first.toDense();
                Eigen::VectorXd b = RHS_level.toDense();
                Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
                
                vh.resize(x.size());
                for(int i = 0; i < x.size(); ++i) {
                    if(std::abs(x(i)) > 1e-10) {
                        vh.set_value(i+1, x(i));
                    }
                }
            } catch(...) {
                // Fall back to relaxation if direct solver fails
                std::cerr << "Direct solver failed, falling back to relaxation" << std::endl;
                vh = relaxation(level, nu1, vh, RHS_level);
            }
        } 
        return vh;
    };


    sparseVector FMGCycle(int level, sparseVector &RHS, int nu1, int nu2, IRestrictionFunc IRes, IInterpolationFunc IInter){
        sparseVector vh(level-1);
        //FMG-1
        if(level != 4){
            generateSparseMatrix(level);
            sparseVector RHS_level = IRes(level, RHS);
            sparseVector v2h = FMGCycle(level/2, RHS_level, nu1, nu2, IRes, IInter);
            //FMG-2
            vh = IInter(level/2, v2h);
        }
        //FMG-3
        vh = VCycle(level, vh, RHS, nu1, nu2, IRes, IInter);
        RHS = (RHS - (discretors[level].first*vh));
        return vh;
    };

    //n是网格大小，R是限制类型，I是插值类型，C是循环类型，nu1和nu2是循环的迭代次数，eps是停止精度
    virtual void Solver(int m, std::string restriction, std::string interpolation, std::string cycle, int nu1, int nu2, double eps, int maxIter){
        IRestrictionFunc IRes;
        IInterpolationFunc IInter;
        
        if(restriction == "fullWeighting") 
            IRes = [this](auto && PH1, auto && PH2) { 
                return fullWeighting(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2)); 
            };
        else if(restriction == "injection") 
            IRes = [this](auto && PH1, auto && PH2) { 
                return injection(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2)); 
            };
        else throw std::invalid_argument("Invalid restriction type");

        if(interpolation == "linear") 
            IInter = [this](auto && PH1, auto && PH2) { 
                return linear(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2)); 
            };
        else if(interpolation == "quadratic") 
            IInter = [this](auto && PH1, auto && PH2) { 
                return quadratic(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2)); 
            };
        else throw std::invalid_argument("Invalid interpolation type");
        
        //初始猜测0
        sol.resize(m-1);
        // v_IG.set_value(m-1, 0);
        sparseVector RHS_level(m-1);
    
        if(cycle == "VCycle"){
            int count = 0;
            generateSparseMatrix(m);
            RHS_level = discretors[m].second;
            do {
                sol = VCycle(m, sol, RHS_level, nu1, nu2, IRes, IInter);
                std::cerr << ++count << " iteration, error = " << errorAnalysis() << std::endl;
            }while (errorAnalysis() > eps && count <= maxIter);
        }
        else if(cycle == "FMGCycle"){
            int count = 0;
            generateSparseMatrix(m);
            RHS_level = discretors[m].second;
            do {
                sol = FMGCycle(m, RHS_level, nu1, nu2, IRes, IInter);
                std::cerr << ++count << " iteration, error = " << errorAnalysis() << std::endl;
            }while (errorAnalysis() > eps && count <= maxIter);
        }
    };

    double errorAnalysis(){
        sparseVector v(this->sol.size());
        double h = 1.0/(this->sol.size()+1);
        for(int i = 1; i <= v.size(); i++)
            v.set_value(i, exact(i*h,i*h));
        double err = 0;
        for(int i = 1; i <= v.size(); i++)
            err += pow(this->sol(i) - v(i), 2);
        return sqrt(err);
    };
};

template<int Dim>
sparseVector Multigrid<Dim>::relaxation(int level,int nu1, sparseVector v0, sparseVector RHS){
    h = 1.0/level;
    // 根据网格层数动态调整松弛权重
    double w = (level <= 8) ? 0.8 : 2.0/3;
    // m = n-1
    if(discretors[level].first.size() != level-1) 
        throw std::invalid_argument("Matrix size not match with level");
    int m = discretors[level].first.size();
    
    sparseMatrix I(m, m);
    for(int i = 1; i <= m; i++) I.set_value(i, i, 1.0);
    sparseMatrix Tw = I - (discretors[level].first*w*h*h/2);
    // std::cerr << "Tw at level " << level << " = " << std::endl;
    // Tw.print();
    // std::cerr << std::endl;
    sparseVector current = RHS*w*h*h/2;
    // std::cerr << "current at level = "<< level << " = " << std::endl;
    // current.print();
    // std::cerr << std::endl;
    
    for(int i = 0; i < nu1; i++) {
        v0 = (Tw*v0) + current;
        // std::cerr << "v0 at level " << level << " after " << i+1 << " iteration = " ;
        // v0.print();
        // std::cerr << std::endl;
    }
    return v0;
}


// //强制特化
// template<int Dim>
// class Multigrid;

template<>
class Multigrid<1> : public Multigrid<0> {
private:
public:

    Multigrid<1>(DynamicFunction& f, DynamicFunction& g, std::string B, DynamicFunction& exact):Multigrid<0>(f, g, B, exact){};

    void constructDirichlet(sparseMatrix& A, sparseVector& b, int n){
        A.set_value(1, 1, 2.0), A.set_value(n-1, n-1, 2.0);
        A.set_value(n-1, n-2, -1.0),A.set_value(1, 2, -1.0);
        b.set_value(1, n*n*exact(0,0)+f(0,0)), 
        b.set_value(n-1, n*n*exact(1,1)+f(1,1));
    }
    void constructNeumann(sparseMatrix& A, sparseVector& b, int n){
        A.set_value(1, 1, 1.0), A.set_value(n-1, n-2, -1.0);
        A.set_value(1, 2, -1.0), A.set_value(n-1, n-1, 1.0);
        b.set_value(1, f(0,0)-n*g(0,0)), 
        b.set_value(n-1, f(1,1)-n*g(1,1));
    }

    //level = n； 生成n-1阶矩阵和右端项
    void generateSparseMatrix(int level) override;
    // void Solver(int m,  std::string, std::string, std::string, int nu1, int nu2, double eps) override;
};

void Multigrid<1>::generateSparseMatrix(int level){
    int m = level;
    sparseMatrix A(m-1, m-1);
    sparseVector b(m-1);
    double h = 1.0/m;
    int indexb = 1;
    for(double i = 2; i <= m-2; i ++)
        b.set_value(++indexb, f(i*h,i*h));
    b.tidyUp();


    for(int i = 2; i <= m-2; i++) for(int index = i-1; index <= i+1; index++)
        if((index-i+2) % 2 != 0) A.set_value(i, index, -1.0);
        else A.set_value(i, index, 2.0);
    switch (conditionType)
    {
    case ConditionType::Dirichlet:
        constructDirichlet(A, b, m);
        break;
    case ConditionType::Neumann:
        constructNeumann(A, b, m);
        break;
    default:
        throw std::invalid_argument("Invalid condition type");
        break;
    }
    A.tidyUp();
    A = A/h/h;

    discretors[m] = std::make_pair(A, b);
}











//在二维里，所有的level等于m，而我们要生成的是（m-1）*（m-1）的矩阵和右端项
template<>
class Multigrid<2> : public Multigrid<0> {
private:
public:
    std::vector<double> x_points;
    std::vector<double> y_points;
    std::vector<std::vector<double>> f_values;
    std::vector<std::vector<double>> g_values;
    // construct condition matrix as follow
    // | g_01    g_02   ... g_0,m-1 |
    // | g_10    g_20   ... g_m-1,0 |
    // | g_m,1   g_m,2  ... g_m,m-1 |
    // | g_1,m   g_2,m  ... g_m-1,m |
    DynamicFunction& gx;
    DynamicFunction& gy;

    Multigrid<2>(DynamicFunction& f, DynamicFunction& g, DynamicFunction& gx, DynamicFunction& gy, 
                 std::string B, DynamicFunction& exact):Multigrid<0>(f, g, B, exact), gx(gx), gy(gy){};

    sparseVector relaxation2D(int level, int nu1, sparseVector v0, sparseVector RHS) {
        h = 1.0/level;
        // 根据网格层数动态调整松弛权重
        double w = (level <= 8) ? 0.8 : 2.0/3;
        // m = n-1
        if(discretors[level].first.size() != (level-1)*(level-1)) 
            throw std::invalid_argument("Matrix size not match with level");
        int m = level;
        
        sparseMatrix I((m-1)*(m-1), (m-1)*(m-1));
        for(int i = 1; i <= (m-1)*(m-1); i++) 
            I.set_value(i, i, 1.0);
        sparseMatrix Tw = I - (discretors[level].first*w*h*h/4);
        
        sparseVector current = RHS*w*h*h/4;
        for(int i = 0; i < nu1; i++) {
            v0 = (Tw*v0) + current;
            // std::cerr << "v0 at level " << level << " after " << i+1 << " iteration = " ;
            // v0.print();
            // std::cerr << std::endl;
        }
        return v0;
    }

    void constructDirichlet2D(sparseMatrix& A, int m){
        // 边界
        for(int k = 2; k <= m-2; k++){
            f_values[1][k] += g_values[0][k-1],
            f_values[m-1][k] += g_values[2][k-1],
            f_values[k][1] += g_values[1][k-1],
            f_values[k][m-1] += g_values[3][k-1];
        }
        // 角点
        f_values[1][1] += (g_values[0][0] + g_values[1][0]),
        f_values[m-1][1] += (g_values[2][0] + g_values[1][m-2]),
        f_values[1][m-1] += (g_values[0][m-2] + g_values[3][0]),
        f_values[m-1][m-1] += (g_values[2][m-2] + g_values[3][m-2]);
    } 
    void constructNeumann2D(sparseMatrix& A, int m){
        for(int k = 2; k <= m-2; k++)
        A.set_value((k-1)*(m-1)+1, (k-1)*(m-1)+1, 3.0), f_values[1][k] += -g_values[0][k-1]*h/h/h,
        A.set_value(k*(m-1), k*(m-1), 3.0), f_values[m-1][k] += g_values[2][k-1]*h/h/h,
        A.set_value(k, k, 3.0), f_values[k][1] += -g_values[1][k-1]*h/h/h,
        A.set_value((m-2)*(m-1)+k, (m-2)*(m-1)+k, 3.0), f_values[k][m-1] += g_values[3][k-1]*h/h/h;

        // 角点
        std::cerr<<g_values[0][0]<<g_values[1][0]<<std::endl;
        A.set_value(1, 1, 2.0), f_values[1][1] += (-g_values[0][0]*h - g_values[1][0]*h),
        A.set_value(m-1, m-1, 2.0), f_values[m-1][1] += (g_values[2][0]*h - g_values[1][m-1]*h)/h/h,
        A.set_value((m-2)*(m-1)+1, (m-2)*(m-1)+1, 2.0), f_values[1][m-1] += (-g_values[0][m-1]*h + g_values[3][0]*h)/h/h,
        A.set_value((m-1)*(m-1), (m-1)*(m-1), 2.0), f_values[m-1][m-1] += (g_values[2][m-1]*h + g_values[3][m-1]*h)/h/h;
    }

    sparseVector fullWeighting2D(int level, const sparseVector& v_fine) {
        // std::cerr << "Enter fullWeighting2D" << std::endl;
        // v_fine.size() = (level-1)*(level-1)
        int coarse_size = level/2 - 1;
        sparseVector v_coarse(coarse_size * coarse_size);
        
        // 二维加权平均限制
        for(int j = 1; j <= coarse_size; j++) {
            for(int i = 1; i <= coarse_size; i++) {
                double sum = 0.0;
                // 3x3邻域加权
                sum += 1.0 * v_fine((2*j-2)*(level-1) + 2*i-1);
                sum += 2.0 * v_fine((2*j-2)*(level-1) + 2*i);
                sum += 1.0 * v_fine((2*j-2)*(level-1) + 2*i+1);
                
                sum += 2.0 * v_fine((2*j-1)*(level-1) + 2*i-1);
                sum += 4.0 * v_fine((2*j-1)*(level-1) + 2*i);
                sum += 2.0 * v_fine((2*j-1)*(level-1) + 2*i+1);
                
                sum += 1.0 * v_fine((2*j)*(level-1) + 2*i-1);
                sum += 2.0 * v_fine((2*j)*(level-1) + 2*i);
                sum += 1.0 * v_fine((2*j)*(level-1) + 2*i+1);
                
                v_coarse.set_value((j-1)*coarse_size + i, sum / 16.0);
            }
        }
        return v_coarse;
    }

    sparseVector injection2D(int level, const sparseVector& v_fine) {
        // v_fine.size() = (level-1)*(level-1)
        int coarse_size = level/2 - 1;
        sparseVector v_coarse(coarse_size * coarse_size);
        
        // 二维Injection
        for(int j = 1; j <= coarse_size; j++) {
            for(int i = 1; i <= coarse_size; i++) {
                v_coarse.set_value((j-1)*coarse_size + i, v_fine((2*j-1)*(level-1) + 2*i-1));
            }
        }
        return v_coarse;
    }

    sparseVector linear2D(int coarse_level, const sparseVector& v_coarse) {
        // std::cerr << "Enter linear2D" << std::endl;
        int fine_level = coarse_level * 2-1;
        coarse_level = coarse_level - 1;

        sparseVector v_fine(fine_level*fine_level);

        //加1是为了索引方便，00处都是0元，无效数据
        std::vector<std::vector<double>> 
        v_f(fine_level+1, std::vector<double>(fine_level+1, 0.0));
    
        std::vector<std::vector<double>>
        v_c(coarse_level+1, std::vector<double>(coarse_level+1, 0.0));

        for(int j = 1; j <= coarse_level; j++) {
            for(int i = 1; i <= coarse_level; i++) {
                v_c[i][j] = v_coarse((j-1)*(coarse_level) + i);
            }
        }

        for(int j = 2; j <= fine_level-1; j++)
            for(int i = 2; i <= fine_level-1; i++) {
                if(i % 2 == 0 && j % 2 == 0)
                    v_f[i][j] = v_c[(i/2-1)][(j/2-1)];
                else if(i % 2 == 0 && j % 2 != 0)
                    v_f[i][j] = 1.0/2*(v_c[(i/2-1)][(j-1)/2] + v_c[(i/2-1)][(j+1)/2]);
                else if(i % 2 != 0 && j % 2 == 0)   
                    v_f[i][j] = 1.0/2*(v_c[(i-1)/2][(j/2-1)] + v_c[(i+1)/2][(j/2-1)]);
                else if(i % 2 != 0 && j % 2 != 0)
                    v_f[i][j] = (v_c[(i-1)/2][(j-1)/2] + v_c[(i-1)/2][(j+1)/2] + v_c[(i+1)/2][(j-1)/2] + v_c[(i+1)/2][(j+1)/2]) / 4.0;
            }
        for(int k = 2; k <= fine_level-1; k++){
            if(k % 2 == 0){
                v_f[1][k] = 3.0/2*v_c[1][k/2] - 1.0/2*v_c[2][(k/2)];
                v_f[k][1] = 3.0/2*v_c[k/2][1] - 1.0/2*v_c[k/2][2];
                v_f[k][fine_level] = 3.0/2*v_c[k/2][coarse_level] - 1.0/2*v_c[(k/2)][coarse_level-1];
                v_f[fine_level][k] = 3.0/2*v_c[coarse_level][k/2] - 1.0/2*v_c[coarse_level-1][k/2];
            }
            else //if(k % 2 != 0)
                v_f[1][k] = 3.0/4*(v_c[1][(k-1)/2] + v_c[1][(k+1)/2]) - 1.0/4*(v_c[2][(k-1)/2]+v_c[2][(k+1)/2]);
                v_f[k][1] = 3.0/4*(v_c[(k-1)/2][1] + v_c[(k+1)/2][1]) - 1.0/4*(v_c[(k-1)/2][2]+v_c[(k+1)/2][2]);
                v_f[k][fine_level] = 3.0/4*(v_c[(k-1)/2][coarse_level] + v_c[(k+1)/2][coarse_level]) - 1.0/4*(v_c[(k-1)/2][coarse_level-1]+v_c[(k+1)/2][coarse_level-1]);
                v_f[fine_level][k] = 3.0/4*(v_c[coarse_level][(k-1)/2] + v_c[coarse_level][(k+1)/2]) - 1.0/4*(v_c[coarse_level-1][(k-1)/2]+v_c[coarse_level-1][(k+1)/2]);
        }
        v_f[1][1] = 3.0/2*v_c[1][1] - 1.0/2*v_c[2][2];
        v_f[fine_level][1] = 3.0/2*v_c[coarse_level][1] - 1.0/2*v_c[coarse_level-1][2];
        v_f[1][fine_level] = 3.0/2*v_c[1][coarse_level] - 1.0/2*v_c[2][coarse_level-1];
        v_f[fine_level][fine_level] = 3.0/2*v_c[coarse_level][coarse_level] - 1.0/2*v_c[coarse_level-1][coarse_level-1];

        for(int j = 1; j <= fine_level; j++) {
            for(int i = 1; i <= fine_level; i++) {
                v_fine.set_value((j-1)*(fine_level) + i, v_f[i][j]);
            }
        }
        return v_fine;
    }

    //level = n； 生成n-1阶矩阵和右端项
    void generateSparseMatrix(int level) override;

    sparseVector VCycle2D(int level, sparseVector vh, sparseVector RHS_level, 
        int nu1, int nu2, IRestrictionFunc IRes, IInterpolationFunc IInter) {
        generateSparseMatrix(level);
        vh = relaxation2D(level, nu1, vh, RHS_level);

        if (level > 4) {

            sparseVector residual = RHS_level - discretors[level].first * vh;
            std::cerr << "Residual norm at level " << level << ": " << residual.norm() << std::endl;
            sparseVector RHS_coarse = IRes(level, residual);
            sparseVector v_2h((level/2-1)*(level/2-1));

            v_2h = VCycle2D(level/2, v_2h, RHS_coarse, nu1, nu2, IRes, IInter);
            vh = vh + IInter(level/2, v_2h);
            vh = relaxation2D(level, nu2, vh, RHS_level);
        }
        else {
            try {
                Eigen::MatrixXd A = discretors[level].first.toDense();
                Eigen::VectorXd b = RHS_level.toDense();
                Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
                
                vh.resize(x.size());
                for(int i = 0; i < x.size(); ++i) {
                    if(std::abs(x(i)) > 1e-10) {
                        vh.set_value(i+1, x(i));
                    }
                }
            } catch(...) {
                std::cerr << "Direct solver failed, falling back to relaxation" << std::endl;
                vh = relaxation(level, nu1, vh, RHS_level);
            }
        } 
        return vh;
    }

    sparseVector FMGCycle2D(int level, sparseVector &RHS, int nu1, int nu2, IRestrictionFunc IRes, IInterpolationFunc IInter){
        sparseVector vh((level-1)*(level-1));
        //FMG-1
        if(level != 4){
            generateSparseMatrix(level);
            sparseVector RHS_level = IRes(level, RHS);
            sparseVector v2h = FMGCycle2D(level/2, RHS_level, nu1, nu2, IRes, IInter);
            //FMG-2
            vh = IInter(level/2, v2h);
        }
        //FMG-3
        vh = VCycle2D(level, vh, RHS, nu1, nu2, IRes, IInter);
        RHS = (RHS - (discretors[level].first*vh));
        return vh;
    };

    void Solver(int m,  std::string restriction, std::string interpolation, std::string cycle, 
                int nu1, int nu2, double eps, int maxIter) override{

        IRestrictionFunc IRes;
        IInterpolationFunc IInter;
        if(restriction == "fullWeighting") 
            IRes = [this](auto && PH1, auto && PH2) { 
                return fullWeighting2D(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2)); 
            };
        else if(restriction == "injection") 
            IRes = [this](auto && PH1, auto && PH2) { 
                return injection2D(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2)); 
            };
        else throw std::invalid_argument("Invalid restriction type");

        if(interpolation == "linear") 
            IInter = [this](auto && PH1, auto && PH2) { 
                return linear2D(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2)); 
            };
        else if(interpolation == "quadratic") 
            IInter = [this](auto && PH1, auto && PH2) { 
                return quadratic(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2)); 
            };
        else throw std::invalid_argument("Invalid interpolation type");
        
        //初始猜测0
        sol.resize((m-1)*(m-1));
        // v_IG.set_value(m-1, 0);
        sparseVector RHS_level((m-1)*(m-1));

        if(cycle == "VCycle"){
            int count = 0;
            generateSparseMatrix(m);
            RHS_level = discretors[m].second;
            do {
                sol = VCycle2D(m, sol, RHS_level, nu1, nu2, IRes, IInter);
                std::cerr << ++count << " iteration, error = " << errorAnalysis2D(m) << std::endl;
            }while (errorAnalysis2D(m) > eps && count <= maxIter);
        }
        else if(cycle == "FMGCycle"){
            int count = 0;
            generateSparseMatrix(m);
            RHS_level = discretors[m].second;
            do {
                sol = FMGCycle2D(m, RHS_level, nu1, nu2, IRes, IInter);
                std::cerr << ++count << " iteration, error = " << errorAnalysis2D(m) << std::endl;
            }while (errorAnalysis2D(m) > eps && count <= maxIter);
        }
    }

    double errorAnalysis2D(int m){
        sparseVector v(this->sol.size());
        double h = 1.0/(this->sol.size()+1);
        for(int i = 1; i <= m-1; i++) for(int j = 1; j <= m-1; j++)
            v.set_value((j-1)*(m-1)+i, exact(i*h,j*h));
        double err = 0;
        for(int i = 1; i <= m-1; i++) for(int j = 1; j <= m-1; j++)
            err += pow(this->sol((i-1)*(m-1)+j) - v((i-1)*(m-1)+j), 2);
        return sqrt(err);
    }

};

//level = m, 生成(m-1)*(m-1)阶矩阵和右端项,但是仍然被记录在discretors[m]中.
void Multigrid<2>::generateSparseMatrix(int level){
    // std::cerr<<"Debug - generateSparseMatrix called with level = "<<level<<std::endl;
    int m = level;
    sparseMatrix A((m-1)*(m-1), (m-1)*(m-1));
    double h = 1.0/m;
    f_values.resize(m+1, std::vector<double>(m+1, 0.0));
    g_values.resize(4, std::vector<double>(m-1, 0.0));
    for(int i = 0; i <= m; i++) 
        x_points.push_back(i*h), y_points.push_back(i*h);
    for(int i = 0; i <= m; i++) for(int j = 0; j <= m; j++)
        f_values[i][j] = f(x_points[i], y_points[j]);
    for(int i = 0 ;i <= m-2; i++){
        g_values[0][i] = exact(x_points[0], y_points[i+1]);
        g_values[1][i] = exact(x_points[i+1], y_points[0]);
        g_values[2][i] = exact(x_points[m], y_points[i+1]);
        g_values[3][i] = exact(x_points[i+1], y_points[m]);
    }
    // std::cerr << "RHS generated with level = " << level << " =" << std::endl;
    // b.print();
    // std::cerr << std::endl;

    //内点
    for(int i = 2; i <= m-2; i++) for(int j = 2; j <= m-2; j++)
        A.set_value((j-1)*(m-1)+i, (j-1)*(m-1)+i, 4.0),
        A.set_value((j-1)*(m-1)+i, (j-1)*(m-1)+i-1, -1.0),
        A.set_value((j-1)*(m-1)+i, (j-1)*(m-1)+i+1, -1.0),
        A.set_value((j-1)*(m-1)+i, (j-2)*(m-1)+i, -1.0),
        A.set_value((j-1)*(m-1)+i, (j)*(m-1)+i, -1.0);
    //边界
    for(int k = 2; k <= m-2; k++)
        A.set_value((k-1)*(m-1)+1, (k-1)*(m-1)+1, 4.0),
        A.set_value((k-1)*(m-1)+1, (k-1)*(m-1)+2, -1.0),
        A.set_value((k-1)*(m-1)+1, (k-2)*(m-1)+1, -1.0),
        A.set_value((k-1)*(m-1)+1, (k)*(m-1)+1, -1.0),

        A.set_value(k*(m-1), k*(m-1), 4.0),
        A.set_value(k*(m-1), k*(m-1)-1, -1.0),
        A.set_value(k*(m-1), (k-1)*(m-1), -1.0),
        A.set_value(k*(m-1), (k+1)*(m-1), -1.0),

        A.set_value(k, k, 4.0),
        A.set_value(k, k-1, -1.0),
        A.set_value(k, k+1, -1.0),
        A.set_value(k, k+m-1, -1.0),

        A.set_value((m-2)*(m-1)+k, (m-2)*(m-1)+k, 4.0),
        A.set_value((m-2)*(m-1)+k, (m-2)*(m-1)+k-1, -1.0),
        A.set_value((m-2)*(m-1)+k, (m-2)*(m-1)+k+1, -1.0),
        A.set_value((m-2)*(m-1)+k, (m-3)*(m-1)+k, -1.0);

    //角点
    A.set_value(1, 1, 4.0), 
    A.set_value(1, 2, -1.0),
    A.set_value(1, m, -1.0);

    A.set_value(m-1, m-1, 4.0), 
    A.set_value(m-1, m-2, -1.0), 
    A.set_value(m-1, 2*(m-1), -1.0);

    A.set_value((m-2)*(m-1)+1, (m-2)*(m-1)+1, 4.0), 
    A.set_value((m-2)*(m-1)+1, (m-2)*(m-1)+2, -1.0), 
    A.set_value((m-2)*(m-1)+1, (m-3)*(m-1)+1, -1.0);

    A.set_value((m-1)*(m-1), (m-1)*(m-1), 4.0), 
    A.set_value((m-1)*(m-1), (m-1)*(m-1)-1, -1.0), 
    A.set_value((m-1)*(m-1), (m-2)*(m-1), -1.0);

    switch (conditionType)
    {
    case ConditionType::Dirichlet:
        constructDirichlet2D(A, m);
        break;
    case ConditionType::Neumann:
        constructNeumann2D(A, m);
        break;
    default:
        throw std::invalid_argument("Invalid condition type");
        break;
    }
    sparseVector b((m-1)*(m-1));
    int index = 0;
    for(int j = 1; j <= m-1; j++) for(int i = 1; i <= m-1; i++)
        b.set_value(++index, f_values[i][j]);
    b.tidyUp();
    A.tidyUp();
    A = A/h/h;
    // A.print();
    // std::cerr << "Matrix & RHS generated with level = " << level << " complete " << std::endl;

    discretors[m] = std::make_pair(A, b);
}