#pragma once
#include <map>
#include "Function.h"
#include "Sparse.h"
#include <functional>
#include <vector>
#include <cmath>

template<int Dim>
class Multigrid{
public:

    // h = 1/n ， A.size() = n-1
    double h;
    DynamicFunction& f,g,exact;
    //网格大小为1/int时对应的矩阵和右端项, 即int-1阶矩阵和右端项
    std::map<int, std::pair<sparseMatrix&, sparseVector&>> discretors;
    ConditionType conditionType;

    // sparseVector(*Restriction)(int , sparseVector& );
    // sparseVector(*Interpolation)(int , sparseVector& );
    // sparseVector(*Cycle)(int , sparseMatrix& A, sparseVector& b, sparseVector& v0, int nu1, int nu2);

    using IRestrictionFunc = std::function<sparseVector(int, sparseVector&)>;
    using IInterpolationFunc = std::function<sparseVector(int, sparseVector&)>;
    using ICycleFunc = std::function<sparseVector(int, sparseMatrix&, sparseVector&, sparseVector&, int nu1, int nu2,
                                              std::string, std::string)>;

    //真解
    sparseVector sol;
    //初始猜测
    sparseVector v_IG;

    Multigrid(DynamicFunction& f, DynamicFunction& g, std::string B, DynamicFunction& exact):f(f), g(g), exact(exact){
        if(B == "Dirichlet") conditionType = Dirichlet;
        else if(B == "Neumann") conditionType = Neumann;
        else throw std::invalid_argument("Invalid boundary condition");
    };

    virtual void generateSparseMatrix(int level) = 0;

    sparseVector relaxation(int level, int nu1, sparseVector& v0);


    //level = n;
    //加权平均, v0.size() = n-1;
    sparseVector fullWeighting(int level, sparseVector& v0){
        sparseVector v(level/2-1);
        v.set_value(1, 1.0/4*(f(0,0) + 2*v0(1) + v0(2)));
        v.set_value(level/2-1, 1.0/4*(v0(level-2) + 2*v0(level-1) + f(1,1)));
        for(int j = 2; j <= level/2-1; j++) 
            v.set_value(j, 1.0/4*(v0(2*j-2) + 2*v0(2*j-1) + v0(2*j)));
        return v;
    };
    //取偶
    sparseVector injection(int level, sparseVector& v0){
        sparseVector v(level/2-1);
        for(int j = 1; j <= level/2; j++)
            v.set_value(j, v0(2*j-1));
        return v;
    };

    //level = n/2;
    //线性插值
    sparseVector linear(int level, sparseVector& v0){
        //v0.size = level - 1;
        sparseVector v(level*2 - 1);
        for(int j = 1; j <= level*2 - 1; j++)
            if(j % 2 != 0) v.set_value(j, v0((j+1)/2));
            else v.set_value(j, 0.5*(v0(j/2) + v0(j/2+1)));
        return v;
    };
    //二次插值
    sparseVector quadratic(int level, sparseVector& v0){
        sparseVector v(level*2 - 1);
        for(int j = 1; j <= level*2 - 1; j++)
            if(j % 2 != 0) v.set_value(j, v0((j+1)/2));
            else if(j == 2) v.set_value(j, 1.0/8*(3*v0(1) + 8*v0(2) - v0(3)));
            else if(j == level*2 - 2) v.set_value(j, 1.0/8*(-1*v0(level*2 - 2) + 6*v0(level*2 - 1) + 3*0v0(level*2 - 3)));
            else v.set_value(j, 1.0/16*(9*(v0(j/2) + v0(j/2+1)) - (v0(j/2-1) + v0(j/2+2))));
        return v;
    };


    //第一次V-Cycle，level = n = finestGridSize
    sparseVector VCycle(int level, sparseVector& v0, int nu1, int nu2, std::string res, std::string inter){
        //VC-1的实现
        v0 = relaxation(level, nu1, v0);

        if(level != 4){
            IRestrictionFunc IRes;
            IInterpolationFunc IInter;
            
            if(res == "fullWeighting") IRes = fullWeighting;
            else if(res == "injection") IRes = injection;
            else throw std::invalid_argument("Invalid restriction type");

            if(inter == "linear") IInter = linear;
            else if(inter == "quadratic") IInter = quadratic;
            else throw std::invalid_argument("Invalid interpolation type");
            

            //接下来是VC-2的实现
            //f^2h
            sparseVector RHS(level/2-1) = IRes(level, discretors[level].second - discretors[level].first*v0);
            sparseVector v_2h(level/2-1);
            v_2h = VCycle(level/2, v_2h, nu1, nu2, res, inter);
            v0 = v0 + IInter(level/2, v_2h);
        }
        v0 = relaxation(level, nu2, v0);
        return v0;
    };


    sparseVector FMGCycle(int level, sparseVector& v);

    //n是网格大小，R是限制类型，I是插值类型，C是循环类型，nu1和nu2是循环的迭代次数，eps是停止精度
    virtual void Solver(int n, std::string, std::string, std::string, int nu1, int nu2, double eps);

    double errorAnalysis(){
        spareseVector v(sol.size());
        for(int i = 1; i <= v.size(); i++)
            v.set_value(i, exact(i,i));
        double err = 0;
        for(int i = 1; i <= v.size(); i++)
            err += pow(v(i)-sol(i), 2);
        return sqrt(err);
    };
};

template<int Dim>
sparseVector Multigrid<Dim>::relaxation(int level,int nu1, sparseVector& v0){
    h = 1.0/(m+1);
    double w = 2.0/3.0;
    // m = n-1
    int m = discretors[level].first.size();
    sparseMatrix I(m, m);
    for(int i = 1; i < m+1; i++) I.set_value(i, i, 1);
    sparseMatrix Tw = I - discretors[level].first*w*h*h/2;
    for(int i = 0; i < nu1; i++) v0 = Tw*v0;
    sparseVector current = discretors[level].second*w*h*h/2;
    for(int i = 0; i < nu1; i++) {
        v0 += current;
        current = Tw*current;
    }
    return v0;
}

template<int Dim>
void Multigrid<Dim>::Solver(int m, 
    std::string restriction, 
    std::string interpolation, 
    std::string cycle, int nu1, int nu2, double eps){
    
    generateSparseMatrix(m);
    v_IG.resize(m-1);
    v_IG.set_value(m-1, 0);

    if(cycle == "VCycle"){
        sol = ICycle(m, v_IG, nu1, nu2, restriction, interpolation);
    }
    else if(cycle == "FMGCycle"){
        sol = FMGCycle(m, v_IG, nu1, nu2, restriction, interpolation);
    }

};

// //强制特化
// template<int Dim>
// class Multigrid;

template<>
class Multigrid<1> : public Multigrid<0> {
    // 直接使用基类成员，无需显式引入
    // 编译器应自动继承所有public/protected成员
private:

public:

    void constructDirichlet(sparseMatrix& A, sparseVector& b, int n){
        A.set_value(1, 1, 2), A.set_value(n-1, n-1, 2);
        b.set_value(1, n*n*g(0,0)+f(0,0)), b.set_value(n-1, n*n*g(n,n)+f(n,n));
    }
    void constructNeumann(sparseMatrix& A, sparseVector& b, int n){
        A.set_value(1, 1, 1), A.set_value(n-1, n-1, 1);
        b.set_value(1, f(0,0)-n*g(0,0)), b.set_value(n-1, f(n,n)-n*g(n,n));
    }

    Multigrid(DynamicFunction& f, DynamicFunction& g, std::string B, DynamicFunction& exact):Multigrid<0>(f, g, B, exact){};
    //level = n； 生成n-1阶矩阵和右端项
    void generateSparseMatrix(int level) override;
    void Solver(int m,  std::string, std::string, std::string, int nu1, int nu2, double eps) override;
};

void Multigrid<1>::generateSparseMatrix(int level){
    int m = level;
    sparseMatrix A(m-1, m-1);
    sparseVector b(m-1);
    double h = 1.0/m;
    for(int i = h; i < 1; i += h)
        b.set_value(i, f(i,i));
    for(int i = 2; i < m-1; i++) for(int index = i-1; index < i+2; index++)
        if(index-i+2 % 2 != 0) A.set_value(i, index, -1);
        else A.set_value(i, index, 2);
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
    discretors[m] = std::make_pair(A, b);
}

