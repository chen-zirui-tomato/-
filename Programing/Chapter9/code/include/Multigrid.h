#pragma once
#include <map>
#include "Function.h"
#include "Sparse.h"
#include "Strategy.h"
#include <vector>
#include <cmath>

template<int Dim>
class MultigridBase{
protected:

    //h = 1/n
    DynamicFunction& f,g;
    //网格大小为1/int时对应的矩阵和右端项
    std::map<int, std::pair<sparseMatrix, sparseVector>> discretors;
    ConditionType conditionType;

    // sparseVector(*Restriction)(int , sparseVector& );
    // sparseVector(*Interpolation)(int , sparseVector& );
    // sparseVector(*Cycle)(int , sparseMatrix& A, sparseVector& b, sparseVector& v0, int nu1, int nu2);

    std::unique_ptr<IRestriction> restriction_;
    std::unique_ptr<IInterpolation> interpolation_;
    std::unique_ptr<ICycle> cycle_;

    void setRestriction(std::unique_ptr<IRestriction>& restriction);
    void setInterpolation(std::unique_ptr<IInterpolation>& interpolation);
    void setCycle(std::unique_ptr<ICycle>& cycle);

    //真解
    sparseVector sol;
public:
    MultigridBase(DynamicFunction& f, DynamicFunction& g, ConditionType B):f(f), g(g), conditionType(B){};
    virtual void generateSparseMatrix(int , sparseMatrix&, sparseVector&) = 0;
    //n是网格大小，R是限制类型，I是插值类型，C是循环类型，nu1和nu2是循环的迭代次数，eps是停止精度
    virtual void Solver(int n, std::unique_ptr<IRestriction>&, std::unique_ptr<IInterpolation>&, std::unique_ptr<ICycle>&, int nu1, int nu2, double eps);
};

template<int Dim>
void MultigridBase<Dim>::setRestriction(std::unique_ptr<IRestriction>& restriction){
    restriction_ = std::move(restriction);
}
template<int Dim>
void MultigridBase<Dim>::setInterpolation(std::unique_ptr<IInterpolation>& interpolation){
    interpolation_ = std::move(interpolation);
}

template<int Dim>
void MultigridBase<Dim>::setCycle(std::unique_ptr<ICycle>& cycle){
    cycle_ = std::move(cycle);
}


//强制特化
template<int Dim>
class Multigrid;

template<>
class Multigrid<1> : public MultigridBase<1>{
private:
    struct MatrixGenerator{
        static void generateGridEquations(Multigrid<1>& mg, sparseMatrix& A, sparseVector&b, int n){
            A = sparseMatrix(n-1, n-1);
            b = sparseVector(n-1);
            for(int i = 1; i < n; i++)
                b.set_value(i, mg.f(i,i));
            for(int i = 2; i < n-1; i++) for(int index = i-1; index < i+2; index++)
                if(index-i+2 % 2 != 0) A.set_value(i, index, -1);
                else A.set_value(i, index, 2);
            switch (mg.conditionType)
            {
            case ConditionType::Dirichlet:
                constructDirichlet(mg, A, b, n);
                break;
            case ConditionType::Neumann:
                constructNeumann(mg, A, b, n);
                break;
            default:
                throw std::invalid_argument("Invalid condition type");
                break;
            }
        }
        static void constructDirichlet(Multigrid<1>& mg, sparseMatrix& A, sparseVector& b, int n){
            A.set_value(1, 1, 2), A.set_value(n-1, n-1, 2);
            b.set_value(1, n*n*mg.g(0,0)+mg.f(0,0)), b.set_value(n-1, n*n*mg.g(n,n)+mg.f(n,n));
        }
        static void constructNeumann(Multigrid<1>& mg, sparseMatrix& A, sparseVector& b, int n){
            A.set_value(1, 1, 1), A.set_value(n-1, n-1, 1);
            b.set_value(1, mg.f(0,0)-n*mg.g(0,0)), b.set_value(n-1, mg.f(n,n)-n*mg.g(n,n));
        }
    };
public:
    Multigrid(DynamicFunction& f, DynamicFunction& g, ConditionType B):MultigridBase<1>(f, g, B){};
    void generateSparseMatrix(int , sparseMatrix&, sparseVector&) override;
    void Solver(int m,  std::unique_ptr<IRestriction>&, std::unique_ptr<IInterpolation>&, std::unique_ptr<ICycle>&, int nu1, int nu2, double eps) override;
};


void Multigrid<1>::generateSparseMatrix(int m, sparseMatrix& A, sparseVector& b){
    MatrixGenerator::generateGridEquations(*this, A, b, m);
    discretors[m] = std::make_pair(A, b);
}

template<int Dim>
void MultigridBase<Dim>::Solver(int m, 
    std::unique_ptr<IRestriction>& restriction, 
    std::unique_ptr<IInterpolation>& interpolation, 
    std::unique_ptr<ICycle>& cycle, int nu1, int nu2, double eps){

    setRestriction(restriction);
    setInterpolation(interpolation);
    setCycle(cycle);

    double w = 2.0/3.0;
    sparseMatrix A;
    sparseVector b, v0;
    generateSparseMatrix(m, A, b);
    sparseMatrix I(m-1, m-1);
    for(int i = 1; i < m; i++) I.set_value(i, i, 1);
    sparseMatrix Tw = I - A*w/m/m/2;




    
    if(restriction_){
        sparseVector fine;
        sparseVector coarse = (*restriction_)(m, fine);
    }
    if(interpolation_){
        sparseVector coarse;
        sparseVector fine = (*interpolation_)(m, coarse);
    }
    if(cycle_){
        sparseVector solution = (*cycle_)(m, A, b, v0, nu1, nu2);
    }

}

