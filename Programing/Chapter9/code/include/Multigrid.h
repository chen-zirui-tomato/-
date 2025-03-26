#pragma once
#include <map>
#include "Function.h"
#include <vector>
#include <cmath>

template<int Dim>
class Multigrid{
private:
    //网格大小为1/int时对应的矩阵和右端项
    std::map<int, std::pair<Matrix&, Vector&>> discretors;

    Vector(*Restriction)(int m, Vector& v);
    Vector(*Interpolation)(int m, Vector& v);
    Vector(*Cycle)(int m, Matrix& A, Vector& b, Vector& v0, int nu1, int nu2);
    Vector sol;
public:
    Multigrid(Function& f, Function& g, BCCond_t B):discretors(f, g, B);
    void Solve(int m, Restriction_t R, Interpolation_t I, Cycle_t C, int nu1, int nu2, double eps);
};

template<int Dim>
Multigrid<Dim>::Multigrid(Function& f, Function& g, BCCond_t B):discretors(f, g, B){
    Restriction = Restrict;
    Interpolation = Interpolate;
    Cycle = Cycle;
}

template<int Dim>
void Multigrid<Dim>::Solve(int max_levels, int nu1, int nu2, double tol){
    sol = Vector(discretors[0].second.size());
    for(int i=0; i<max_levels; i++){
        Vector v0 = discretors[i].second;
        Matrix A = discretors[i].first;
        Vector b = discretors[i].second;
        Vector v = Cycle(i, A, b, v0, nu1, nu2);
        discretors[i+1].second = v;
        discretors[i+1].first = A;
    }
}