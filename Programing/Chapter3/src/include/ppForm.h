#include <iostream>
#include "Function.h"
#include <vector>
#include <algorithm>
#include <cmath>

class ppForm{
private:
    const Function& f; //插值函数
    int nodeCount; //插值节点的数量, n = nodeCount - 1
    int n;
    std::vector<double> nodes; //插值节点的位置, nodes【0】= x0
    std::vector<double> values; //插值节点的函数值, values【0】= f(x0)
    SplineType splineType; //插值类型
    std::vector<double> M; //存储二阶导数值
    std::vector<double> miu; //存储miu值
    std::vector<double> nabuda; //存储nabuda值
    std::vector<double> b; //存储三阶差商值
    std::vector<std::vector<double>> A; //存储三次样条的三对角矩阵
public:
    ppForm(const Function& f, double a, double b,int nodeCount, SplineType splineType = SplineType:: Natural);

    ppForm(const Function& f, std::vector<double> nodes, SplineType splineType = SplineType:: Natural);

    ppForm(const Function& f, double a, double b,int nodeCount, bool isLinear);

    void useFunction();
    
    double LinearInterpolation(double x);

    void computeNaturalSpline();

    void computeCompleteSpline();

    void computePeriodicSpline();


    double calMiu(int i);

    double calNabuda(int i);

    double calCubicDivided(int i);

    double calQuadraticDivided(int i);

    void caltridiagonal();

    double getValue(double x);
};


