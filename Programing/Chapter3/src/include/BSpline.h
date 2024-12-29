#include <iostream>
#include "Function.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>

class BSpline{
private:
    int k; //degree of bspline
    const Function& f;
    int nodeCount;
    int n; //n是插值点个数减一
    std::vector<double> nodes; //插值点
    std::vector<double> values;
    std::vector<double> P; //控制点
    std::vector<double> t; //节点
    std::vector<std::vector<double>> A; //系数矩阵
    SplineType splineType; //插值类型
    
public:
    BSpline(int k,const Function& f, double a, double b, int nodeCount, SplineType splineType = SplineType:: Natural);

    BSpline(int k, const Function& f, std::vector<double> nodes, SplineType splineType = SplineType:: Natural);

    double constructBSpline(int i,int k, double x);

    void constructMAtrix();

    void addrange();

    void addNatrualCondition();

    void addCompleteCondition();

    void addPeriodicCondition();

    void addEndCondition();
    
    void addArbitraryCondition();

    void solve();

    double getValue(double x);

    // void QRsolver(){ //QR分解
    //     std::vector<std::vector<double>> Q; //存储Q矩阵
    //     Q.resize(n+k);
    //     for(int i = 0; i <= n+k-1; i++){
    //         Q[i].resize(n+k);
    //         for(int j = 0; j <= n+k-1; j++){
    //             if(i == j) Q[i][j] = 1.0;
    //             else Q[i][j] = 0;
    //         }
    //     }

    //     //构造Q和R矩阵
    //     for(int l = 0; l < n+k-1; l++){
    //         std::vector<double> a; //a是A的第k列
    //         for(int i = l; i <= n+k-1; i++)
    //             a.push_back(A[i][l]);

    //         std::vector<double> u;
    //         u.resize(n+k-1,0.0);
    //         for(int i = 0; i <= n+k-1; i++)
    //             if(i == 0) u[i] = (a[i] - getTwoNorm(a));
    //             else u[i] = a[i];
    //         double uTemp = getTwoNorm(u);
    //         for(int i = 0; i <= n+k-1; i++)
    //             u[i] /= uTemp;

    //         std::vector<std::vector<double>> uuT = vectorMultiple(u, u);  //u乘uT 
    //         std::vector<std::vector<double>> H; //Householder矩阵
    //         H.resize(n+k);
    //         for(int i = 0; i <= n+k-1; i++){
    //             H[i].resize(n+k);
    //             for(int j = 0; j <= n+k-1; j++)
    //                 if(i == j) H[i][j] = 1.0;
    //                 else H[i][j] = 0.0;
    //         }
    //         for(int i = l; i <= n+k-1; i++)
    //             for(int j = l; j <= n+k-1; j++)
    //                 if(i == j) H[i][j] = 1.0 - 2*uuT[i][j];
    //                 else H[i][j] = -2*uuT[i][j];
            
    //         A = matrixMultiple(H, A); //更新A矩阵
    //         Q = matrixMultiple(H, Q); //更新Q矩阵
    //     }
    //     for(int i = 0; i <= n+k-1; i++) {
    //         std::cerr << std::endl;
    //         for(int j = 0; j <= n+k-1; j++) std::cerr << A[i][j] << " ";
    //     }
    //     std::cerr << std::endl;

    //     //两边同时乘以Q的转置
    //     Q = getTranspose(Q);
    //     values = matrixMultipleVector(Q, values); //更新values
    // }

    // double getTwoNorm(std::vector<double> v){ //计算二范数
    //     double sum = 0.0;
    //     for(size_t i = 0; i < v.size(); i++)
    //         sum += v[i]*v[i];
    //     return sqrt(sum);
    // }

    // std::vector<std::vector<double>> vectorMultiple(std::vector<double> a, std::vector<double> b){ //向量乘法
    //     if(a.size() != b.size()) std::cout << "Error" << std::endl;
    //     std::vector<std::vector<double>> result;
    //     result.resize(a.size());
    //     for(size_t i = 0; i < a.size(); i++){
    //         result[i].resize(b.size());
    //         for(size_t j = 0; j < b.size(); j++)
    //             result[i][j] = a[i]*b[j];
    //     }
    //     return result;
    // }

    // std::vector<std::vector<double>> matrixMultiple
    //     (std::vector<std::vector<double>> A, std::vector<std::vector<double>> B){ //矩阵乘法
    //     if(A[0].size() != B.size()) std::cout << "Error" << std::endl;
    //     std::vector<std::vector<double>> result;
    //     result.resize(A.size());
    //     for(size_t i = 0; i < A.size(); i++){
    //         result[i].resize(B[0].size());
    //         for(size_t j = 0; j < B[0].size(); j++)
    //             for(size_t k = 0; k < A[0].size(); k++)
    //                 result[i][j] += A[i][k]*B[k][j];
    //     }
    //     return result;
    // }

    // std::vector<double> matrixMultipleVector(std::vector<std::vector<double>> A, std::vector<double> b){ //矩阵与向量乘法
    //     if(A[0].size() != b.size()) std::cout << "Error" << std::endl;
    //     std::vector<double> result;
    //     result.resize(A.size());
    //     for(size_t i = 0; i < A.size(); i++)
    //         for(size_t j = 0; j < b.size(); j++)
    //             result[i] += A[i][j]*b[j];
    //     return result;
    // }

    // std::vector<std::vector<double>> getTranspose(std::vector<std::vector<double>> A){ //矩阵转置
    //     std::vector<std::vector<double>> result;
    //     result.resize(A[0].size());
    //     for(size_t i = 0; i < A[0].size(); i++){
    //         result[i].resize(A.size());
    //         for(size_t j = 0; j < A.size(); j++)
    //             result[i][j] = A[j][i];
    //     }
    //     return result;
    // }

    // std::vector<double>  solveUpperTriangular(){ //解上三角矩阵Ax=P
    //     std::vector<double> P(n+k-1,0.0);
    //     for(int i = n+k-2; i >= 0; i--){
    //         double sum = 0.0;
    //         for(int j = i+1; j <= n+k-1; j++)
    //             sum += A[i][j]*P[j];
    //         P[i] = (values[i] - sum)/A[i][i];
    //     }
    //     return P;
    // }
};