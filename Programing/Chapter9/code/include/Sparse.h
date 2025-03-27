#ifndef Sparse_H
#define Sparse_H

#include <iostream>
#include <vector>
//=================================================
//这里所有调用的i，j都代指数学上的下标，与数组下标作区分！！！
//=================================================

class sparseVector {
public:
    int size;
    std::vector<std::pair<int, double>> data;

    void set_value(int index, double value);

    sparseVector();
    sparseVector(int);
    sparseVector(const sparseVector& other);
    const double operator()(int) const;
    sparseVector& operator=(const sparseVector& other);
    sparseVector& operator=(sparseVector&& other);
    ~sparseVector();

};

class sparseMatrix {
public:
    int row, col;
    //每一行用一个向量存储对应非0位置的元素及其值
    std::vector<std::vector<std::pair<int, double>>> data;

    sparseMatrix();
    sparseMatrix(int row, int col);
    sparseMatrix(const sparseMatrix& other);
    sparseMatrix& operator=(const sparseMatrix& other);
    sparseMatrix& operator=(sparseMatrix&& other);
    ~sparseMatrix();

    //实现matrix(i,j)=value!=0
    double& operator()(int i, int j);
    //读取matrix(i,j)
    const double& operator()(int i, int j) const;

    sparseMatrix Triverse() const;
    //消除0元
    void tidyUp();
    void set_value(int row, int col, double value);
    void print() const;

    sparseMatrix operator+(const sparseMatrix& other) const;
    sparseMatrix operator-(const sparseMatrix& other) const;
    sparseMatrix operator*(const sparseMatrix& other) const;
    sparseMatrix operator*(double scalar) const;
    sparseMatrix operator/(double scalar) const;  

    sparseVector operator*(const sparseVector& other) const;
    
    // 提取对角矩阵
    sparseMatrix extractDiagonal() const;
    
    // 加权Jacobi迭代求解 Ax = b
    // 参数: b - 右端项, x0 - 初始猜测, omega - 松弛因子(0 < omega < 1)
    // max_iter - 最大迭代次数, tol - 收敛容差
    // 返回: 解向量和实际迭代次数
    std::pair<sparseVector, int> weightedJacobi(
        const sparseVector& b, 
        const sparseVector& x0, 
        double omega = 0.7, 
        int max_iter = 1000, 
        double tol = 1e-6) const;
};


#endif