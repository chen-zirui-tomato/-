#ifndef Sparse_H
#define Sparse_H

#include <vector>
//=================================================
//这里所有调用的i，j都代指数学上的下标，与数组下标作区分！！！
//=================================================

class sparseVector {
public:
    int Size;
    std::vector<std::pair<int, double>> data;
    sparseVector();
    //生成长度为int的零向量
    sparseVector(int);
    sparseVector(const sparseVector& other);

    void set_value(int index, double value);
    void tidyUp();
    [[nodiscard]] int size() const;
    void resize(int newSize);
    const double operator()(int) const;
    sparseVector operator+(const sparseVector& other) const;
    sparseVector operator-(const sparseVector& other) const;
    sparseVector operator*(double scalar);
    sparseVector operator/(double scalar);
    sparseVector operator=(const sparseVector& other);
    sparseVector operator=(sparseVector&& other);
    ~sparseVector();
};

class sparseMatrix {
public:
    int row, col;
    //每一行用一个向量存储对应非0位置的元素及其值
    std::vector<std::vector<std::pair<int, double>>> data;

    sparseMatrix();
    //生成row行的零矩阵
    sparseMatrix(int row, int col);
    sparseMatrix(const sparseMatrix& other);
    sparseMatrix operator=(const sparseMatrix& other);
    sparseMatrix operator=(sparseMatrix&& other);
    ~sparseMatrix();

    //实现matrix(i,j)=value!=0
    double& operator()(int i, int j);
    //读取matrix(i,j)
    double operator()(int i, int j) const;

    [[nodiscard]] sparseMatrix Triverse() const;
    //消除0元
    void tidyUp();
    [[nodiscard]] int size() const;
    void set_value(int row, int col, double value);
    void print() const;

    sparseMatrix operator+(const sparseMatrix& other) const;
    sparseMatrix operator-(const sparseMatrix& other) const;
    sparseMatrix operator*(const sparseMatrix& other) const;
    sparseMatrix operator*(double scalar);
    sparseMatrix operator/(double scalar) ;
  
    sparseVector operator*(const sparseVector& other);
};


#endif