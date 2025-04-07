#ifndef Sparse_H
#define Sparse_H

#include <vector>
#include <Eigen/Dense>
//=================================================
//这里所有调用的i，j都代指数学上的下标，与数组下标作区分！！！
//=================================================

class sparseVector {
public:
    int Size;
    std::vector<std::pair<int, double>> data;
    sparseVector();
    
    // Convert to Eigen dense vector
    Eigen::VectorXd toDense() const {
        Eigen::VectorXd dense(Size);
        dense.setZero();
        for (const auto& elem : data) {
            dense(elem.first - 1) = elem.second;  // Convert to 0-based index
        }
        return dense;
    }
    //生成长度为int的零向量
    sparseVector(int);
    sparseVector(const sparseVector& other);

    void set_value(int index, double value);
    void tidyUp();
    void print() const;
    [[nodiscard]] double norm() const;
    [[nodiscard]] int size() const;
    void resize(int newSize);
    const double operator()(int) const;
    sparseVector operator+(const sparseVector& other) const;
    sparseVector operator-(const sparseVector& other) const;
    sparseVector operator*(double scalar) const;
    sparseVector operator/(double scalar) const;
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
    
    // Convert to Eigen dense matrix
    [[nodiscard]] Eigen::MatrixXd toDense() const {
        Eigen::MatrixXd dense(row, col);
        dense.setZero();
        for (int i = 0; i < row; ++i) {
            for (const auto& elem : data[i]) {
                dense(i, elem.first - 1) = elem.second;  // Convert to 0-based index
            }
        }
        return dense;
    }
    //生成row行的零矩阵
    sparseMatrix(int row, int col);
    sparseMatrix(const sparseMatrix& other);
    sparseMatrix operator=(const sparseMatrix& other);
    sparseMatrix operator=(sparseMatrix&& other);
    ~sparseMatrix();

    // //实现matrix(i,j)=value!=0
    // double& operator()(int i, int j);
    //读取matrix(i,j)
    double operator()(int i, int j) const;

    [[nodiscard]] sparseMatrix Triverse() const;
    //获得D
    [[nodiscard]] sparseMatrix Diagonal() const;
    //获得D^{-1}
    [[nodiscard]] sparseMatrix InverseDiagonal() const;
    //消除0元
    void tidyUp();
    [[nodiscard]] int size() const;
    void set_value(int row, int col, double value);
    void print() const;

    sparseMatrix operator+(const sparseMatrix& other) const;
    sparseMatrix operator-(const sparseMatrix& other) const;
    sparseMatrix operator*(const sparseMatrix& other) const;
    sparseMatrix operator*(double scalar) const;
    sparseMatrix operator/(double scalar) const;
  
    sparseVector operator*(const sparseVector& other) const;
};


#endif