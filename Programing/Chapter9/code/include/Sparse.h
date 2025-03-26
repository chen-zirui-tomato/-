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

    sparseVector(int);
    sparseVector(const sparseVector& other);
    const double operator()(int) const;
    sparseVector& operator=(const sparseVector& other);
    sparseVector& operator=(sparseVector&& other);
    ~sparseVector();

};

sparseVector::sparseVector(int size) : size(size) {
    data.resize(size);
}

sparseVector::sparseVector(const sparseVector& other) {
    size = other.size;
    data = other.data;
}

void sparseVector::set_value(int index, double value){
    if(value != 0)
        for(auto& p : data)
        if(p.first == index){
            p.second = value;
            return;
        }
        else
            data.emplace_back(index, value);
    else
        for(auto it = data.begin(); it != data.end();) {
            if(it -> first == index) it = data.erase(it);
            else ++it;
            return;
        }
}

const double sparseVector::operator()(int index) const{
    if(index <= 0 || index > size)
        throw std::out_of_range("Index out of range");
    for(auto& p : data) if(p.first == index)
        return p.second;
    return 0;
}

sparseVector & sparseVector::operator = (const sparseVector & rhs){
    if(this == &rhs) return *this;
    sparseVector copy(rhs);
    std::swap(*this, copy);
    return *this;
}

sparseVector & sparseVector::operator = (sparseVector && rhs){
    std::swap(data, rhs.data);
    size = rhs.size;
    return *this;
}

sparseVector::~sparseVector(){}


class sparseMatrix {
public:
    int row, col;
    //每一行用一个向量存储对应非0位置的元素及其值
    std::vector<std::vector<std::pair<int, double>>> data;

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
};

sparseMatrix::sparseMatrix(int row, int col) : row(row), col(col) {
    data.resize(row);
}

sparseMatrix::sparseMatrix(const sparseMatrix& other) {
    row = other.row;
    col = other.col;
    data = other.data;
}

sparseMatrix & sparseMatrix::operator = (const sparseMatrix & rhs){
    if(this == &rhs) return *this;
    sparseMatrix copy(rhs);
    std::swap(*this, copy);
    return *this;
}

sparseMatrix & sparseMatrix::operator = (sparseMatrix && rhs){
    std::swap(data, rhs.data);
    row = rhs.row;
    col = rhs.col;
    return *this;
}

sparseMatrix::~sparseMatrix(){}

sparseMatrix sparseMatrix::Triverse() const {
    sparseMatrix result(col, row);
    for(int i = 1; i <= row; i++) for(auto& p : data[i-1]) result.set_value(p.first, i, p.second);
    return result;
}

void sparseMatrix::tidyUp(){
    for(int i = 1; i <= row; i++)
        for(auto it = data[i-1].begin(); it!= data[i-1].end(); )
            if(it->second == 0) data[i-1].erase(it++);
            else ++it;
}

void sparseMatrix::set_value(int row, int col, double value){
    //row-1是因为矩阵的行数从1开始，而数组的下标从0开始；
    //对于列的位置则不做改变，保持数学上的一致性。之后的操作同理。
    if(value != 0)
        for(auto& p : data[row-1])
        if(p.first == col){
            p.second = value;
            return;
        }
        else{
            data[row-1].emplace_back(col, value);
            return;
        }
    else
        for(auto it = data[row-1].begin(); it!= data[row-1].end();) {
            if(it -> first == col) it = data[row-1].erase(it);
            else ++it;
            return;
        }
}

double& sparseMatrix::operator()(int i, int j) {
    if(i <= 0 || i > row || j <= 0 || j > col)
        throw std::out_of_range("Index out of range");
    for(auto& p : data[i-1]) if(p.first == j)
        return p.second;
    data[i-1].emplace_back(j, 0.0);
    return data[i-1].back().second;
        
}

const double& sparseMatrix::operator()(int i, int j) const {
    for(auto& p : data[i-1])
        if(p.first == j)
            return p.second;
    return 0;
}

sparseMatrix sparseMatrix::operator+(const sparseMatrix& other) const {    
    if(row != other.row || col != other.col)
        throw std::invalid_argument("Matrix + size not match");
    sparseMatrix result(row, col); 
    for(int i = 1; i <= row; i++) for(auto& p : data[i-1]) for(auto& q : other.data[i-1]) if(p.first == q.first)
        result.set_value(i, p.first, p.second + q.second);
    return result;
}

sparseMatrix sparseMatrix::operator-(const sparseMatrix& other) const {
    if(row != other.row || col != other.col)
        throw std::invalid_argument("Matrix - size not match");
    sparseMatrix result(row, col); 
    for(int i = 1; i <= row; i++) for(auto& p : data[i-1]) for(auto& q : other.data[i-1]) if(p.first == q.first)
        result.set_value(i, p.first, p.second - q.second);
    return result;
}

sparseMatrix sparseMatrix::operator*(const sparseMatrix& other) const {
    if(col != other.row)
        throw std::invalid_argument("Matrix * size not match");
    sparseMatrix result(row, other.col); 
    sparseMatrix otherT = other.Triverse();
    //固定行
    for(int i = 1; i <= row; i++) 
    //遍历列
    for(int j = 1; j <= other.col; j++) {
        double sum = 0;
        //对于i行每个非零元，找other的j列中对应非0元计算乘积
        for(auto& p : data[i-1]) 
            for(auto& q : otherT.data[j-1]) if(p.first == q.first){
                sum += p.second * q.second;
                break;
            }
        result.set_value(i, j, sum);
    }
    return result;
}

sparseMatrix sparseMatrix::operator*(double scalar) const {
    sparseMatrix result(row, col);
    for(int i = 1; i <= row; i++) for(auto& p : data[i-1])
        result.set_value(i, p.first, p.second * scalar);
    return result;
}

sparseMatrix sparseMatrix::operator/(double scalar) const {
    if(scalar == 0)
        throw std::invalid_argument("Divide by zero");
    sparseMatrix result(row, col);
    for(int i = 1; i <= row; i++) for(auto& p : data[i-1])
        result.set_value(i, p.first, p.second / scalar);
    return result;
}

void sparseMatrix::print() const {
    for(int i = 1; i <= row; i++){
        for(auto& p : data[i-1])
            std::cout << p.second << " ";
        std::cout << std::endl;
    }
}

sparseVector sparseMatrix::operator*(const sparseVector& other) const {
    if(col != other.size)
        throw std::invalid_argument("Matrix * size not match");
    sparseVector result(row);
    for(int i = 1; i <= row; i++){
        double sum = 0;
        for(auto& p : data[i-1]) 
            sum += p.second * other(p.first);
        if(sum != 0) result.set_value(i, sum);
    }
    return result;
}



#endif