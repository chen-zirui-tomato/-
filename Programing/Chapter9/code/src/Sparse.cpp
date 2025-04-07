#include "../include/Sparse.h"
#include <cmath>
#include <iostream>
#include <map>

sparseVector::sparseVector():Size(-1), data(){}

sparseVector::sparseVector(int size) : Size(size) {}

sparseVector::sparseVector(const sparseVector& other) {
    Size = other.Size;
    data = other.data;
}

int sparseVector::size() const{
    return Size;
}

void sparseVector::resize(int size){
    Size = size;
}

void sparseVector::tidyUp(){
    for(auto it = data.begin(); it!= data.end(); )
        if(it->second == 0) {
            // std::cerr<< "Debug - sparseVector::tidyUp called with index: " << it->first << std::endl;
            it = data.erase(it);
        }
        else ++it;
}

void sparseVector::print() const{
    for(auto& p : data)
        std::cout << p.first << " " << p.second << std::endl;
}

double sparseVector::norm() const{
    double sum = 0;
    for(auto& p : data)
        sum += p.second * p.second;
    return sqrt(sum);
}

void sparseVector::set_value(int index, double value){
    if(index <= 0 || index > Size){
        std::cerr << "Debug - sparseVector::set_value called with index: " << index << ", value: " << value << std::endl;
        throw std::out_of_range("Vector set value index out of range");
    }
    for(auto& p : data) {
        if(p.first == index){
            p.second = value;
            return;
        }
    }
    data.emplace_back(index, value);
}

const double sparseVector::operator()(int index) const{
    if(index <= 0 || index > Size){
        std::cerr << "this size = " << Size << " index = " << index << std::endl;
        throw std::out_of_range("Vector() index out of range");
    }
    for(auto& p : data) if(p.first == index)
        return p.second;
    return 0;
}

sparseVector sparseVector::operator+(const sparseVector& other) const {
    if(Size != other.Size)
        throw std::invalid_argument("Vector + size not match");
    
    std::map<int, double> temp_result;

    for(const auto& elem : data) {
        temp_result[elem.first] += elem.second;
    }
    for(const auto& elem : other.data) {
        temp_result[elem.first] += elem.second;
    }
    
    sparseVector result(Size);
    for(const auto& pair : temp_result) {
        if(pair.second != 0) {
            result.data.emplace_back(pair.first, pair.second);
        }
    }
    
    return result;
}

sparseVector sparseVector::operator-(const sparseVector& other) const {
    if(Size != other.Size)
        throw std::invalid_argument("Vector - size not match");
    
    std::map<int, double> temp_result;

    for(const auto& elem : data) {
        temp_result[elem.first] += elem.second;
    }
    for(const auto& elem : other.data) {
        temp_result[elem.first] -= elem.second;
    }
    
    sparseVector result(Size);
    for(const auto& pair : temp_result) {
        if(pair.second != 0) {
            result.data.emplace_back(pair.first, pair.second);
        }
    }
    
    return result;
}

sparseVector sparseVector::operator*(double scalar) const {
    sparseVector result(Size);
    for(auto& p : data)
        result.set_value(p.first, p.second * scalar);
    result.tidyUp();
    return result;
}

sparseVector sparseVector::operator/(double scalar) const {
    if(scalar == 0)
        throw std::invalid_argument("Vector / 0");
    sparseVector result(Size);
    for(auto& p : data)
        result.set_value(p.first, p.second / scalar);
    result.tidyUp();
    return result;
}

sparseVector sparseVector::operator = (const sparseVector & rhs){
    if(this == &rhs) return *this;
    sparseVector copy(rhs);
    std::swap(*this, copy);
    return *this;
}

sparseVector sparseVector::operator = (sparseVector && rhs){
    std::swap(data, rhs.data);
    Size = rhs.Size;
    return *this;
}

sparseVector::~sparseVector(){

}

//==========================================================================

sparseMatrix::sparseMatrix():row(-1), col(-1), data(){}

sparseMatrix::sparseMatrix(int row, int col) : row(row), col(col) {
    data.resize(row);
}

sparseMatrix::sparseMatrix(const sparseMatrix& other) {
    row = other.row;
    col = other.col;
    data = other.data;
}

sparseMatrix sparseMatrix::operator = (const sparseMatrix & rhs){
    if(this == &rhs) return *this;
    sparseMatrix copy(rhs);
    std::swap(*this, copy);
    return *this;
}

sparseMatrix sparseMatrix::operator = (sparseMatrix && rhs){
    std::swap(data, rhs.data);
    row = rhs.row;
    col = rhs.col;
    return *this;
}

sparseMatrix::~sparseMatrix(){}

sparseMatrix sparseMatrix::Triverse() const {
    sparseMatrix result(col, row);
    for(int i = 1; i <= row; i++) 
    for(auto& p : data[i-1]) 
        result.set_value(p.first, i, p.second);
    return result;
}

sparseMatrix sparseMatrix::Diagonal() const {
    sparseMatrix result(row, col);
    for(int i = 1; i <= row; i++) {
        result.set_value(i, i, data[i-1][i].second);
    }
    return result;
}

sparseMatrix sparseMatrix::InverseDiagonal() const {
    sparseMatrix result(row, col);
    for(int i = 1; i <= row; i++) {
        if(data[i-1][i].second == 0) 
            throw std::invalid_argument("Matrix is singular");
        result.set_value(i, i, 1.0/data[i-1][i].second);
    }
    return result;
}

int sparseMatrix::size() const {
    return row;
}

void sparseMatrix::tidyUp(){
    for(int i = 1; i <= row; i++)
        for(auto it = data[i-1].begin(); it!= data[i-1].end(); )
            if(it->second == 0) data[i-1].erase(it++);
            else ++it;
}

void sparseMatrix::set_value(int i, int j, double value){
    //row-1是因为矩阵的行数从1开始，而数组的下标从0开始；
    //对于列的位置则不做改变，保持数学上的一致性。之后的操作同理。
    if(i <= 0 || i > row || j <= 0 || j > col)
        throw std::out_of_range("Matrix set value index out of range");
    if(value != 0)
    for(auto& p : data[i-1]){
        if(p.first == j){
            p.second = value;
            return;
        }
    }
    data[i-1].emplace_back(j, value);
}

// double& sparseMatrix::operator()(int i, int j) {
//     if(i <= 0 || i > row || j <= 0 || j > col)
//         throw std::out_of_range("Index out of range");
//     for(auto& p : data[i-1]) if(p.first == j)
//         return p.second;
//     data[i-1].emplace_back(j, 0.0);
//     return data[i-1].back().second;
        
// }

double sparseMatrix::operator()(int i, int j) const {
    if(i <= 0 || i > row || j <= 0 || j > col)
        throw std::out_of_range("Matrix() Index out of range");
    for(auto& p : data[i-1])
        if(p.first == j)
            return p.second;
    return 0;
}

sparseMatrix sparseMatrix::operator+(const sparseMatrix& other) const {    
    if(row != other.row || col != other.col)
        throw std::invalid_argument("Matrix + size not match");
    
    sparseMatrix result(row, col);
    for(int i = 1; i <= row; i++) {
        std::map<int, double> row_result;
        for(auto& elem : data[i-1]) {
            row_result[elem.first] += elem.second;
        }
        for(auto& elem : other.data[i-1]) {
            row_result[elem.first] += elem.second;
        }
        for(auto& pair : row_result) {
            if(pair.second != 0) {
                result.data[i-1].emplace_back(pair.first, pair.second);
            }
        }
    }
    return result;
}

sparseMatrix sparseMatrix::operator-(const sparseMatrix& other) const {
    if(row != other.row || col != other.col){
        std::cerr << "this size: " << row << " " << col << " other size: " << other.row << " " << other.col <<std::endl;
        throw std::invalid_argument("Matrix - size not match");
    }
    sparseMatrix result(row, col);
    for(int i = 1; i <= row; i++) {
        std::map<int, double> row_result;
        for(auto& elem : data[i-1]) {
            row_result[elem.first] += elem.second;
        }
        for(auto& elem : other.data[i-1]) {
            row_result[elem.first] -= elem.second;
        }
        for(auto& pair : row_result) {
            if(pair.second != 0) {
                result.data[i-1].emplace_back(pair.first, pair.second);
            }
        }
    }
    return result;
}

// sparseMatrix sparseMatrix::operator*(const sparseMatrix& other) const {
//     if(col != other.row)
//         throw std::invalid_argument("Matrix * size not match");
//     sparseMatrix result(row, other.col); 
//     sparseMatrix otherT = other.Triverse();
//     //固定行
//     for(int i = 1; i <= row; i++) 
//     //遍历列
//     for(int j = 1; j <= other.col; j++) {
//         double sum = 0;
//         //对于i行每个非零元，找other的j列中对应非0元计算乘积
//         for(auto& p : data[i-1]) 
//             for(auto& q : otherT.data[j-1]) if(p.first == q.first){
//                 sum += p.second * q.second;
//                 break;
//             }
//         result.set_value(i, j, sum);
//     }
//     result.tidyUp();
//     return result;
// }

sparseMatrix sparseMatrix::operator*(double scalar) const {
    sparseMatrix result(row, col);
    for(int i = 1; i <= row; i++) 
    for(auto& p : data[i-1])
        result.set_value(i, p.first, p.second * scalar);
    result.tidyUp();
    return result;
}

sparseMatrix sparseMatrix::operator/(double scalar) const {
    if(scalar == 0)
        throw std::invalid_argument("Divide by zero");
    sparseMatrix result(row, col);
    for(int i = 1; i <= row; i++) for(auto& p : data[i-1])
        result.set_value(i, p.first, p.second / scalar);
    result.tidyUp();
    return result;
}

void sparseMatrix::print() const {
    for(int i = 1; i <= row; i++){
        for(auto& p : data[i-1])
            std::cout << i << " " << p.first << " " << p.second << " ";
        std::cout << std::endl;
    }
}

sparseVector sparseMatrix::operator*(const sparseVector& other) const {
    if(col != other.Size){
        std::cerr << "this size: " << col << " other size: " << other.Size <<std::endl;
        throw std::invalid_argument("Matrix * vector size not match");
    }
    
    std::map<int, double> vecMap;
    for(const auto& p : other.data) {
        vecMap[p.first] = p.second;
    }
    
    sparseVector result(row);
    for(int i = 1; i <= row; i++) {
        double sum = 0;
        for(const auto& p : data[i-1]) {
            auto it = vecMap.find(p.first);
            if(it != vecMap.end()) {
                sum += p.second * it->second;
            }
        }
        if(sum != 0) {
            result.set_value(i, sum);
        }
    }
    return result;
}


