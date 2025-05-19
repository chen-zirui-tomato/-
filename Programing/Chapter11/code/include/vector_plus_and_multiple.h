#pragma once
#include <vector>
#include <algorithm>


inline std::vector<double> operator+(
    const std::vector<double>& a,
    const std::vector<double>& b
) {

    std::vector<double> result(a.size()); 
    for (int i = 0; i < a.size(); ++i) {  
        result[i] = a[i] + b[i];
    }
    return result;
}

inline std::vector<double> operator-(
    const std::vector<double>& a,
    const std::vector<double>& b
) {
    std::vector<double> result(a.size()); 
    for (int i = 0; i < a.size(); ++i) {  
        result[i] = a[i] - b[i];
    }
    return result;
}

// 复合加法赋值运算符 (v1 += v2)
inline std::vector<double>& operator+=(
    std::vector<double>& a,
    const std::vector<double>& b
) {
    for (int i = 0; i < a.size(); ++i) {  
        a[i] += b[i];
    }
    return a;
}

// 数乘运算符 (vector * scalar)
inline std::vector<double> operator*(const std::vector<double>& vec, double scalar) {
    std::vector<double> result(vec.size());
    std::transform(vec.begin(), vec.end(), result.begin(),
                   [scalar](double x) { return x * scalar; });
    return result;
}

// 数乘运算符 (scalar * vector)
inline std::vector<double> operator*(double scalar, const std::vector<double>& vec) {
    return vec * scalar; // 复用上述实现
}

// 复合赋值运算符 (vector *= scalar)
inline std::vector<double>& operator*=(std::vector<double>& vec, double scalar) {
    std::transform(vec.begin(), vec.end(), vec.begin(),
                   [scalar](double x) { return x * scalar; });
    return vec;
}