#include <stdexcept>
#include <vector>
#include <algorithm>


inline std::vector<double> operator+(
    const std::vector<double>& a,
    const std::vector<double>& b
) {
    // 强制维度校验（六维）
    if (a.size() != 6 || b.size() != 6) {
        throw std::invalid_argument("向量必须为六维 + ");
    }

    std::vector<double> result(6); 
    for (int i = 0; i < 6; ++i) {  
        result[i] = a[i] + b[i];
    }
    return result;
}

inline std::vector<double> operator-(
    const std::vector<double>& a,
    const std::vector<double>& b
) {
    // 强制维度校验（六维）
    if (a.size() != 6 || b.size() != 6) {
        throw std::invalid_argument("向量必须为六维 - ");
    }

    std::vector<double> result(6); 
    for (int i = 0; i < 6; ++i) {  
        result[i] = a[i] - b[i];
    }
    return result;
}

// 复合加法赋值运算符 (v1 += v2)
inline std::vector<double>& operator+=(
    std::vector<double>& a,
    const std::vector<double>& b
) {
    if (a.size() != 6 || b.size() != 6) {
        throw std::invalid_argument("向量必须为六维 +=");
    }

    a[0] += b[0];
    a[1] += b[1];
    a[2] += b[2];
    a[3] += b[3];
    a[4] += b[4];
    a[5] += b[5];
    
    return a;
}

// 数乘运算符 (vector * scalar)
inline std::vector<double> operator*(const std::vector<double>& vec, double scalar) {
    if(vec.size() != 6){
        throw std::invalid_argument("向量必须为六维 * ");
    }
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
    if(vec.size() != 6){
        throw std::invalid_argument("向量必须为六维 *= ");
    }
    std::transform(vec.begin(), vec.end(), vec.begin(),
                   [scalar](double x) { return x * scalar; });
    return vec;
}