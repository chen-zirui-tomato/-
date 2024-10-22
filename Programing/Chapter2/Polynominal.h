#include <vector>
#include <iostream>
#include <cmath>

class Polynominal {
protected:
/**
 * @brief 存放多项式系数
 */
std::vector<double> coeff;

/**
 * @brief 多项式次数
 */
int n;

public:
/**
 * @brief 默认构造函数
 */
Polynominal()=default;

/**
 * @brief 常数
 */
Polynominal(double a): coeff{a},n(0){};

/**
 * @brief 单项式
 */
Polynominal(double a, double b): coeff{a,b},n(1){};

/**
 * @brief 一般多项式
 */
Polynominal(std::vector<double> coeff, int n)
            :coeff(coeff),n(n){};
/**
 * @brief 析构函数
 */   
~Polynominal(){};

Polynominal operator+(const Polynominal& other);
Polynominal operator*(const Polynominal& other);
double calculate(double x) const;
void print() const;
};

Polynominal Polynominal::operator+(const Polynominal& other) {
    int maxDegree = std::max(n, other.n);
    std::vector<double> newCoeff(maxDegree + 1, 0.0);
    
    for (int i = 0; i <= n; ++i) {
        newCoeff[i] += coeff[i];
    }
    for (int i = 0; i <= other.n; ++i) {
        newCoeff[i] += other.coeff[i];
    }

    return Polynominal(newCoeff, maxDegree);
}

Polynominal Polynominal::operator*(const Polynominal& other) {
    if (n == 0 || other.n == 0) {
        return Polynominal(0.0, 0);
    }
    int newDegree = n + other.n;
    std::vector<double> newCoeff(newDegree + 1, 0.0);

    for (int i = 0; i <= n; ++i) {  
        for (int j = 0; j <= other.n; ++j) {
            newCoeff[i + j] += coeff[i] * other.coeff[j];
        }
    }

    return Polynominal(newCoeff, newDegree);
}

void Polynominal::print() const {
    for (int i = 0; i <= n; ++i) {
        std::cout << coeff[i] << "x^" << i;
        if (i < n) std::cout << " + "; // 输出时加上加号
    }
    std::cout << std::endl;
}

double Polynominal::calculate(double x) const{
    double result = 0.0;
    for(int i = n; i >= 0 ; --i){
        result += coeff[i] * pow(x, i);
    }
    return result;
}