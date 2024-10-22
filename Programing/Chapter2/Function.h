#ifndef FUNCTION
#define FUNCTION
#include <iostream>
#include <vector>

class Function{
public:
    virtual double operator()(double x) const = 0;
    virtual double derivative(double x) const {
        double h = 1e-6;
        return (operator()(x+h) - operator()(x-h))/(2*h);
    }
    virtual ~Function(){}
};

class Interpolation{
protected:
    // const Function& F;
    /**
     * @param ponits vertor of dictinct points
     */
    const std::vector<double> points;

    /**
     * @param values verctor of values at points
     */
    const std::vector<double> values;

    /**
     * @param n the size of points 
     */
    const int n = points.size();

    /**
     * @param x a point in domain
     */
    const double x;

    /**
     * @brief differ matrix
     */
    std::vector<std::vector<double>> differs;

public:
    Interpolation(std::vector<double> points, std::vector<double> values, double x):
                    points(points), values(values), x(x){};
    /**
     * @brief 求解函数
     */
    double solve();

    /**
     * @brief 生成差商函数
     */
    virtual void generatediffers() = 0;

    ~Interpolation(){}
};

double Interpolation::solve(){
    // std::cerr << "solving..." << std::endl;
    // for(int i = 0; i < n; i++){
    //     std::cout << "i = " << i << "  points[i] = " << points[i] << "  f[i] = " << f[i] << std::endl;
    // }
    // for(int i = 0; i < n; i++){
    //     for(int j = 0; j < n; j++){
    //         std::cout << "differs[" << i << "][" << j << "] = " << differs[i][j] << std::endl;
    //     }}
    generatediffers();
    double fx = 0.0;
    double h ;
    double x = this->x;
    for(int i = n ; i >= 2; --i){
        h = x - points[i - 2];
        fx = h * (fx + differs[i - 1][i - 1]);
        // std::cout <<i<< "  h = " << h << " fx = " << fx << "differ(i-1) = " << differ(i-1) << std::endl;
    }
    fx += differs[0][0];
    // std::cout << "final fx = " << fx << std::endl;
    // std::cout << std::endl;
    return fx ;
};
#endif