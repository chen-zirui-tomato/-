#include "Newtoninterpolation.h"
#include <cmath>
#include <fstream>

class f :public Function{
public:
    double operator()(double x) const{
        return 1.0 / (1.0 + 25.0 * x * x);
    }
};

int main(){
    std::cout << "use chebyshev points" << std::endl;
    /**
     * @param func is the function to be solved
     */
    f func ;
    /**
     * @param n is the size of points
     */
    int n ;
    std::cout << "please enter the value of n : " << std::endl;
    std::cin >> n ;

    /**
     * @param points is the vector of point
     */
    std::vector<double> points;

    /**
     * @param
     */
    std::vector<double> fx;

    double pi = acos(-1);
    // std::cout << pi;
    for(int k = 1 ; k <= n ; k++){
        points.push_back(cos((2.0*k - 1.0)/(2.0*n)*pi));
    }

    for(int k = 1; k <= n; k++){
        fx.push_back(func(points[k-1]));
    }
    
    double xStart = -1;
    double xEnd = 1;
    double step = 0.05;

    //把计算结果写入文件
    std::ofstream outfile("outputC.dat");
    /**
     * @param xvalue is the vector of x-axis values
     * @param yvalue is the vector of y-axis values
     */
    std::vector<double> xvalue, yvalue;
    for(double x = xStart; x <= xEnd; x += step) {
        xvalue.push_back(x);
        NewtonInterpolation newtonInterpolation(points, fx, x);
        yvalue.push_back(newtonInterpolation.solve());
        // std::cout << x << " " << yvalue.back() << std::endl;
        outfile << x << " " << yvalue.back() << std::endl;
    }
    outfile.close();
    std::cout << "Interpolation finished! Generated outputC.dat." << std::endl;
    return 0;
}
