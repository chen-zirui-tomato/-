#include "Newtoninterpolation.h"
#include <fstream>

class f :public Function{
public:
    double operator()(double x) const {
        return 1/(1 + x * x);
    }
};

int main(){
    
    std::cout << "use Newton interpolation :" << std::endl;
    /**
     * @param func is the function to be solved
     */
    f func ;

    std::vector<double> fx;

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
    for(int i = 0 ; i <= n ; i++){
        points.push_back(-5.0 + 10.0*i/n);
    }

    for(int i = 0 ; i <= n ; i++){
        fx.push_back(func(points[i]));
    }

    double xStart = -5;
    double xEnd = 5;
    double step = 0.1;

    //把计算结果写入文件
    std::ofstream outfile("outputB.dat");
    /**
     * @param xvalue is the vector of x-axis values
     * @param yvalue is the vector of y-axis values
     */
    std::vector<double> xvalue, yvalue;
    for(double x = xStart; x <= xEnd; x += step) {
        xvalue.push_back(x);
        NewtonInterpolation newtonInterpolation(points, fx, x);
        yvalue.push_back(newtonInterpolation.solve());
        //std::cout << x << " " << yvalue.back() << std::endl;
        outfile << x << " " << yvalue.back() << std::endl;
    }
    outfile.close();
    std::cout << "Interpolation finished! Generate outputB.dat" << std::endl;
    return 0;
}