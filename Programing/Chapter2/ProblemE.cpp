#include "Newtoninterpolation.h"
#include <fstream>

int main(){
    std::vector<double> points = {0, 6, 10, 13, 17, 20, 28};
    std::vector<double> values1 = {6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7};
    std::vector<double> values2 = {6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89};

    double xStart = 0;
    double xEnd = 28;
    double step = 0.1;
    std::ofstream outfile1("outputE1.dat");
    for(double x = xStart; x <= xEnd; x += step) {
        //  if (x == 6 || x == 10 || x == 13 || x == 17 || x == 20) {
        //     continue; // 跳过这些值
        // }
        std::vector<double> xvalue, yvalue;
        xvalue.push_back(x);
        NewtonInterpolation newtonInterpolation(points, values1, x);
        yvalue.push_back(newtonInterpolation.solve());
        //std::cout << x << " " << yvalue.back() << std::endl;
        outfile1 << x << " " << yvalue.back() << std::endl;
    }
    outfile1.close();
    std::cout << "Interpolation1 finished! Generate outputE1.dat" << std::endl;

    std::ofstream outfile2("outputE2.dat");
    for(double x = xStart; x <= xEnd ; x += step) {
        //  if (x == 6 || x == 10 || x == 13 || x == 17 || x == 20) {
        //     continue; // 跳过这些值
        // }
        std::vector<double> xvalue, yvalue;
        xvalue.push_back(x);
        NewtonInterpolation newtonInterpolation(points, values2, x);
        yvalue.push_back(newtonInterpolation.solve());
        //std::cout << x << " " << yvalue.back() << std::endl;
        outfile2 << x << " " << yvalue.back() << std::endl;
    }
    outfile2.close();
    std::cout << "Interpolation2 finished! Generate outputE2.dat" << std::endl;

    NewtonInterpolation A1(points, values1, 43);
    NewtonInterpolation A2(points, values2, 43);
    std::cout << "A1(43) = " << A1.solve() << std::endl;
    std::cout << "A2(43) = " << A2.solve() << std::endl;


    return 0;

}