#include "Beziercurve.h"
#include <cmath>
#include <fstream>

class F1 : public Function {
public:
    double operator()(double x)const {
        return 2.0/3.0*(sqrt(3-x*x) + sqrt(fabs(x)));
    }
};

class F2 : public Function {
public:
    double operator()(double x)const {
        return 2.0/3.0*(-sqrt(3-x*x) + sqrt(fabs(x)));
    }
};


int main() {

    double xStart = -sqrt(2.9);
    double xEnd = sqrt(2.9);
    int m = 40;
    double steps = 0.01;
    std::vector<double> points;
    std::vector<double> f1values;
    std::vector<double> f1values_;
    std::vector<double> f2values;
    std::vector<double> f2values_;
    F1 f1;
    F2 f2;

    //跟书上一样的m+1个点。
    for(int i = 0; i <= m; ++i){
        points.push_back((xStart + i*(xEnd-xStart)/m));
        f1values.push_back(f1(points[i]));
        f1values_.push_back(f1.derivative(points[i]));
        f2values.push_back(f2(points[i]));
        f2values_.push_back(f2.derivative(points[i]));
    // std::cout << points[i]<<" " <<values[i]<< " " <<values_[i]<< std::endl;
    }   

    // BezierCurve bezier(points, values, values_, m+1);
    // bezier.segment();
    // std::cout << "Bezier curve constructed." << std::endl;
    // double A = bezier.solveP(xStart);
    // std::cout << "P(0.5) = " << A << std::endl;


    BezierCurve bezier1(points, f1values, f1values_, m+1);
    bezier1.segment();
    BezierCurve bezier2(points, f2values, f2values_, m+1);
    bezier2.segment();
    std::cout << "Bezier curve constructed." << std::endl;

        
    std::ofstream outfile("outputF1.dat");
    std::vector<double> xvalue,y1value, y2value;
    for(double x = xStart; x <= xEnd; x += steps){
        xvalue.push_back(x);
        y1value.push_back(bezier1.solveP(x));
        outfile << x << " " << y1value.back() << std::endl;
    }
    outfile.close();

    points.clear();
    std::ofstream outfile2("outputF2.dat");
    for(double x = xStart; x <= xEnd; x += steps){
        y2value.push_back(bezier2.solveP(x));
        outfile2 << x << " " << y2value.back() << std::endl;
    }
    outfile2.close();
    std::cout << "Interpolation finished! Generated outputF.dat." << std::endl;
    return 0;
}