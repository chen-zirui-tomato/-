#include "include/ppForm.cpp"
#include "include/BSpline.cpp"
#include <fstream>


int main(){
    double error = 0.0;
    std::ofstream outfile("bin/testA.dat");
    double xStart = -1, xEnd = 1, dx = 0.01;

    class f1 :public Function{
    public:
        double operator()(double x) const{
            return 1/(1+25*x*x);
        }
    };
    //画图
    
    ppForm obj(f1(), -1, 1, 11, SplineType::Complete);
    for(double i = xStart+dx; i < xEnd; i += dx){
        outfile << i << " " << obj.getValue(i) << std::endl;
        if(std::fabs(obj.getValue(i) - f1()(i)) > error) error = std::fabs(obj.getValue(i) - f1()(i));
    }

    // ppForm p(f1(), -1, 1, 21, true);
    // for(double i = xStart; i <= xEnd; i += dx){
    //     outfile << i << " " << p.LinearInterpolation(i) << std::endl;
    // }
    // BSpline obj(3, f1(), -1, 1, 41, SplineType::Complete);
    // for(double i = xStart+dx; i < xEnd; i += dx){
    //     outfile << i << " " << obj.getValue(i) << std::endl;
    // }

    std::cerr << "Errormax: " << error << std::endl;
    outfile.close();
    return 0;
}