#include "include/BSpline.cpp"
#include <fstream>


int main(){
    std::ofstream outfile("bin/testB.dat");
    double xStart = -5, xEnd = 5, dx = 0.01;

    class f2 :public Function{
    public:
        double operator()(double x) const{
            return 1/(1+x*x);
        }
    };
    
    BSpline obj(3, f2(), -5, 5, 11, SplineType::Complete);
    // BSpline obj(2, f2(), -4.5, 4.5, 10);
    for(double i = xStart+dx; i < xEnd; i += dx){
        outfile << i << " " << obj.getValue(i) << std::endl;
    }
    double error = std::fabs(obj.getValue(-3.5) - f2()(-3.5));
    std::cout << "Error: " << error << std::endl;
    outfile.close();
    return 0;
}