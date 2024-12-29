#include "include/BSpline.cpp"
#include <fstream>


int main(){
    std::ofstream outfile("bin/testM.dat");
    double xStart = -2, xEnd = 2, dx = 0.01;

    class f1 :public Function{
    public:
        double operator()(double x) const{
            return exp(x);
        }
    };

    class f2 :public Function{
    public:
        double operator()(double x) const{
            return sin(x);
        }
    };
    
    class f3 :public Function{
    public:
        double operator()(double x) const{
            return exp(x)*sin(x);
        }
    };

    BSpline obj(3, f3(), -5, 5, 81, SplineType::Natural);
    // BSpline obj(2, f2(), -4.5, 4.5, 10);
    for(double i = xStart+dx; i < xEnd; i += dx){
        outfile << i << " " << obj.getValue(i) << std::endl;
    }
    outfile.close();
    return 0;
}