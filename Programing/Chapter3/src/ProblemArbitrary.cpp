#include "include/BSpline.cpp"
#include <fstream>


int main(){
    std::ofstream outfile("bin/testAR.dat");
    double xStart = -5, xEnd = 5, dx = 0.01;

    class f :public Function{
    public:
        double operator()(double x) const{
            return 1/(1+x*x);
        }
    };
    
    BSpline obj(7, f(), -5, 5, 81);
    for(double i = xStart+dx; i < xEnd; i += dx)
        outfile << i << " " << obj.getValue(i) << std::endl;
        // std::cout << i << " " << obj.getValue(i) << std::endl;
    outfile.close();
    return 0;
}