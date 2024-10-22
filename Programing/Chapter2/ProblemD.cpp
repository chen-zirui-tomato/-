#include "Hermiteinterpolation.h"
#include <cmath>
#include <fstream>

int main(){

    std::vector<double> points = {0, 0, 3, 3, 5, 5, 8, 8, 13, 13};
    std::vector<double> f = {0, 0, 225, 225, 383, 383, 623, 623, 993, 993};
    std::vector<double> points_ = {75, 77, 80, 74, 72};

    double h=1e-6;
   
    int k = 10;
    Hermiteinterpolation hermite(points, f, points_, k);
    hermite.solve();
    std::cout << "Displacement =" << hermite.solve() << std::endl;
    Hermiteinterpolation hermite1(points, f, points_, k+h);
    double D = (hermite1.solve()-hermite.solve())/h;
    std::cout << "Speed =" << D << std::endl;

    return 0;
}