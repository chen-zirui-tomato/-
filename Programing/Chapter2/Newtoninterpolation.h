#include "Function.h"




class NewtonInterpolation : public Interpolation{
public:
    NewtonInterpolation(std::vector<double> points, std::vector<double> values, double x) 
                        : Interpolation(points, values, x) {};
    void generatediffers() override;
};

void NewtonInterpolation::generatediffers() {
    differs.resize(n, std::vector<double>(n, 0.0));
    for(int i = 0; i < n; i++){
        differs[i][0] = values[i];
    }
    for(int j = 1;j < n;j++){
        for( int i = j; i < n; i++){
            differs[i][j] = (differs[i][j-1] - differs[i-1][j-1])
                            /(points[i]-points[i-j]);
        }
    }
}
