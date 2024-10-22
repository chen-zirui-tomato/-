#include "Function.h"




class Hermiteinterpolation : public Interpolation{
private:
/**
 * @param
 */
std::vector<double> values_;

public:
    Hermiteinterpolation(std::vector<double> points, std::vector<double> values,std::vector <double> values_, double x) 
                        : Interpolation(points, values, x) ,values_(values_){};
    void generatediffers() override;
};

void Hermiteinterpolation::generatediffers(){
    differs.resize(n, std::vector<double>(n, 0.0));
 
    for(int i = 0 ;i < n;i++){
        differs[i][0] = values[i] ;
    }
    int s = 0;
    for(int i = 1;i < n;i++){
        if(points[i] != points[i-1]){
                differs[i][1] = (differs[i][0] - differs[i-1][0])
                                    /(points[i]-points[i-1]);
        }
        else{
            differs[i][1] = values_[s]; 
            ++s;
        }
    }
    for(int j = 2;j < n;j++){
        for( int i = j; i < n; i++){
            differs[i][j] = (differs[i][j-1] - differs[i-1][j-1])
                            /(points[i]-points[i-j]);
        }
    }
}
