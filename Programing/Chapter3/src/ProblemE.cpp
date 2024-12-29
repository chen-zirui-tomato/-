#include "include/BSpline.cpp"
#include "include/Point2d.h"
#include <fstream>
// #include <boost/math/quadrature/trapezoidal.hpp>

double map(double x){
    return 1/(1+exp(-4*x));
};


int main(){
    double pi = acos(-1.0);
    std::ofstream outfile1("bin/testE1.dat");
    std::ofstream outfile2("bin/testE2.dat");
    std::ofstream outfile3("bin/testE3.dat");
    std::ofstream outfile4("bin/testE4.dat");
    std::ofstream outfile5("bin/testE5.dat");
    std::ofstream outfile6("bin/testE6.dat");
    std::ofstream outfile7("bin/testE7.dat");
    double dx = 0.01;

    // class f1 :public Function{
    // public:
    //     double operator()(double x) const{
    //         return 2.0/3.0*(sqrt(3-x*x) + sqrt(fabs(x)));
    //     }
    // };
    // class f2 :public Function{
    // public:
    //     double operator()(double x) const{
    //         return 2.0/3.0*(-sqrt(3-x*x) + sqrt(fabs(x)));
    //     }
    // };

    // BSpline obj1(3, f1(), -1.7, 1.7, 160, SplineType::Natural);
    // for(double i = xStart+dx; i < xEnd; i += dx){
    //     outfile1 << i << " " << obj1.getValue(i) << std::endl;
    // }
    // BSpline obj2(3, f2(), -1.7, 1.7, 160, SplineType::Natural);
    // for(double i = xStart+dx; i < xEnd; i += dx){
    //     outfile2 << i << " " << obj2.getValue(i) << std::endl;
    // }

    /////////////////////////////////////////////////
    class f1x:public Function{
    public:
        double operator()(double x) const{
            return sqrt(3.0)*sin(x);
        }
    };
    class f1y:public Function{
    public:
        double operator()(double x) const{
            return 2.0/3.0*(sqrt(3.0)*cos(x) + sqrt(sqrt(3.0)*fabs(sin(x))));
        }
    };
    BSpline obj1x(3, f1x(), 0, pi, 40, SplineType::Complete);
    BSpline obj1y(3, f1y(), 0, pi, 40, SplineType::Complete);
    for(double i = 0.0+dx; i <= pi; i += dx)
        outfile1 << obj1x.getValue(i) << " " << obj1y.getValue(i) << std::endl;
        // outfile1 << f1x()(i) << " " << f1y()(i) << std::endl;
    for(double i = 0.0+dx; i <= pi; i += dx)
        outfile2 << -obj1x.getValue(i) << " " << obj1y.getValue(i) << std::endl;
        // outfile2 << -f1x()(i) << " " << f1y()(i) << std::endl;
    //////////////////////////////////////////////

    // //做cumulative chordal插值点
    // double result = boost::math::quadrature::trapezoidal(f, 0.0, pi-0.2);
    // std::cout << "积分结果: " << result << std::endl;
    // std::vector<Point2d> points;
    // for(int i = 0; i<= N; ++i){
    //     double t = i*pi/N;
    //     points.push_back(Point2d(f1x()(t), f1y()(t)));
    // }
    int N = 40;
    std::vector<double> t;
    t.push_back(0.0);
    for(int i = 1; i<= N; ++i){
        t.push_back(pi*map(pi*i/N - (pi)/2));
    }   
    t.push_back(pi);
    std::vector<double> ti;
    for(int i = 0; i<= N; ++i){
        ti.push_back(t[i]);
    }
    BSpline obj2x(3, f1x(), ti, SplineType::Complete);
    BSpline obj2y(3, f1y(), ti, SplineType::Complete);
    for(double i = 0.0+dx; i < pi-dx; i += dx)
        outfile3 << obj2x.getValue(i) << " " << obj2y.getValue(i) << std::endl;
    for(double i = 0.0+dx; i < pi-dx; i += dx)
        outfile4 << -obj2x.getValue(i) << " " << obj2y.getValue(i) << std::endl;

    //////////////////////////////////////////////////////////////
    class f2x:public Function{
    public:
        double operator()(double x) const{
            return sin(x) + x*cos(x);
        }
    };
    class f2y:public Function{
    public:
        double operator()(double x) const{
            return cos(x) - x*sin(x);
        }
    };
    class f2s:public Function{
    public:
        double operator()(double x) const{
            return 2.0*sqrt(1+x*x/4)*x/2 
                + 2.0*log(sqrt(1+x*x/4)+x/2) ;
        }
        Function* operator+(double value) {
            class AddedFunction : public Function{
            public:
                const f2s& baseFunction;
                double addtionalValue;

                AddedFunction(const f2s& f, double v): baseFunction(f), addtionalValue(v) {}

                double operator()(double x) const override{
                    return baseFunction(x) + addtionalValue;
                }
            };
            return new AddedFunction(*this, value);
        }
    };
    
    BSpline obj3x(3, f2x(), 0.0, 6.0*pi, 160, SplineType::Complete);
    BSpline obj3y(3, f2y(), 0.0, 6.0*pi, 160, SplineType::Complete);
    for(double i = 0.0+dx; i <= 6.0*pi; i += dx)
        outfile5 << obj3x.getValue(i) << " " << obj3y.getValue(i) << std::endl;

    ////////////////////////////////////////////////////////////////
    int N2 = 160; //控制点个数
    // std::vector<Point2d> points2;
    // for(int i = 0; i<= N2; ++i){
    //     double t = i*6*pi/N2;
    //     points2.push_back(Point2d(f2x()(t), f2y()(t)));
    // }
    std::vector<double> t2;
    t2.push_back(0.0);
    double length = f2s()(6*pi);
    double interval = length/N2;
    // std::cout << "Length" << f2s()(6*pi) << std::endl;
    for(int i = 1; i< N2; ++i){
        f2s function;
        Function* newfunction = function + -interval*i;
        double temp = newfunction->bisection_method(0.0, 6*pi);
        t2.push_back(temp);
        // std::cout << "t2[" << i << "] = " << temp << std::endl;
        delete newfunction;
    }


    BSpline obj4x(3, f2x(), t2, SplineType::Complete);
    BSpline obj4y(3, f2y(), t2, SplineType::Complete);
    for(double i = 0.1; i < pi*5.9; i += dx)
        outfile6 << obj4x.getValue(i) << " " << obj4y.getValue(i) << std::endl;
        // outfile6 << f2x()(i) << " " << f2y()(i) << std::endl;
    
    ///////////////////////////////////////////////////////////////////
    class f3x:public Function{
    public:
        double operator()(double x) const{
            return sin(cos(x))*cos(sin(x));
        }
    };
    class f3y:public Function{
        public:
        double operator()(double x) const{
            return sin(cos(x))*sin(sin(x));
        }
    };
    class f3z:public Function{
        public:
        double operator()(double x) const{
            return cos(cos(x));
        }
    };

    BSpline obj5x(3, f3x(), 0.0, 2*pi, 160, SplineType::Periodic);
    BSpline obj5y(3, f3y(), 0.0, 2*pi, 160, SplineType::Periodic);
    BSpline obj5z(3, f3z(), 0.0, 2*pi, 160, SplineType::Periodic);
    for(double i = 0.0+dx; i <= 2*pi; i += dx)
        outfile7 << obj5x.getValue(i) << " " << obj5y.getValue(i) << " " << obj5z.getValue(i) << std::endl;
        //  outfile7 << f3x()(i) << " " << f3y()(i) << " " << f3z()(i) << std::endl;

    // double length3 = 2*pi;
    // int N3 = 10;

    outfile1.close();
    outfile2.close(); 
    outfile3.close();
    outfile4.close();
    outfile5.close();
    outfile6.close();
    outfile7.close();
    return 0;
}