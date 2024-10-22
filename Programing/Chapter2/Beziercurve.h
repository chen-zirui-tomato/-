#include "Function.h"
#include <vector>
#include <cmath>

class eachBezierCurve  {
protected:
// 用于存放q[4]
std::vector<double> q;
// 用于存放区间范围
double a, b;

public:
    // 在构造函数中初始化父类
    eachBezierCurve(std::vector<double> points, std::vector<double> values, std::vector<double> values_ )
        :  a(points[0]), b(points[1]) {
            q.resize(4);
            // 初始化q,n = 3
            q[0] = (values[0]);
            q[1] = (values[0]+values_[0]/3.0);
            q[2] = (values[1]-values_[0]/3.0);
            q[3] = (values[1]);
            // std::cerr << "points[0] = " << points[0] << std::endl;
            // std::cerr << "values_[0] = " << values_[0] << std::endl;
            // std::cerr << "q[0] = " << q[0] << std::endl;
            // std::cerr << "q[1] = " << q[1] << std::endl;
            // std::cerr << "q[2] = " << q[2] << std::endl;
            // std::cerr << "q[3] = " << q[3] << std::endl;
        };

    //计算bezier曲线
    double solve(double x);

    //计算组合数
    int combination(int n, int k);

    //计算Bernstein多项式的值
    double bernstein(int n, int k, double t);

    //更新q
    void updateq(double x);

    //判断x是否在区间范围内
    bool isInsegment(double x);
    
    ~eachBezierCurve() {};
};


double eachBezierCurve::solve(double x) {
    // std::cerr << "use solve" << std::endl;
    // std::cerr << "x = " << x << std::endl;
    double t = (x - a) / (b - a);
    double eps = 1e-6;
    t = std::max(eps, std::min(1 - eps, t));

    double result = 0.0;
    // for(int i = 0; i < 4; ++i){
    //     std::cerr << "q[" << i << "] = " << q[i] << std::endl;
    // }
    updateq(t);
    for(int i = 0; i < 4; ++i){
        result += q[i] ;
    }
    // std::cerr << "resultsolve = " << result << std::endl << std::endl;

    return result;
}

int eachBezierCurve::combination(int n, int k) {
    // std::cerr << "use combination" << std::endl;
    if (k > n || k < 0) return 0; // 不合法的情况
    if (k > n - k) k = n - k; // C(n, k) = C(n, n-k)
    
    int result = 1;
    for (int i = 0; i < k; ++i) {
        result = result * (n - i) / (i + 1);
        // std::cerr << "resultcom = " << result << std::endl;
    }
    return result;
}

double eachBezierCurve::bernstein(int n, int k, double x) {
    // std::cerr << "calculation bernstein" << std::endl;
    double result = combination(n, k) * pow(x, k) * pow(1 - x, n - k);
    // std::cerr << "n = " << n << std::endl;
    // std::cerr << "k = " << k << std::endl;
    // std::cerr << "x = " << x << std::endl;
    // std::cerr << "resultbern = " << result << std::endl;
    return result;
}

void eachBezierCurve::updateq(double x) {
    for(int i = 0; i <= 3; ++i){
        // std::cerr << "update " << i << "times " << std::endl;
        q[i] = q[i] * bernstein(3,i,x);
        // std::cerr << "q[" << i << "] = " << q[i] << std::endl;

    }
}

bool eachBezierCurve::isInsegment(double x) {
    // std::cerr << "use isInsegment" << std::endl;
    return x >= a && x <= b;
}





class BezierCurve {
protected:
    // 用于存放控制点
    std::vector<double> points;
    // 用于存放函数值
    std::vector<double> values;
    // 用于存放导数值
    std::vector<double> values_;
    // 控制点数量 -1
    const int n ; 
    //计算曲线点
    // const double x;
    //用于存放eachBezierCurve
    std::vector<eachBezierCurve> eachBezierCurves;

public:
    BezierCurve(std::vector<double> points, std::vector<double> values, std::vector<double> values_, int n)
        : points(points), values(values), values_(values_) , n(n) {};        


    //对Bezier曲线分段求值
    void segment();

    //计算
    double solveP(double x);

    ~BezierCurve() {};
};

void BezierCurve::segment() {
    // std::cerr << "use segment" << std::endl;
    for(int i = 0; i < n; ++i){
        std::vector<double> p, v, v_;
        p.push_back(points[i]);
        p.push_back(points[i+1]);
        v.push_back(values[i]);
        v.push_back(values[i+1]);
        v_.push_back(values_[i]);
        eachBezierCurve segment(p, v, v_);
        eachBezierCurves.push_back(segment);
    }
}

double BezierCurve::solveP(double x) {
    // std::cerr<<"use solveP"<<std::endl;
    for(auto& segment : eachBezierCurves){
        if(segment.isInsegment(x)){
            return segment.solve(x);
        }
    }
    return 0.0; 
}





