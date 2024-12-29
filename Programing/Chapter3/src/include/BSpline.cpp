#include "BSpline.h"

BSpline::BSpline(int k,const Function& f, double a, double b, int nodeCount, SplineType splineType)
        : k(k), f(f), nodeCount(nodeCount), splineType(splineType){
            n = nodeCount - 1;
            if(a > b) std::swap(a,b);
            double h = (b-a)/n;
            for(int i = 0;i <= n;i++){
                nodes.push_back(a+i*h);
                values.push_back(f(a+i*h));
            }
            //初始化节点
            for(int i = k;i >= 1; i--)
                // t.push_back(a);
                t.push_back(a - i*h);
            for(int i = 0; i <= n; i++)
                t.push_back(nodes[i]);
            for(int i = 1; i <= k; i++)
                // t.push_back(b);
                t.push_back(b + i*h);
            // for(size_t i = 0; i < t.size(); i++) std::cerr << t[i] << " ";
            // std::cerr << std::endl;
            constructMAtrix();
        };

BSpline::BSpline(int k, const Function& f, std::vector<double> nodes, SplineType splineType)
        : k(k), f(f), nodes(nodes), splineType(splineType){
            n = nodes.size() - 1;
            double h = (nodes[n] - nodes[0])/n;
            for (int i = 0; i <= n; i++) {
                values.push_back(f(nodes[i]));
            }
            for(int i = k;i >= 1; i--)
                t.push_back(nodes[0] - i*h);
            for(int i = 0; i <= n; i++)
                t.push_back(nodes[i]);
            for(int i = 1; i <= k; i++)
                t.push_back(nodes[n] + i*h);
            constructMAtrix();
        };
    
double BSpline::constructBSpline(int i,int k, double x){ //构造bsplineiN，下标i，上标k, 参数x
        if(k == 0)  
            if(t[i-1] < x && x <= t[i]) return 1.0;
            else return 0.0;     
        else
            return (x - t[i-1])/(t[i + k-1] - t[i-1])*constructBSpline(i, k-1, x) 
                    + (t[i+k] - x)/(t[i+k] - t[i])*constructBSpline(i+1, k-1, x); 
    }

void BSpline::constructMAtrix(){
        //n+1行，n+N列
        A.resize(n+k);
        for(int i = 0; i <= n+k-1; i++)
            A[i].resize(n+k);
        for(int i = 0; i <= n+k-1; i++) for(int j = 0; j <= n+k-1; j++) {
            if(i == 0 && j == n+k-1) continue;
            if(i <= j && j <= i + k - 1 && i <= n ) 
                A[i][j] = constructBSpline(j+1, k, t[k+i]);
            else A[i][j] = 0.0;
        }
        addrange();
        solve();
    }

void BSpline::addrange(){ //增加边值条件
    if(k == 3)
        switch (splineType)
        {
        case SplineType::Natural:
            addNatrualCondition();
            break;
        case SplineType::Complete:
            addCompleteCondition();
            break;
        case SplineType::Periodic:
            addPeriodicCondition();
            break;
        }
    if(k == 2)
        addEndCondition();
    if(k == 7)
        addArbitraryCondition();
    }

void BSpline::addNatrualCondition(){
        for(int i = n+1; i <= n+k-1; i++){
            values.push_back(0.0);
        }
        std::cerr << std::endl;
        A[n+1][0] = t[4] - t[1];
        A[n+1][1] = 1/(t[1] - t[4]) + 1/(t[2] - t[5]);
        A[n+1][2] = t[5] - t[2];
        A[n+2][n+k-3] = t[n+4] - t[n+1];
        A[n+2][n+k-2] = 1/(t[n+1] - t[n+4]) + 1/(t[n+2] - t[n+5]);
        A[n+2][n+k-1] = t[n+5] - t[n+2];
        
    }

void BSpline::addCompleteCondition(){
    values.push_back(f.derivative(nodes[0]));
    values.push_back(f.derivative(nodes[n]));
    A[n+1][0] = -3*constructBSpline(2,2,nodes[0])/(t[3] - t[1]);
    A[n+1][1] = 3*(constructBSpline(2,2,nodes[0])/(t[3] - t[1]) 
                 - constructBSpline(3,2,nodes[0])/(t[4] - t[2]));
    A[n+1][2] = 3*(constructBSpline(3,2,nodes[0])/(t[4] - t[2]) );
    A[n+2][n+k-3] = -3*constructBSpline(n+2,2,nodes[n])/(t[n+3] - t[n+1]);
    A[n+2][n+k-2] = 3*(constructBSpline(n+2,2,nodes[n])/(t[n+3] - t[n+1]) 
                     - constructBSpline(n+3,2,nodes[n])/(t[n+4] - t[n+2]));
    A[n+2][n+k-1] = 3*(constructBSpline(n+3,2,nodes[n])/(t[n+4] - t[n+2]));
}

void BSpline::addPeriodicCondition(){
    values.push_back(0);
    values.push_back(0);
    A[n+1][0] = -3*constructBSpline(2,2,nodes[0])/(t[3] - t[1]);
    A[n+1][1] = 3*(constructBSpline(2,2,nodes[0])/(t[3] - t[1]) 
                 - constructBSpline(3,2,nodes[0])/(t[4] - t[2]));
    A[n+1][2] = 3*(constructBSpline(3,2,nodes[0])/(t[4] - t[2]) );
    A[n+1][n+k-3] = 3*constructBSpline(n+2,2,nodes[n])/(t[n+3] - t[n+1]);
    A[n+1][n+k-2] = -3*(constructBSpline(n+2,2,nodes[n])/(t[n+3] - t[n+1]) 
                     - constructBSpline(n+3,2,nodes[n])/(t[n+4] - t[n+2]));
    A[n+1][n+k-1] = -3*(constructBSpline(n+3,2,nodes[n])/(t[n+4] - t[n+2]));

    A[n+2][0] = t[4] - t[1];
    A[n+2][1] = 1/(t[1] - t[4]) + 1/(t[2] - t[5]);
    A[n+2][2] = t[5] - t[2];
    A[n+2][n+k-3] = -t[n+4] + t[n+1];
    A[n+2][n+k-2] = -1/(t[n+1] - t[n+4]) - 1/(t[n+2] - t[n+5]);
    A[n+2][n+k-1] = -t[n+5] + t[n+2];
}

void BSpline::addEndCondition(){
    values.push_back(f(nodes[0]-0.5)+f(nodes[0]+0.5));
    A[n+1][0] = 0.2;
    A[n+1][1] = 0.8;
    A[n+1][n+k-2] = 0.8;
    A[n+1][n+k-1] = 0.2;
}

void BSpline::addArbitraryCondition(){
    values.push_back(2);
    values.push_back(2);
    A[n+1][0] = 1;
    A[n+1][1] = 0;
    A[n+1][2] = 0;
    A[n+1][n+k-3] = 0;
    A[n+1][n+k-2] = 0;
    A[n+1][n+k-1] = 1;

    A[n+2][0] = 0;
    A[n+2][1] = 1;
    A[n+2][2] = 0;
    A[n+2][n+k-3] = 1;
    A[n+2][n+k-2] = 0;
    A[n+2][n+k-1] = 0;

    /////////////////////////////////////////////////////////////////////////
    values.push_back(f.derivative(nodes[0]));
    values.push_back(f.derivative(nodes[n]));
    A[n+3][0] = 0;
    A[n+3][1] = 0 ;
    A[n+3][2] = 1;
    A[n+4][n+k-3] = 0;
    A[n+4][n+k-2] = 0;
    A[n+4][n+k-1] = 1;

    /////////////////////////////////////////////////////////////////////////////
    values.push_back(1);
    values.push_back(1);
    A[n+5][0] = 0.4;
    A[n+5][1] = 0.4;
    A[n+5][2] = 0.3;
    A[n+6][n+k-3] = 0.3;
    A[n+6][n+k-2] = 0.4;
    A[n+6][n+k-1] = 0.4;
}

void BSpline::solve() {
        Eigen::MatrixXd eigenA(n+k, n+k);
        for (size_t i = 0; i < A.size(); ++i) 
            for (size_t j = 0; j < A[0].size(); ++j) 
                eigenA(i, j) = A[i][j];
        
        Eigen::VectorXd eigenValues(n+k);
        for (size_t i = 0; i < values.size(); ++i) 
            eigenValues(i) = values[i];
        
        // Eigen::VectorXd P_solution = eigenA.colPivHouseholderQr().solve(eigenValues);
        Eigen::VectorXd P_solution = eigenA.householderQr().solve(eigenValues);
        P.resize(P_solution.size());
        for (int i = 0; i < P_solution.size(); ++i) 
            P[i] = P_solution(i);

        // for(size_t i = 0; i < A.size(); i++) {
        //     std::cerr << std::endl;
        //     for(size_t j = 0; j < A[0].size(); j++) std::cerr << A[i][j] << " ";
        // }
        // std::cerr << std::endl;
        // for(size_t i = 0; i < values.size(); i++)
        //     std::cerr << values[i] << " ";
        // std::cerr << std::endl;
        // for(size_t i = 0; i < P.size(); i++)
        //     std::cerr << P[i] << " ";
        // std::cerr << std::endl;
    }

double BSpline::getValue(double x){
        double sum = 0.0;
        int i;
        for(i = 0; i <= n+k-1; i++) if(t[i] <= x && x <= t[i+1]) break;
        for(int j = 2; j >= -1; --j)
            sum += constructBSpline(i-j, k, x)*P[i-j-1];
        return sum;
    }
