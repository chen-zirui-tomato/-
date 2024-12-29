#include "ppForm.h"

ppForm::ppForm(const Function& f, double a, double b,int nodeCount, SplineType splineType)
        : f(f), nodeCount(nodeCount), splineType(splineType){
        //调用三次样条
        n = nodeCount - 1;
        if(a > b) std::swap(a,b);
        double h = (b-a)/n;
        for(int i = 0;i <= n;i++){
            nodes.push_back(a+i*h);
            values.push_back(f(a+i*h));
        }
        useFunction();
    };  

ppForm::ppForm(const Function& f, std::vector<double> nodes, SplineType splineType)
        : f(f), nodes(nodes), splineType(splineType){
        //任意节点以vector形式传入
        n = nodes.size() - 1;
        for (int i = 0; i <= n; i++) {
            values.push_back(f(nodes[i]));
        }
        useFunction();
    };

ppForm::ppForm(const Function& f, double a, double b,int nodeCount, bool isLinear)
        : f(f), nodeCount(nodeCount){
        //调用一次样条
        (void)isLinear;
        n = nodeCount - 1;
        if(a > b) std::swap(a,b);
        double h = (b-a)/n;
        for(int i = 0;i <= n;i++){
            nodes.push_back(a+i*h);
            values.push_back(f(a+i*h));
        }
    };

void ppForm::useFunction(){
        switch (splineType)
        {
        case SplineType::Natural:
            computeNaturalSpline();
            break;
        case SplineType::Complete:
            computeCompleteSpline();
            break;
        case SplineType::Periodic:
            computePeriodicSpline();
            break;
        }
    }

double ppForm::LinearInterpolation(double x){ //一次样条插值
        int i;
        for(i = 0; i <= n-1; i++) 
            if(nodes[i] <= x && x <= nodes[i+1]) break;
        double h = nodes[i+1] - nodes[i];
        double h1 = x - nodes[i];
        double h2 = nodes[i+1] - x;
        return (h1/h*values[i+1] + h2/h*values[i]);
    }

void ppForm::computeNaturalSpline(){ //计算自然三次样条 M(0)=0,M(n)=0
        //int n = nodeCount - 1 ;
        
        //miu.resize(nodeCount - 3); 
        // miu.resize(n-2);
        for(int i = 2;i <= n-1;i++)
            //接受i，返回miu(i)，存放在miu[i-2]中
            miu.push_back(calMiu(i));
        //是否需要反转miu？？

        //nabuda.resize(nodeCount - 3);
        // nabuda.resize(n-2);
        for(int i = 1;i <= n-2;i++)
            //接受i，返回nabuda(i)，存放在nabuda[i-1]中
            nabuda.push_back(calNabuda(i));

        //A.resize(nodeCount - 2); //A[0][0] = A(1,1) ~ A[n-2][n-2] = A(n-1,n-1)
        A.resize(n-1);
        for(int i = 0;i <= n-2 ; i++){ 
            A[i].resize(n-1);
            for(int j = 0;j <= n-2; j++){
                if (i == j) A[i][j] = 2.0;
                else if (i == j-1) A[i][j] = nabuda[i];
                else if (i == j+1) A[i][j] = miu[j];
                else A[i][j] = 0.0;
            }
        }
        // for(int i = 0;i <= n-2 ; i++) {std::cerr<<std::endl;
        // for(int j = 0;j <= n-2; j++) std::cerr << A[i][j] << " ";}
        // std::cerr<<std::endl;

        //b.resize(nodeCount - 2)
        // b.resize(n-1);
        for(int i = 1; i <= n-1; i++)
            //接受i,返回di,存放在b[i-1]中
            b.push_back(calCubicDivided(i));
        // std::reverse(b.begin(),b.end());
        // for(int i = 0;i <= n-2 ; i++) std::cerr << b[i] << " " << std::endl;
        //M.resize(nodeCount);
        // M.resize(n+1);
        M.push_back(0.0); //M(n）= 0
        caltridiagonal();
        M.push_back(0.0); //M(0) = 0
        std::reverse(M.begin(),M.end());
        // for(int i=0;i<=n;i++) std::cerr <<"M["<<i<<"]="<< M[i] << " ";
        // std::cerr<<std::endl;  
    } 

void ppForm::computeCompleteSpline(){ //计算完全三次样条

        // miu.resize(n);
        for(int i = 1; i <=n-1; i++) 
            miu.push_back(calMiu(i));
        miu.push_back(1.0);

        // nabuda.resize(n);
        nabuda.push_back(1.0);
        for(int i = 1; i <= n-1; i++)
            nabuda.push_back(calNabuda(i));

        A.resize(n+1);
        for(int i =0; i <= n ;i++){
            A[i].resize(n+1);
            for(int j = 0; j <= n; j++){
                if(i == j) A[i][j] = 2.0;
                else if(i == j-1) A[i][j] = nabuda[i];
                else if(i == j+1) A[i][j] = miu[j];
                else A[i][j] = 0.0;
            }
        }

        // b.resize(n+1);
        double temp = 6/(nodes[1] - nodes[0])*(calQuadraticDivided(0)
                                            - f.derivative(nodes[0]));
        b.push_back(temp);
        for(int i = 1; i<=n ;i++)
            b.push_back(calCubicDivided(i));
        b.push_back(6/(nodes[n] - nodes[n-1])*(f.derivative(nodes[n]) 
                                            - calQuadraticDivided(n-1)));

        // M.resize(n+1);
        caltridiagonal();
        std::reverse(M.begin(),M.end());
    }

void ppForm::computePeriodicSpline(){ //计算周期三次样条

        // miu.resize(n);
        //miu1存在miu[0]中
        for(int i = 1; i <= n; i++)
            miu.push_back(calMiu(i));
        
        // nabuda.resize(n);
        for(int i = 1; i <= n; i++)
            nabuda.push_back(calNabuda(i));
        
        A.resize(n);
        for(int i = 0; i <= n-1; i++){
            A[i].resize(n);
            for(int j = 0; j <= n-1; j++){
                if(i == j) A[i][j] = 2.0;
                else if(i == j-1) A[i][j] = nabuda[i];
                else if(i == j+1) A[i][j] = miu[j];
                else A[i][j] = 0.0;
            }
        }
        A[0][n-1] = miu[0];
        A[n-1][0] = nabuda[n-1];

        // b.resize(n);
        for(int i = 1; i <= n; i++)
            b.push_back(calCubicDivided(i));
        
        // M.resize(n);
        caltridiagonal();
        std::reverse(M.begin(),M.end());
    }

double ppForm::calMiu(int i){ //计算miu，接受i，返回miu(i)
        return (nodes[i] - nodes[i-1])/(nodes[i+1] - nodes[i-1]);
    }

double ppForm::calNabuda(int i){ //计算nabuda，接受i，返回nabuda(i)
        return (nodes[i+1] - nodes[i])/(nodes[i+1] - nodes[i-1]);
    }

double ppForm::calCubicDivided(int i){ //计算三阶差商, 接受i，返回b[i] = d(i+1)
        return 6/(nodes[i+1] - nodes[i-1])*(calQuadraticDivided(i) - calQuadraticDivided(i-1));
    }

double ppForm::calQuadraticDivided(int i){ //计算二阶差商, 接受i，返回f[xi+1,xi]
        // std::cerr << "\\"<<(values[i+1] - values[i])/(nodes[i+1] - nodes[i]) <<std::endl;
        return (values[i+1] - values[i])/(nodes[i+1] - nodes[i]);
    }

void ppForm::caltridiagonal(){// 计算三次样条的三对角矩阵
        for(int i = 1 ; i <= n-2 ; i++){ //从第二行开始改变主对角线的值,i表示第i+1行，
        //miu【0】=miu(2), nabuda【0】=nabuda(1)
            A[i][i] = A[i][i] - nabuda[i-1]*miu[i-1]/A[i-1][i-1];
        // for(int i = 0;i <= n-2 ; i++) {std::cerr<<std::endl;
        // for(int j = 0;j <= n-2; j++) std::cerr << A[i][j] << " ";}
        //     std::cerr<<std::endl;
        }
        double temp = b[n-2]/A[n-2][n-2];
        M.push_back(temp); //M(n-1)
        for(int i = n-3 ; i >= 0; i--) {//求解M，后插进M，再逆置
            temp = (b[i] - nabuda[i]*temp)/A[i][i];
            M.push_back(temp);
        }
    }

double ppForm::getValue(double x){ //插值函数的计算
        int i;
        //nodes[0] = x0, nodes[n] = nodes[nodeCount-1]= xn 
        for(i=0; i <= n-1 ;i++) if(nodes[i] <= x && x <= nodes[i+1]) break;
        double h = nodes[i+1] - nodes[i];
        double h1 = nodes[i+1] - x;
        double h2 = x - nodes[i];
        return pow(h1,3)/6/h*M[i] + pow(h2,3)/6/h*M[i+1]
                +(values[i] - M[i]/6*pow(h,2))*h1/h
                +(values[i+1] - M[i+1]/6*pow(h,2))*h2/h;
        // double result = 1/6/h*(M[i]*pow(nodes[i+1],3) - M[i+1]*pow(nodes[i],3)
        //                         + (6*values[i] - M[i]*pow(h,2))*nodes[i+1] 
        //                         - (6*values[i+1] - M[i+1]*pow(h,2))*nodes[i]);
        // result += x/6/h*(3*pow(nodes[i],2)*M[i+1] - 3*pow(nodes[i+1],2)*M[i] 
        //                         + 6*values[i+1] - M[i+1]*pow(h,2) 
        //                         - 6*values[i] + M[i]*pow(h,2));
        // result += pow(x,2)/2/h*(M[i]*nodes[i+1] - M[i+1]*nodes[i]);
        // return result += pow(x,3)/6/h*(M[i+1] - M[i]);
    };