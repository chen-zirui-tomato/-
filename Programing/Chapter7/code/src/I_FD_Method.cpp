#include "I_FD_Method.h"

I_FD_Method::I_FD_Method(int m, DynamicFunction& f, DynamicFunction& b, DynamicFunction& gx, DynamicFunction& gy,
     DynamicFunction& centerFunc, std::string condition_type, std::vector<double> center, double radius):
    FD_Method::FD_Method(m, f, b, gx, gy, condition_type), centerFunc(centerFunc), delta(0.01*h), radius(radius)
    {
        // code here
        this->center[0] = center[0];
        this->center[1] = center[1];
        
    }

void I_FD_Method::re_construct_equations(){
    for(int i=0; i<m+2; i++) for(int j=0; j<m+2; j++) {
        double p[2] = {x_points[i], y_points[j]};
        if(get_Length(p, this->center) < radius){
            mark_in[std::make_pair(i,j)] = true;
            for(int k = 0; k < m*m; k++) 
                A((j-1)*m+i-1,k) = 0, A(k, (j-1)*m+i-1) = 0;
            f_values[i][j] = 0;
            }
        else if(radius < get_Length(p, this->center) < radius + h) 
            mark_bound[std::make_pair(i,j)] = true;
        }
    for(int i = 0; i < m+2; i++) for(int j = 0; j < m+2; j++) if(mark_bound[std::make_pair(i,j)])
        for(int dx : {-1,0,1}) for(int dy : {-1,0,1}) 
            if(abs(dx + dy) == 1 && mark_in[std::make_pair(i+dx,j+dy)]){
                Sign sign1 = Sign::null, sign2 = Sign::null;
                if(dy == 1) (sign1 == Sign::null) ? sign1 = Sign::up : sign2 = Sign::up;
                else if(dy == -1) (sign1 == Sign::null) ? sign1 = Sign::down : sign2 = Sign::down;
                else if(dx == +1) (sign1 == Sign::null) ? sign1 = Sign::right : sign2 = Sign::right;
                else if(dx == -1) (sign1 == Sign::null) ? sign1 = Sign::left : sign2 = Sign::left;
                re_construct_RHS(i,j,sign1,sign2);
            }
    }



void I_FD_Method::re_construct_RHS(int i, int j, Sign sign1, Sign sign2){
    if(sign2 == Sign::null) judge(i, j, sign1);
    else{
        judge(i, j, sign1);
        judge(i, j, sign2);
    }
}

void I_FD_Method::judge(int i, int j, Sign sign){
    double y;
    switch(sign){
        case Sign::up:
            findBoundaryNeighbor(i, j, static_cast<double>(j), 1, sign); //u_ik,u_i,j+1处的边界条件值为y
            break;
        case Sign::down:
            findBoundaryNeighbor(i, j, static_cast<double>(j), -1, sign); //ik
            break;
        case Sign::right:
            findBoundaryNeighbor(i, j, static_cast<double>(i), 1, sign); //kj
            break;
        case Sign::left:
            findBoundaryNeighbor(i, j, static_cast<double>(i), -1, sign); //kj
            break;
        default:
            throw std::invalid_argument("sign should not be null");
    }
}

void I_FD_Method::findBoundaryNeighbor(int i, int j, double k, int sign_, Sign sign){
    double length = 0;
    bool flag = false;
    double p[2] = {x_points[i], y_points[j]};
    
    if(get_Length(p, this->center) < radius) throw std::invalid_argument("p should not be in the circle");
    if(static_cast<double> (i) == k ) flag = true;
    do{
        length += delta*sign_;
        if(flag) p[0] += delta*sign_; 
        else p[1] += delta*sign_;

        if(get_Length(p, this->center) <= radius) break;
    }while(1);
    length = fabs(length);
    double values = (1.0-h/length)*f(x_points[i], y_points[j])
                    +h/length*b(p[0], p[1]);
    switch(condition_type)
    {
        case ConditionType::Dirichlet:
          
        build_Dirichlet_condition(i, j, values, sign);
            break;
        case ConditionType::Neumann:
            build_Neumann_condition(i, j, length, sign);
            break;
    }
}

void I_FD_Method::build_Dirichlet_condition(int i, int j, double y, Sign sign){
    f_values[i][j] += y;
}

void I_FD_Method::build_Neumann_condition(int i, int j, double length, Sign sign){
    double cos, sin;
    double a[2];
    double ax, ay, bx, by, cx, cy, Ldx, Ldy;
    switch(sign){
        case Sign::up:
            cos = (y_points[j] - center[1] + length)/radius;
            sin = (x_points[i] - center[0])/radius;
        
            if(!mark_in[std::make_pair(i-1,j)]) 
                ax = -h, a[0] = x_points[i-1];
            else ax = h, a[0] = x_points[i+1];
            ay = -length, a[1] = y_points[j];
            bx = 0, by = -length, cx = 0, cy = 1-length;
            Ldx = x_points[i], Ldy = y_points[j] + length;
            break;
        case Sign::down:
            cos = (y_points[j] - center[1] - length)/radius;
            sin = (x_points[i] - center[0])/radius;
            if(!mark_in[std::make_pair(i-1,j)]) 
                ax = -h, a[0] = x_points[i-1];
            else ax = h, a[0] = x_points[i+1];
            ay = length, a[1] = y_points[j];
            bx = 0, by = length, cx = 0, cy = length-1;
            Ldx = x_points[i], Ldy = y_points[j] - length;
            break;
        case Sign::right:
            cos = (y_points[j] - center[1])/radius;
            sin = (x_points[i] - center[0] + length)/radius;
            if(!mark_in[std::make_pair(i,j-1)]) 
                ay = h, a[1] = y_points[j-1];
            else ay = -h, a[1] = y_points[j+1];
            ax = -length, a[0] = x_points[i];
            bx = -length, by = 0, cx = 1-length, cy = 0;
            Ldx = x_points[i] + length, Ldy = y_points[j];
            break;
        case Sign::left:
            cos = (y_points[j] - center[1])/radius;
            sin = (x_points[i] - center[0] - length)/radius;
            if(!mark_in[std::make_pair(i,j-1)]) 
                ay = h, a[1] = y_points[j-1];
            else ay = -h, a[1] = y_points[j+1];
            ax = length, a[0] = x_points[i];
            bx = length, by = 0, cx = length-1, cy = 0;
            Ldx = x_points[i] - length, Ldy = y_points[j];
            break;
    }
    Eigen::Matrix<double, 3, 3> G;
    G << 1, 1, 1,
        ax, bx, cx,
        ay, by, cy;
    Eigen::Vector<double, 3> Gb;
    Gb << 0, cos, sin;
    Eigen::Vector<double, 3> Ga = G.partialPivLu().solve(Gb);
    
    f_values[i][j] += (centerFunc(Ldx, Ldy)-Ga[0]*f(a[0], a[1])-Ga[1]*f(x_points[i], y_points[j]))/Ga[2];

}

void I_FD_Method::I_solver(){
    construct_A_matrix();
    construct_RHS_vector();
    re_construct_equations();
    solver();
}
