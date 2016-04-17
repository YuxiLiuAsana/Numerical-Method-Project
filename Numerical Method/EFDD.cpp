//
//  EFDD.cpp
//  Finite Different Method
//
//  Created by LiuYuxi on 12/13/15.
//  Copyright © 2015 LiuYuxi. All rights reserved.
//

#include "EFDD.hpp"
EFDD::EFDD():FD()
{

}

EFDD::~EFDD()
{

}

EFDD::EFDD(double _S0, double _K, double _r,double _q, double _T,double _sigma, int _M, double _alpha, double _T_div)
{
    S0 = _S0;
    K = _K;
    r = _r;
    q = _q;
    T = _T;
    sigma = _sigma;
    M = _M;
    M1 = M;
    alpha = _alpha;
    T_div = _T_div;
    t_div = (T-T_div)*sigma*sigma/2.0;
    x_compute = log(S0/K)+log(1-q);
    t_final = T * sigma * sigma/2.0;
    alpha1 = alpha;
    double dt1 = t_div/M;
    double dx = sqrt(dt1/alpha);
    x_left = log(S0/K) + (r - sigma*sigma/2.0)*T - 3 * sigma * sqrt(T);
    x_right = log(S0/K)+ (r - sigma*sigma/2.0)*T + 3 * sigma * sqrt(T);
    //cout<<x_left<<","<<x_right<<endl;
    N1 = ceil((x_compute - x_left)/dx);
    //cout<<N1<<endl;
    N2 = ceil((x_right-x_compute)/dx);
    N = N1 + N2;
    x_left = x_compute - N1 * dx;
    x_right = x_compute + N2 * dx;
    a = r/sigma/sigma - 0.5;
    b = pow((r/sigma/sigma + 0.5),2);
    
    
    double dt2 = alpha * dx * dx;
    M2 = ceil((t_final-t_div)/dt2);
    dt2 = (t_final-t_div)/(double) M2;
    alpha2 = dt2/dx/dx;

}

//边界条件修改
double EFDD::f(int i)
{
    double x = x_left + (x_right - x_left)/(double)N * (double)i;
    return K * exp(a * x)*(exp(x)-1) * (exp(x)>1);
}


double EFDD::g_left_old(int i)
{
    return 0;
}

double EFDD::g_right_old(int i)
{
    double t = t_div / (double )M * i;
    return K * exp(a * x_right + b * t) * (exp(x_right)-exp(-2 * r * t / sigma/sigma));
}

double EFDD::g_left_new(int i)
{
    return 0;
}

double EFDD::g_right_new(int i )
{
    double t = (t_final - t_div) / M2 * i;
    return K * exp(a * (x_right - log(1-q)) + b * t) * (exp(x_right - log(1-q)) - exp(-2 * r * t / sigma /sigma));
}

void EFDD::ForwardEuler_old()
{
    M1 = M;
    alpha = alpha1;
    data.clear();
    vector<double> m;
    for(int i = 0 ; i < N+1; i ++)
        m.push_back(0);
    for(int i = 0; i < M+1; i++)
        data.push_back(m);
    for(int i = 0; i < N+1; i++)
        data[0][i]=f(i);
    for(int i = 0 ; i < M+1; i++)
    {
        data[i][0] = g_left_old(i);
        data[i][N] = g_right_old(i);
    }
    ForwardEuler();
}

void EFDD::BackwardEuler_old(string method)
{
    M1 = M;
    alpha = alpha1;
    data.clear();
    vector<double> m;
    for(int i = 0 ; i < N+1; i ++)
        m.push_back(0);
    for(int i = 0; i < M+1; i++)
        data.push_back(m);
    for(int i = 0; i < N+1; i++)
        data[0][i]=f(i);
    for(int i = 0 ; i < M+1; i++)
    {
        data[i][0] = g_left_old(i);
        data[i][N] = g_right_old(i);
    }
    BackwardEuler(method);
}

void EFDD::CrankNicolson_old(string method)
{
    M1 = M;
    alpha = alpha1;
    data.clear();
    vector<double> m;
    for(int i = 0 ; i < N+1; i ++)
        m.push_back(0);
    for(int i = 0; i < M+1; i++)
        data.push_back(m);
    for(int i = 0; i < N+1; i++)
        data[0][i]=f(i);
    for(int i = 0 ; i < M+1; i++)
    {
        data[i][0] = g_left_old(i);
        data[i][N] = g_right_old(i);
    }
    CrankNicolson(method);
}

void EFDD::ForwardEuler_new()
{
    ForwardEuler_old();
    alpha = alpha2;
   
    vector<double> b  = data[M];
    M = M2;
    data.clear();
    for(int i = 0; i < M2 + 1; i++)
    {
        data.push_back(b);
        data[i][0] = g_left_new(i);
        data[i][N] = g_right_new(i);
    }
    ForwardEuler();
   
}
void EFDD::BackwardEuler_new(string method)
{
    BackwardEuler_old(method);
    alpha = alpha2;
    
    vector<double> b  = data[M];
    M = M2;
    data.clear();
    for(int i = 0; i < M2 + 1; i++)
    {
        data.push_back(b);
        data[i][0] = g_left_new(i);
        data[i][N] = g_right_new(i);
    }
    BackwardEuler(method);

}
void EFDD::CrankNicolson_new(string method)
{
    CrankNicolson_old(method);
    alpha = alpha2;
    
    vector<double> b  = data[M];
    M = M2;
    data.clear();
    for(int i = 0; i < M2 + 1; i++)
    {
        data.push_back(b);
        data[i][0] = g_left_new(i);
        data[i][N] = g_right_new(i);
    }
    CrankNicolson(method);

}

double EFDD::u_value()
{
    return data[M2][N1];
}

double EFDD::Option_value()
{
    
    return exp(-a * log(S0/K)- b * t_final) * data[M2][N1];
}

double EFDD::Delta()
{
    double dx = (x_right-x_left)/N;
    double v1 = exp(-a * (log(S0/K) + dx) - b * t_final) * data[M2][N1 + 1];
    double v_1 =exp(-a * (log(S0/K) - dx) - b * t_final) * data[M2][N1 - 1];
    double s1 = K * exp(log(S0/K)+dx);
    double s_1 = K * exp(log(S0/K)-dx);
    return (v1 - v_1)/(s1 - s_1);
}

double EFDD::Gamma()
{
    double dx = (x_right-x_left)/N;
    double v1 = exp(-a * (log(S0/K) + dx) - b * t_final) * data[M2][N1 + 1];
    double v_1 =exp(-a * (log(S0/K) - dx) - b * t_final) * data[M2][N1 - 1];
    double v0 = exp(-a * (log(S0/K) ) - b * t_final) * data[M2][N1 ];
    double s1 = K * exp(log(S0/K)+dx);
    double s_1 = K * exp(log(S0/K)-dx);
    return ((S0-s_1)*v1 - (s1-s_1)* v0 + (s1 - S0)* v_1)/((S0 - s_1)* (s1 - S0)* ((s1 - s_1)/2.0));
}

double EFDD::Theta()
{
    double dt = (t_final-t_div)/(double)M2;
    return -(Option_value()-exp(-a * log(S0/K)- b * (t_final-dt)) * data[M2-1][N1])/(2*dt/sigma/sigma);
}