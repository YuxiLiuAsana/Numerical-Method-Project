//
//  EFD.cpp
//  Finite Different Method
//
//  Created by LiuYuxi on 12/11/15.
//  Copyright © 2015 LiuYuxi. All rights reserved.
//

#include "EFD.hpp"
EFD::EFD(double _S0, double _K, double _r, double _q, double _T, double _sigma, int _M, double _alpha )
{
    
    S0 = _S0;
    K = _K;
    x_compute = log(S0/K);
    r = _r;
    q = _q;
    T = _T;
    sigma = _sigma;
    alpha = _alpha;
    t_final = sigma * sigma * T/2.0;
    
    
    //修改domain
    //************************************************************************************************
    x_left = log(S0/K) + (r - q - sigma * sigma / 2.0) * T - 3 * sigma * sqrt(T);
    x_right = log(S0/K) + (r - q - sigma * sigma / 2.0) * T + 3 * sigma * sqrt(T);
    //******************************************************************************************
    
    
    
    
    
    M = _M;
    N = floor((x_right - x_left)/sqrt(t_final/M/alpha));
    alpha = (t_final/M)/((x_right-x_left)*(x_right - x_left)/N/N);
    
    a = (r-q)/sigma/sigma - 0.5;
    b = pow((r-q)/sigma/sigma+0.5,2)+2.0 * q/sigma/sigma;
    data.clear();
    vector<double> s;
    for(int i = 0; i < N + 1; i ++)
        s.push_back(0);
    for(int i = 0; i < M + 1; i++)
        data.push_back(s);
    for(int i = 0; i < N + 1; i ++)
    {
        data[0][i] = f(i);
    }
    for(int j = 0; j < M + 1; j++)
    {
        data[j][0] = g_left(j);
        data[j][N] = g_right(j);
    }
}

EFD::EFD():FD()
{

}

EFD::~EFD()
{

}













//*****************************************************************************************************************

// 边界条件修改
double EFD::f(int i)
{
    double x = x_left + i * (x_right - x_left) / N;
    //call
    //return K * exp(a * x ) * (exp(x)-1) * (1 < exp(x));
    //put
    return K * exp(a * x) * (1 - exp(x)) * (1 > exp(x));
}

// 边界条件修改
double EFD::g_left(int j)
{
    double t = t_final/(double) M * j;
    //call
    //return 0;
    
    //put
    return K * exp(a * x_left + b * t)*(exp(- 2 * r * t/sigma /sigma) - exp(x_left - 2 * q * t / sigma/sigma));
}


//边界条件修改
double EFD::g_right(int j)
{
    double t = t_final / (double)M * (double)j;
    
    //call
    //return -K * exp(a * x_right + b * t)*(exp(- 2 * r * t/sigma /sigma) - exp(x_right - 2 * q * t / sigma/sigma));
    
    //put
    return 0;
}

//*****************************************************************************************************************































double EFD::ErrorPointwise1()
{
    int i = (log(S0/K)-x_left)/(x_right-x_left)*N;

    double xi = x_left + i * (x_right-x_left)/(double)N;

    double xi1 = x_left + (i+1) * (x_right-x_left)/(double)N;

    double si = K * exp(xi);

    double si1 = K * exp(xi1);

    double vi = exp(-a*xi - b * t_final)* data[M][i];

    double vi1 = exp(-a*xi1 - b * t_final) * data[M][i+1];
   // cout<<data[M][i]<<","<<data[M][i+1]<<endl;
    double va = ((si1-S0)*vi+(S0-si)*vi1)/(si1-si);
    //cout<<va;
    Option_Pricing temp(S0,K,r,q,T,sigma);
    double ve = temp.E_BSM_Put_Price();
    return abs(ve-va);
}

double EFD::ErrorPointwise2()
{
    int i = (log(S0/K)-x_left)/(x_right-x_left)*N;
    //cout<<i<<endl;
    double xi = x_left + i * (x_right-x_left)/N;
   // cout<<xi<<endl;
    double xi1 = x_left + (i+1) * (x_right-x_left)/N;
    
    double ul =((xi1-log(S0/K))*data[M][i]+(log(S0/K)-xi)*data[M][i+1])/(xi1-xi);
    
    double va = exp(-a*log(S0/K)-b*t_final)*ul;
    //cout<<va<<endl;
    Option_Pricing temp(S0,K,r,q,T,sigma);
    double ve = temp.E_BSM_Put_Price();
    return abs(ve-va);

}

double EFD::ErrorRMS()
{
    double sum = 0;
    int count = 0;
    for(int i = 0; i < N+1; i++)
    {
        double xi = x_left + i * (x_right-x_left)/N;
        double  va = exp(-a * xi -b * t_final)*data[M][i];
        Option_Pricing temp(K*exp(xi),K,r,q,T,sigma);
        double ve = temp.E_BSM_Put_Price();
        if(ve/S0>0.00001)
        {
            sum += (ve-va)*(ve-va)/ve/ve;
            count++;
        }
    }
    sum = sum/(double)count;
    sum = sqrt(sum);
    return sum;
}

double EFD::Delta()
{
    int i = (log(S0/K)-x_left)/(x_right-x_left)*N;
    double xi = x_left + i * (x_right-x_left)/N;
    double xi1 = x_left + (i +1)* (x_right-x_left)/N;
    double si = K * exp(xi);
    double si1 = K * exp(xi1);
    double vi = exp(-a * xi - b * t_final) * data[M][i];
    double vi1 = exp(-a * xi1 - b * t_final) * data[M][i+1];
    return (vi1-vi)/(si1-si);
}

double EFD::Gamma()
{
    int i = (log(S0/K)-x_left)/(x_right-x_left)*N;
    double x = x_left + (i-1)*(x_right-x_left)/N;
    double xi = x_left + i * (x_right-x_left)/N;
    double xi1 = x_left + (i +1)* (x_right-x_left)/N;
    double xi2 = x_left + (i + 2) * (x_right-x_left)/N;
    double s = K * exp(x);
    double si = K * exp(xi);
    double si1 = K * exp(xi1);
    double si2 = K * exp(xi2);
    Option_Pricing ti2(si2,K,r,q,T,sigma);
    Option_Pricing ti1(si1,K,r,q,T,sigma);
    Option_Pricing ti(si,K,r,q,T,sigma);
    Option_Pricing t(s,K,r,q,T,sigma);
    cout<<ti2.E_BSM_Put_Price()<<endl;
    cout<<ti1.E_BSM_Put_Price()<<endl;
    cout<<ti.E_BSM_Put_Price()<<endl;
    cout<<t.E_BSM_Put_Price()<<endl;
    
    double v = exp(-a * x - b * t_final) * data[M][i-1];
    double vi = exp(-a * xi - b * t_final) * data[M][i];
    double vi1 = exp(-a * xi1 - b * t_final) * data[M][i+1];
    double vi2 = exp(-a * xi2 - b * t_final) * data[M][i+2];
    //cout<<vi2<<endl<<vi1<<endl<<vi<<endl<<v<<endl;
    
    return((vi2-vi1)/(si2-si1)-(vi-v)/(si-s))/((si2+si1)/2-(si+s)/2);
    
}

double EFD::Theta()
{
    int i = (log(S0/K)-x_left)/(x_right-x_left)*N;
    double xi = x_left + i * (x_right-x_left)/N;
    double xi1 = x_left + (i +1)* (x_right-x_left)/N;
    double si = K * exp(xi);
    double si1 = K * exp(xi1);
    double t = 2 * (t_final-t_final/(double)M)/sigma/sigma;
    Option_Pricing ti1(si1,K,r,q,t,sigma);
    Option_Pricing ti(si,K,r,q,t,sigma);
    cout<<ti1.E_BSM_Put_Price()<<endl;
    cout<<ti.E_BSM_Put_Price()<<endl;
    
    double vit = exp(-a*xi-b*(t_final-t_final/M))*data[M-1][i];
    double vi1t= exp(-a*xi1-b*(t_final-t_final/M))*data[M-1][i+1];
    double vst = ((si1-S0)*vit + (S0-si)*vi1t)/(si1-si);
    double vi = exp(-a*xi-b*(t_final))*data[M][i];
    double vi1= exp(-a*xi1-b*(t_final))*data[M][i+1];
    double vso =((si1-S0)*vi + (S0-si)*vi1)/(si1-si);
    //cout<<vi1t<<endl<<vit<<endl;
    return (vso-vst)/(T/M);
}

void EFD :: ForwardEuler()
{
    
    for(int i = 1; i < M + 1; i++)
    {
        for( int j = 1; j< N ; j++)
        {
            data[i][j] = alpha * data[i - 1][j - 1] + ( 1 - 2 * alpha) * data[i-1][j] + alpha * data[i - 1][j + 1];
            
            
        }
    }
}

void EFD :: BackwardEuler(string method)
{
    Matrix<double> a;
    
    for(int i = 0; i < N-1;i++)
    {
        vector<double> temp;
        for(int j = 0; j < N-1; j++)
            temp.push_back(0);
        a.Data.push_back(temp);
    }
    for(int i = 0; i < N-1; i++)
    {
        a.Data[i][i] = 1.0 + 2.0 * alpha;
    }
    for(int i = 0; i < N-2; i++)
    {
        a.Data[i][i+1] = -alpha;
        a.Data[i+1][i] = -alpha;
    }
    
    for(int i = 1; i <= M; i++)
    {
        vector<double> b;
        for(int j = 1; j < N; j ++)
            b.push_back(data[i-1][j]);
        b[0] += alpha * data[i][0];
        b[N-2] += alpha * data[i][N];
        vector<double> result;
        
        if(method == "LU")
            result = a.linear_solve_Lu_no_pivoting(b);
        else if(method == "SOR")
        {
            vector<double> x0;
            for(int j = 1; j < N; j ++)
                x0.push_back(0);
            result = a.SORP(b,1.2, x0);
        }
        
        for(int j = 1; j < N; j++ )
        {
            data[i][j] = result[j-1];

        }
    }
    
}

void EFD :: CrankNicolson(string method)
{
    Matrix<double> a,b1;
    
    for(int i = 0; i < N-1;i++)
    {
        vector<double> temp;
        for(int j = 0; j < N-1; j++)
        {
            temp.push_back(0);
        }
        a.Data.push_back(temp);
        b1.Data.push_back(temp);
    }
    for(int i = 0; i < N-1; i++)
    {
        a.Data[i][i] = 1.0 + alpha;
        b1.Data[i][i] = 1.0 - alpha;
    }
    for(int i = 0; i < N-2; i++)
    {
        a.Data[i][i+1] = -alpha/2.0;
        a.Data[i+1][i] = -alpha/2.0;
        b1.Data[i][i+1] = alpha/2.0;
        b1.Data[i+1][i] = alpha/2.0;
    }
    
    for(int i = 1; i <= M; i++)
    {
        vector<double> b ;
        for(int j = 1; j < N; j ++)
            b.push_back(data[i-1][j]);
        b = b1*b;
        
        b[0] += alpha /2.0 * (data[i][0] + data[i-1][0]);
        b[N-2] += alpha /2.0 * (data[i][N] + data[i-1][N]);
        vector<double> result;
        
        if(method == "LU")
        {
            result = a.linear_solve_Lu_no_pivoting(b);
            for(int j = 1; j < N; j++ )
            {
                data[i][j] = result[j-1];
                           }
        }
        else if(method == "SOR")
        {
            vector<double> x0;
            for(int j = 1; j < N; j ++)
                x0.push_back(0);
            
            //修改tol
            double tol = 1e-8;
            
            vector<double> x_old = x0;
            vector<double> x_new = x0;
            double w = 1.2;
            
            while(true)
            {
                x_old = x_new;
                for(int j = 0; j < x0.size(); j++)
                {
                    double sum1 = 0, sum2 = 0;
                    for(int k = 0; k < j ; k++)
                    {
                        sum1 += a.Data[j][k] * x_new[k];
                    }
                    for(int k = j+1; k< x0.size(); k++)
                    {
                        sum2 += a.Data[j][k] * x_old[k];
                    }
                    x_new[j] = (1-w) * x_old[j] - w/a.Data[j][j] * (sum1 + sum2 ) + w * b[j]/a.Data[j][j];
                }
                
                double sum = 0;
                for(int j = 0; j < x_old.size(); j ++)
                {
                    sum += (x_new[j] - x_old[j]) * (x_new[j] - x_old[j]);
                }
                sum = sqrt(sum);
                if(sum < tol)
                    break;
            }
            for(int j = 1; j < N; j ++)
                data[i][j] = x_new[j-1];
        }
        
        
    }
}

double EFD::Value_approximate()
{
    int i = (log(S0/K)-x_left)/(x_right-x_left)*(double)N;
    double xi = x_left + i * (x_right-x_left)/N;
    double xi1 = x_left + (i+1) * (x_right-x_left)/N;
    double si = K * exp(xi);
    double si1 = K * exp(xi1);
    double vi = exp(-a*xi - b * t_final)* data[M][i];
    double vi1 = exp(-a*xi1 - b * t_final) * data[M][i+1];
    double va = ((si1-S0)*vi+(S0-si)*vi1)/(si1-si);
    return va;
}



