//
//  AFD.cpp
//  Finite Different Method
//
//  Created by LiuYuxi on 12/11/15.
//  Copyright © 2015 LiuYuxi. All rights reserved.
//

#include "AFD.hpp"
AFD::AFD(double _S0, double _K, double _r, double _q, double _T, double _sigma, int _M, double _alpha ):EFD(_S0,_K,_r,_q,_T,_sigma,_M,_alpha)
{
    for(int j = 0; j < M + 1; j++)
    {
        data[j][0] = g_left(j);
        data[j][N] = g_right(j);
    }
    for(int j = 0; j < N + 1; j ++)
    {
        data[0][j] = f(j);
    }
}

AFD::AFD():EFD()
{
    
}

AFD::~AFD()
{
    
}


























//改这部分 边界条件

//****************************************************************************************************
double AFD::f(int i)
{
    //call
    double x = x_left + i * (x_right - x_left) / N;
    //return K * exp(a * x ) * (exp(x)-1) * (1 < exp(x));
    
    //put
    return K * exp(a * x) * (1- exp(x)) * (1 > exp(x));
}


// 边界条件修改
double AFD::g_left(int j)
{
    double t = t_final / (double)M * (double)j;
    //call
    //return 0;
    
    //put
    return K * exp(a * x_left + b * t) * (1-exp(x_left));
}


//边界条件修改
double AFD::g_right(int j)
{
    double t = t_final / (double)M * (double)j;
    //call
    //return K * exp(a * x_right + b * t)*(exp(x_right)-1);
    //put
    return 0;
    
}
// 提前执行修改
double AFD::EarlyExercisePremium(int i, int j)
{
    //call
    //double m = -K * exp(a * (x_left + (x_right - x_left)*j/(double)N) + b * t_final/(double)M * i) * (1 - exp(x_left+(x_right - x_left)*j/(double)N)) * (1 < exp(x_left+(x_right - x_left)*j/(double)N)) ;
    //put
    double m = K * exp(a * (x_left + (x_right - x_left)*j/(double)N) + b * t_final/(double)M * i) * (1 - exp(x_left+(x_right - x_left)*j/(double)N)) * (1 > exp(x_left+(x_right - x_left)*j/(double)N)) ;
    return m;
}


//**************************************************************************************************************


























void AFD :: ForwardEuler()
{
    
    for(int i = 1; i < M + 1; i++)
    {
        for( int j = 1; j< N ; j++)
        {
            data[i][j] = alpha * data[i - 1][j - 1] + ( 1 - 2 * alpha) * data[i-1][j] + alpha * data[i - 1][j + 1];
            
            
            if(data[i][j] < EarlyExercisePremium(i,j))
                data[i][j] = EarlyExercisePremium(i,j);
        }
    }
}

void AFD :: BackwardEuler(string method)
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
                x0.push_back(EarlyExercisePremium(i, j));
            result = a.SORP(b,1.2, x0);
        }
        
        for(int j = 1; j < N; j++ )
        {
            data[i][j] = result[j-1];
            if(data[i][j] < EarlyExercisePremium(i, j))
                data[i][j] = EarlyExercisePremium(i, j);
        }
    }
    
}

void AFD :: CrankNicolson(string method)
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
                if(data[i][j] < EarlyExercisePremium(i, j))
                    data[i][j] = EarlyExercisePremium(i, j);
            }
        }
        else if(method == "SOR")
        {
            vector<double> x0;
            for(int j = 1; j < N; j ++)
                x0.push_back(EarlyExercisePremium(i, j));
            
            //修改tol
            //***********************************************************************************************
            double tol = 1e-8;
            //***********************************************************************************************
            
            
            
            
            
            
            vector<double> x_old = x0;
            vector<double> x_new = x0;
            
            
            
            
            
            
            
            
            
            
            // 修改w
            //**********************************************************************************************
            double w = 1.2;
            //**********************************************************************************************
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
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
                    x_new[j] = max((1-w) * x_old[j] - w/a.Data[j][j] * (sum1 + sum2 ) + w * b[j]/a.Data[j][j], EarlyExercisePremium(i, j+1));
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

double AFD::ErrorPointwise1()
{
    int i = (log(S0/K)-x_left)/(x_right-x_left)*(double)N;
    double xi = x_left + i * (x_right-x_left)/N;
    double xi1 = x_left + (i+1) * (x_right-x_left)/N;
    double si = K * exp(xi);
    double si1 = K * exp(xi1);
    double vi = exp(-a*xi - b * t_final)* data[M][i];
    double vi1 = exp(-a*xi1 - b * t_final) * data[M][i+1];
    double va = ((si1-S0)*vi+(S0-si)*vi1)/(si1-si);
    return abs(ve-va);
}

double AFD::ErrorPointwise2()
{
    int i = (log(S0/K)-x_left)/(x_right-x_left)*N;

    double xi = x_left + i * (x_right-x_left)/N;

    double xi1 = x_left + (i+1) * (x_right-x_left)/N;
    
    double ul =((xi1-log(S0/K))*data[M][i]+(log(S0/K)-xi)*data[M][i+1])/(xi1-xi);
    
    double va = exp(-a*log(S0/K)-b*t_final)*ul;
 
    
    return abs(ve-va);

    
}

double AFD::VarianceReductionValue(string Method)
{
       if(Method == "F")
    {
        ForwardEuler();
    }
    else if(Method == "BL")
    {
        BackwardEuler("LU");
    }
    else if(Method == "BS")
    {
        BackwardEuler("SOR");
    }
    else if(Method == "CL")
    {
        CrankNicolson("LU");
    }
    else if(Method == "CS")
    {
        CrankNicolson("SOR");
    }
    int i = (log(S0/K)-x_left)/(x_right-x_left)*N;
    cout<<data[M][i]<<","<<data[M][i+1]<<",";
    
    double xi = x_left + i * (x_right-x_left)/(double)N;
    
    double xi1 = x_left + (i+1) * (x_right-x_left)/(double)N;
    
    double si = K * exp(xi);
    
    double si1 = K * exp(xi1);
    
    double vi = exp(-a*xi - b * t_final)* data[M][i];
    
    double vi1 = exp(-a*xi1 - b * t_final) * data[M][i+1];
    // cout<<data[M][i]<<","<<data[M][i+1]<<endl;
    double v_american_approximate = ((si1-S0)*vi+(S0-si)*vi1)/(si1-si);
    
    EFD e(S0,K,r,q,T,sigma, M, alpha);
    if(Method == "F")
    {
        e.ForwardEuler();
    }
    else if(Method == "BL")
    {
        e.BackwardEuler("LU");
    }
    else if(Method == "BS")
    {
        e.BackwardEuler("SOR");
    }
    else if(Method == "CL")
    {
        e.CrankNicolson("LU");
    }
    else if(Method == "CS")
    {
        e.CrankNicolson("SOR");
    }

    xi = e.x_left + i * (e.x_right-e.x_left)/(double)e.N;
    xi1 = e.x_left + (i+1) * (e.x_right-e.x_left)/(double)e.N;
        
    si = e.K * exp(xi);
        
    si1 = e.K * exp(xi1);
        
    vi = exp(-e.a*xi - e.b * e.t_final)* e.data[M][i];
        
    vi1 = exp(-e.a*xi1 - e.b * e.t_final) * e.data[M][i+1];
    
    double v_european_approximate = ((si1-e.S0)*vi+(e.S0-si)*vi1)/(si1-si);
    
    Option_Pricing temp(S0,K,r,q,T,sigma);
    
    //注意修改 call 和 put
    double v_european_exact = temp.E_BSM_Put_Price();
    
    cout<<v_american_approximate<<","<<v_european_approximate<<","<<v_american_approximate-v_european_approximate + v_european_exact<<endl;
    return v_american_approximate-v_european_approximate + v_european_exact;
    
}

double AFD::ErrorPointwiseVR(string Method)
{
    return abs(VarianceReductionValue(Method)-ve);
}

//call 和 put 不一样
vector<vector<double>> AFD::EarlyExerciseDomain()
{
    vector<double> as;
    vector<double> bs;
    for(int i = 0 ; i < M+1; i++)
    {
        for(int j = 0; j < N+1; j++)
        {
            if(data[i][j] != EarlyExercisePremium(i, j))
            {
                as.push_back(T - 2 * i * t_final / (double)M / sigma / sigma);
                
                bs.push_back(K *(exp(x_left + (j-1)/(double)N * (x_right-x_left))+exp(x_left + j/(double)N * (x_right-x_left)))/2);
                break;
            }
            
        }
    }
    vector<vector<double>> s;
    s.push_back(as);
    s.push_back(bs);
    return s;
}

double AFD::Value_approximate()
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


