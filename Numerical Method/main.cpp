//
//  main.cpp
//  Finite Different Method
//
//  Created by LiuYuxi on 12/8/15.
//  Copyright © 2015 LiuYuxi. All rights reserved.
//

#include <iostream>

#include "AFD.hpp"
#include "EFDD.hpp"
#include "BEFD.hpp"
#include <iomanip>
template<typename T>
vector<vector<T>> Read_Matrix(char *s)
{
    vector<vector<T>> ret;
    vector<double> p;
    ifstream f1(s);
    string value,value1;
    while ( f1.good())
    {
        ofstream f2("b.csv");
        if(!f2)
        {
            ret.clear();
            return ret;
        }
        
        p.clear();
        f2.clear();
        getline (f1, value);
        f2<<value;
        f2.close();
        ifstream f3("b.csv");
        while(f3.good())
        {
            
            getline(f3,value1,',');
            char s[100] ;
            strcpy(s,value1.c_str());
            p.push_back(atof(s));
            
        }
        ret.push_back(p);
    }
    ret.pop_back();
    return ret;
}

template<typename T>
vector<T> Read_vector(char *s)
{
    vector<T> ret;
    ifstream f1(s);
    string value,value1;
    ret.clear();
    while ( f1.good())
    {
        ofstream f2("b.csv");
        if(!f2)
        {
            ret.clear();
            return ret;
        }
        
        
        f2.clear();
        getline (f1, value);
        f2<<value;
        f2.close();
        ifstream f3("b.csv");
        while(f3.good())
        {
            
            getline(f3,value1,',');
            char s[100] ;
            strcpy(s,value1.c_str());
            ret.push_back(atof(s));
            
        }
        
    }
    ret.pop_back();
    return ret;
}

template<typename T>
void Write_Matrix(Matrix<T>& source, char *s)
{
    ofstream f1(s);
    if(!f1) return;
    f1.clear();
    for(int i = 0; i < source.Data.size(); i++)
    {
        f1<<setiosflags(ios::scientific)<<setprecision(20)<<source.Data[i][0];
        for(int j = 1; j< source.Data[i].size();j++)
            f1<<','<<setiosflags(ios::scientific)<<setprecision(20)<<source.Data[i][j];
        f1<<endl;
    }
    f1.close();
}

template<typename T>
void Write_vector(vector<T>& source, bool col,char *s)
{
    ofstream f1(s);
    if(!f1) return;
    f1.clear();
    f1<<setiosflags(ios::scientific)<<setprecision(11)<<source[0];
    for(int i = 1; i < source.size(); i++)
    {
        if(col == true)
            f1<<endl;
        else 
            f1<<',';
        f1<<setiosflags(ios::scientific)<<setprecision(11)<<source[i];
    }
    f1<<endl;
    f1.close();
}

int main()
{
    
   /* //**********************************problem 2 **********************************************
    //改这个部分 初值
    //*****************************************************************
    /*double S0 = 27;
    double K = 30;
    double r = 0.03;
    double q = 0.01;
    double T = 5.0/12.0;
    double sigma = 0.35;*/
    //*****************************************************************
    //然后去改European 和 American 的boundary condition, European 要根据Call 还是 put 来修改
    
    cout<<setprecision(12);
    //int M = 1;
    //fixed computational domain
    //别忘记修改 alpha
    /*for(int i =0 ;i < 4; i++)
    {
        M = M * 4;
        AFD test(S0,K,r,q,T,sigma,M,0.5);
        cout<< test.alpha<<","<< test.N<<","<<test.x_left<<","<<test.x_right<<","<<test.x_compute<<","<<test.t_final<<","<<test.t_final/(double)test.M<<","<<(test.x_right-test.x_left)/(double)test.N<<endl;
    }*/
    
    //Solution and Variance Reduction Method
    /*M = 1;
    for(int i = 0; i < 4; i++)
    {
        M = M * 4;
        AFD test(S0,K,r,q,T,sigma, M ,4);//修改alpha
        //注意他要我用什么方法，修改“CS”，修改tol 和 w in SOR 或者 Projected SOR， in both American and European cpp
        test.VarianceReductionValue("CS");
    }*/
    
    //alpha_temp = 0.5, M = 4, N = 12
    /*AFD test(S0,K,r,q,T,sigma,4, 0.5);//注意修改M
    test.ForwardEuler();//注意修改方法，w，tol， 以及是否是Project SOR
    cout<<test;*///  看他给M的顺序，是否要反过来
    
    /*EFD test(S0,K,r,q,T,sigma, 4,0.5);
    test.ForwardEuler();
    cout<<test.Delta()<<endl;
    cout<<test.Gamma()<<endl;
    cout<<test.Theta()<<endl;*/
    
    /*int M = 4;
    for(int i = 0; i < 3; i++)
    {
        M = M * 4;
        EFD test(S0,K,r,q,T,sigma, M,0.5);
        test.ForwardEuler();
        cout<<test.Value_approximate()<<","<<test.Delta()<<","<<test.Gamma()<<","<<test.Theta()<<endl;
        
    }
    Option_Pricing test(S0,K,r,q,T,sigma);
    cout<<test.E_BSM_Put_Price()<<","<<test.E_BSM_Put_Delta()<<","<<test.E_BSM_Put_Gamma()<<","<<test.E_BSM_Put_Theta()<<endl;*/
    
    
    //change of domain
    /*M = 1;
    for(int i = 0 ; i < 4; i++)
    {
        M = M * 4;
        AFD test(S0,K,r, q,T,sigma,M,0.5);//修改domain在EFD.cpp
        test.VarianceReductionValue("F");
    }*/
     
     //implied volatility
    /*
     double sigma0 = 0.1;//修改sigma
     double sigma1 = 0.4;//修改sigma
     //cout<<sigma0<<","<<sigma1<<",";
     for(int i = 2; i < 11; i++)
     {
     AFD test0(S0, K, r,q, T, sigma0, 64,0.4);//注意题里面alpha和M的值
     AFD test1(S0, K ,r,q,T, sigma1, 64,0.4);
     test1.ForwardEuler();//注意他让用的方法
     test0.ForwardEuler();
     double temp1 = test1.Value_approximate()-test1.ve;
         //cout<<temp1<<",";
     double temp0 = test0.Value_approximate()-test0.ve;
         //cout<<temp0<<",";
     double sigma_new = sigma1 - (temp1 * (sigma1 - sigma0)/(temp1 - temp0));
     sigma0 = sigma1;
     sigma1 = sigma_new;
         
     cout<<sigma1<<endl;
     
     }*/
    
    /*AFD test(S0,K,r,q,T,sigma,4,0.5);
    test.ForwardEuler();
    //Option_Pricing s(S0,K,r,q,T,sigma);
    test.Theta();
    //cout<<s.E_BSM_Put_Price()<<endl;
    //cout<<s.E_BSM_Put_Delta()<<endl;
    //cout<<s.E_BSM_Put_Gamma()<<endl;*/
    //cout<<s.E_BSM_Put_Theta()<<endl;
    
    
    
    //Barrier Option
    
    
    
    
    
    //problem 3 ***********************************************************************
    
    //修改这部分的值
    //*************************************************************************************************************
    /*double S0 = 52;
    double K = 55;
    double D = 45;
    double U = 62;
    double r = 0.02;
    double q = 0.005;
    double T = 5.0/12.0;
    double sigma = 0.25;*/
    //************************************************************************************************************
     
     
    
     
     
     // 修改边界条件在BEFD.cpp中
     
     //domain
    /*int M = 1;
    for(int i = 0; i < 4; i++)
    {
        M = M * 4;
        BEFD test(S0, K ,D,U,r,q,T,sigma,M, 0.5);
        cout<<test.alpha<<","<<test.N<<","<<test.x_left<<","<<test.x_right<<","<<test.x_compute<<","<<test.t_final<<","<<test.t_final/(double) test.M<<","<<(test.x_right - test.x_left)/(double) test.N<<endl;
    }*/
    
    //solution
     /*int M = 1;
    for(int i = 0; i < 4; i++)
    {
        M = M * 4;
        BEFD test(S0, K, D, U, r, q, T, sigma, M ,4);
        test.CrankNicolson("LU");//注意看他的方法是什么
        test.ErrorPointwise1();
        cout<<test.Delta()<<","<<test.Gamma()<<","<<test.Theta()<<endl;
        
    }*/
    
    //int M = 256;
    /*BEFD test(S0,K,D,U,r,q,T,sigma, M ,0.5);
    test.ForwardEuler();//注意方法
    cout<<test;*/
    /*for(int i = 0; i < test.N + 1; i ++)
    {
        double x = test.x_left + (test.x_right - test.x_left)/(double)test.N * i;
        double t = test.T*test.sigma * test.sigma * 0.5;
        cout<< exp(-test.a * x - test.b * t)* test.data[M][i]<<",";
    }
    cout<<endl;*/
    
    /*BEFD test(S0,K,D,U,r,q,T,sigma, M ,0.5);
    test.ForwardEuler();
    test.Theta();*/
    
    
    
    
    
    
    
    
    
    
    
    
    //*************************problem one************************
    /*double S0 = 47;
    double K = 50;
    double sigma = 0.25;
    double q = 0.01;
    double r = 0.03;
    double T = 7.0/12.0;
    
    Option_Pricing test(S0,K,r,q,T,sigma);*/
    //test.A_BT_Put_Price(200);
    //cout<<test.A_BT_Put_Delta(200)<<endl;
    //cout<<test.A_BT_Put_Gamma(200)<<endl;
    //cout<<test.A_BT_Put_Theta(200)<<endl;
    /*int step = 5;
    while(true)
    {
        step = step * 2;
        double a = test.A_BT_Put_Price(step);
        double b = test.A_BT_Put_Price(2 * step);
        if(abs(b-a)<1e-4)
            break;
    }
    cout<<step * 2;*/
    //test.A_BT_Put_Price(5120);
    //cout<<test.A_BT_Put_Delta(5120)<<endl;
    //cout<<test.A_BT_Put_Gamma(5120)<<endl;
    //cout<<test.A_BT_Put_Theta(5120)<<endl;
    //test.E_BT_Put_Price(5120);
    //double dt = T/5120.0;
    //double sup = S0 * exp(2* sigma * sqrt(dt));
    //double sdown = S0 * exp(-2*sigma * sqrt(dt));
    //cout<<test.E_BSM_Put_Price();
    //Option_Pricing testu(sdown,K,r,q,T-2*dt,sigma);
    //cout<<testu.E_BSM_Put_Price();
    //cout<<test.E_BT_Put_Delta(5120)<<endl;
    //cout<<test.E_BT_Put_Gamma(5120)<<endl;
    //cout<<test.E_BT_Put_Theta(5120)<<endl;
    
    
    
    //**************************checking***************
    double So = 40;
    double K = 41;
    double r = 0.03;
    double q = 0.01;
    double T = 1;
    double sigma = 0.3;
    Option_Pricing test(So,K,r,q,T,sigma);
    cout<<test.A_BT_Call_Price(1);
    
    
    
    
}
