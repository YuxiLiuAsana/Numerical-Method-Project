//
//  EFDD.hpp
//  Finite Different Method
//
//  Created by LiuYuxi on 12/13/15.
//  Copyright © 2015 LiuYuxi. All rights reserved.
//

#ifndef EFDD_hpp
#define EFDD_hpp

#include "FD.hpp"
class EFDD:public FD
{
public:
    double S0;
    double K;
    double r;
    double q;
    double T;
    double sigma;
    double T_div;
    int M1 ;
    int M2;
    int N1;
    int N2;
    double a;
    double b;
    double x_compute;
    double t_div;
    EFDD(double,double, double, double, double ,double, int, double, double);
    EFDD();
    ~EFDD();
    double f(int i);
    double g_left_old(int i);
    double g_right_old(int i);
    double g_left_new(int i);
    double g_right_new(int i);
    double alpha1;
    double alpha2;
    void ForwardEuler_old();
    void ForwardEuler_new();
    void BackwardEuler_old(string method);
    void BackwardEuler_new(string method);
    void CrankNicolson_old(string method);
    void CrankNicolson_new(string method);
    double u_value();
    double Option_value();
    double Delta();
    double Theta();
    double Gamma();

    
};

#endif /* EFDD_hpp */
