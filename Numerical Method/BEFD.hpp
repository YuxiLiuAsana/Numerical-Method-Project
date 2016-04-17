//
//  BEFD.hpp
//  Finite Different Method
//
//  Created by LiuYuxi on 12/15/15.
//  Copyright © 2015 LiuYuxi. All rights reserved.
//

#ifndef BEFD_hpp
#define BEFD_hpp

#include <iostream>
#include "FD.hpp"
#include "Option_pricing.hpp"

class BEFD:public FD
{
public:
    double S0;
    double K;
    double r;
    double q;
    double T;
    double U;
    double D;
    double sigma;
    double a;
    double b;
    double x_compute;
    BEFD();
    BEFD(double, double, double, double, double, double, double, double, int, double);
    ~BEFD();
    
    double f(int i);//i represents the index
    double g_left(int j);
    double g_right(int j);
    double ErrorPointwise1();
    double ErrorPointwise2();
    double ErrorRMS();
    double Delta();
    double Gamma();
    double Theta();
    void ForwardEuler();
    void BackwardEuler(string method);
    void CrankNicolson(string method);
    double Value_approximate();
    
    
};


#endif /* BEFD_hpp */
