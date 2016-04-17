//
//  BSFD.hpp
//  Finite Different Method
//
//  Created by LiuYuxi on 12/11/15.
//  Copyright © 2015 LiuYuxi. All rights reserved.
//

#ifndef BSFD_hpp
#define BSFD_hpp

#include <iostream>
#include "FD.hpp"
#include "Option_pricing.hpp"

class EFD:public FD
{
public:
    double S0;
    double K;
    double r;
    double q;
    double T;
    double sigma;
    double a;
    double b;
    double x_compute;
    EFD();
    EFD(double, double, double, double, double, double, int, double);
    ~EFD();
    
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

#endif /* BSFD_hpp */
