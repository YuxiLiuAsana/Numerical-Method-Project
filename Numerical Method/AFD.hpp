//
//  AFD.hpp
//  Finite Different Method
//
//  Created by LiuYuxi on 12/11/15.
//  Copyright © 2015 LiuYuxi. All rights reserved.
//

#ifndef AFD_hpp
#define AFD_hpp

#include <iostream>
#include "EFD.hpp"
class AFD:public EFD
{
public:
    //change the initial value
    double ve = 2.45;
    AFD();
    AFD(double S0, double K, double r, double q, double T,double sigma, int M, double alpha);
    ~AFD();
    double f(int i);//i represents the index
    double g_left(int j);
    double g_right(int j);
    double EarlyExercisePremium(int i, int j);
    void ForwardEuler();
    void BackwardEuler(string method);
    void CrankNicolson(string method);
    double ErrorPointwise1();
    double ErrorPointwise2();
    double VarianceReductionValue(string Method);
    double ErrorPointwiseVR(string Method);
    vector<vector<double>> EarlyExerciseDomain();
    double Value_approximate();
};

#endif /* AFD_hpp */
