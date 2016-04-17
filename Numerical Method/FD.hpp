//
//  FD.hpp
//  Finite Different Method
//
//  Created by LiuYuxi on 12/8/15.
//  Copyright © 2015 LiuYuxi. All rights reserved.
//

#ifndef FD_hpp
#define FD_hpp

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "Matrix.h"
#include "Matrix.cpp"
using namespace std;
class FD
{
public:
    double x_left;// index 0
    double x_right; // index N
    double t_final; // index M
    int M;
    int N;
    double alpha;
    vector<vector<double>> data;
    FD();
    FD( double, double, double, int, int, double);
    ~FD();
    
    double f( int );
    double g_left( int );
    double g_right( int );
    double u_exact( int, int );
    void ForwardEuler();
    void BackwardEuler(string method);
    void CrankNicolson(string method);
    double MaxPointwiseError();
    double ErrorRMS();
    friend ostream& operator<<(ostream &out, FD& source);
    
};

#endif /* FD_hpp */
