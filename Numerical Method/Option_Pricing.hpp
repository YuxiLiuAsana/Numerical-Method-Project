//
//  Option_Pricing.hpp
//  NLA
//
//  Created by LiuYuxi on 10/2/15.
//  Copyright © 2015 LiuYuxi. All rights reserved.
//

#ifndef Option_Pricing_hpp
#define Option_Pricing_hpp

#include <iostream>
#include<cmath>
#include<fstream>
#include <boost/math/distributions/normal.hpp>
using namespace std;
class Option_Pricing
{
public:
    double So;
    double K;
    double r;
    double q;
    double T;
    double sigma;
    double d1;
    double d2;
   
    
public:
    Option_Pricing();
    Option_Pricing(double, double, double, double, double, double);
    Option_Pricing(Option_Pricing&);
    ~Option_Pricing();
    
    Option_Pricing& operator= (Option_Pricing& );
    
    double N(double d);//cumulative normal distribution
    double n(double d);//percentage normal distribution
    
    /*
     The following functions are named in this way:
     E/A: European option or American Option
     
     BSM/BT/ABT/BBSM/BBSMRE:
     BSM: Black Scholes model
     BT: Binomial Tree model
     ABT: Average Biniomial Tree model
     BBSM: Binomial Black-Scholes
     BBSMRE: Binomial Black-Scholes model with Richardson Extrapolation
     TT : Trinomial Tree
     
     Call/Put: put or call option
     
     Price/Delta/Gamma/Theta/Rho/Vega: The corresponding objective
    
    */
    double A_ABT_Put_Price(int );
    double A_ABT_Put_Delta(int );
    double A_ABT_Put_Gamma(int );
    double A_ABT_Put_Theta(int );
    
    double A_BBSM_Put_Price(int );
    
    double A_BBSMRE_Put_Price(int );
    
    double A_BT_Put_Price(int );
    double A_BT_Put_Delta(int );
    double A_BT_Put_Gamma(int );
    double A_BT_Put_Theta(int );
    
    double A_BT_Call_Price(int );
    
    
    
    double A_TBSM_Put_Price(int );
    double A_TBSM_Put_Delta(int );
    double A_TBSM_Put_Gamma(int );
    double A_TBSM_Put_Theta(int );
    
    double A_TBSMRE_Put_Price(int );
    double A_TBSMRE_Put_Delta(int );
    double A_TBSMRE_Put_Gamma(int );
    double A_TBSMRE_Put_Theta(int );
    
    double A_TT_Put_Price(int );
    double A_TT_Put_Delta(int );
    double A_TT_Put_Gamma(int );
    double A_TT_Put_Theta(int );
  

    
    double E_ABT_Put_Price(int );
    double E_ABT_Put_Delta(int );
    double E_ABT_Put_Gamma(int );
    double E_ABT_Put_Theta(int );
    
    vector<double> E_BADT_Pricer(int );
    
    double E_BBSM_Put_Price(int );
    
    double E_BBSMRE_Put_Price(int );
    
    
    double E_BSM_Call_Price();
    double E_BSM_Put_Price();
    double E_BSM_Put_Delta();
    double E_BSM_Put_Gamma();
    double E_BSM_Put_Theta();
    
  
    double E_BT_Call_Price(int );
    
    double E_BT_Put_Price(int );
    double E_BT_Put_Delta(int );
    double E_BT_Put_Gamma(int );
    double E_BT_Put_Theta(int );
    
    vector<double> E_TADT_Pricer(int );

    double E_TT_Put_Price(int );
    double E_TT_Put_Delta(int );
    double E_TT_Put_Gamma(int );
    double E_TT_Put_Theta(int );
    
    double E_TBSM_Put_Price(int );
    double E_TBSM_Put_Delta(int );
    double E_TBSM_Put_Gamma(int );
    double E_TBSM_Put_Theta(int );
    
    double E_TBSMRE_Put_Price(int );
    double E_TBSMRE_Put_Delta(int );
    double E_TBSMRE_Put_Gamma(int );
    double E_TBSMRE_Put_Theta(int );
    
    
    
    
    
    
    
};



#endif /* Option_Pricing_hpp */
