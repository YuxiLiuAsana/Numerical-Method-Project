//
//  Option_Pricing.cpp
//  NLA
//
//  Created by LiuYuxi on 10/2/15.
//  Copyright © 2015 LiuYuxi. All rights reserved.
//

//so,k,r,q,t,sigma

#include "Option_Pricing.hpp"
#define PI 3.14159265358979323846

Option_Pricing::Option_Pricing():K(40), So(41) , q(0.01), sigma(0.3), r(0.03), T(1)
{
    d1 = 1/sigma/sqrt(T)*(log(So/K) + (r - q + sigma * sigma / 2) * T);
    d2 = d1 - sigma * sqrt(T);
}

Option_Pricing::Option_Pricing(double So1, double K1, double r1, double q1, double T1, double sigma1):K(K1), So(So1) , q(q1), sigma(sigma1), r(r1), T(T1)
{
    d1 = 1/sigma/sqrt(T)*(log(So/K) + (r - q + sigma * sigma / 2) * T);
    d2 = d1 - sigma * sqrt(T);
}

Option_Pricing::Option_Pricing(Option_Pricing& source):K(source.K), So(source.So) , q(source.q), sigma(source.sigma), r(source.r), T(source.T)
{
    d1 = 1/sigma/sqrt(T)*(log(So/K) + (r - q + sigma * sigma / 2) * T);
    d2 = d1 - sigma * sqrt(T);
}

Option_Pricing::~Option_Pricing()
{

}

Option_Pricing& Option_Pricing::operator= (Option_Pricing& source)
{
    K = source.K;
    So = source.So;
    r = source.r;
    q = source.q;
    T = source.T;
    sigma = source.sigma;
    d1 = source.d1;
    d2 = source.d2;
    return *this;
}

double Option_Pricing::N(double d)//cumulative normal distribution
{
    boost::math::normal_distribution<double> myNormal;
    return cdf(myNormal,d);
}

double Option_Pricing::n(double d)//percentage normal distribution
{
    boost::math::normal_distribution<double> myNormal;
    return pdf(myNormal,d);
}

double Option_Pricing::E_BSM_Put_Price()
{
    return N(-d2) * K * exp(-r * T) - N(-d1) * So * exp( -q * T);
}


double Option_Pricing::E_BSM_Put_Delta()
{
    return exp(-q * T) * (N(d1) - 1);
}


double Option_Pricing::E_BSM_Put_Gamma()
{
    return n(d1)*exp(-q * T)/So/sigma/sqrt(T);
}

double Option_Pricing::E_BSM_Put_Theta()
{
    return -So * sigma * exp(-q * T- d1 * d1 / 2)/(2*sqrt(2 * T * PI)) - q * So * exp( - q * T) * N( -d1 ) + r * K * exp( - r * T) * N(-d2);
}

double Option_Pricing::E_BT_Put_Price(int step)
{
    double price[step+1][2];
    double temp = So * exp(sigma * step * sqrt(T/step));
    for(int i = 0; i < step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(-2 * sigma * sqrt(T/step));
    }
    for(int i = step ; i >= 1; i--)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown;
            pup = (exp((r-q) * T/step) - exp(- sigma * sqrt(T/step)))/(exp(sigma * sqrt(T/step)) - exp(-sigma * sqrt(T/step)));
            pdown =1-pup;
            price[j][1] = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pdown);
            price[j][0] = price[j][0]/exp(sigma * sqrt(T/step));
        }
        //code fixed part***********************
        if(i == 3)
            cout<<price[0][1]<<","<<price[1][1]<<","<<price[2][1]<<endl;
        if(i == 2)
            cout<<price[0][1]<<","<<price[1][1]<<endl;
        if(i == 1)
            cout<<price[0][1]<<endl;
    }
    return price[0][1];
    
}

double Option_Pricing::E_BT_Put_Delta(int step)
{
    double price[step+1][2];
    double temp = So * exp(sigma * step * sqrt(T/step));
    for(int i = 0; i < step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(-2 * sigma * sqrt(T/step));
    }
    for(int i = step ; i >= 2; i--)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown;
            pup = (exp((r-q) * T/step) - exp(- sigma * sqrt(T/step)))/(exp(sigma * sqrt(T/step)) - exp(-sigma * sqrt(T/step)));
            pdown =1-pup;
            price[j][1] = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pdown);
            price[j][0] = price[j][0]/exp(sigma * sqrt(T/step));
        }
    }
    return (price[0][1]-price[1][1])/(price[0][0]-price[1][0]);

}

double Option_Pricing::E_BT_Put_Gamma(int step)
{
    double price[step+1][2];
    double temp = So * exp(sigma * step * sqrt(T/step));
    for(int i = 0; i < step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(-2 * sigma * sqrt(T/step));
    }
    for(int i = step ; i >= 3; i--)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown;
            pup = (exp((r-q) * T/step) - exp(- sigma * sqrt(T/step)))/(exp(sigma * sqrt(T/step)) - exp(-sigma * sqrt(T/step)));
            pdown =1-pup;
            price[j][1] = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pdown);
            price[j][0] = price[j][0]/exp(sigma * sqrt(T/step));
        }
    }
    return ((price[0][1]-price[1][1])/(price[0][0]-price[1][0])-(price[1][1]-price[2][1])/(price[1][0]-price[2][0]))/(price[0][0] - price[2][0] )*2;
}

double Option_Pricing::E_BT_Put_Theta(int step)
{
    double price[step+1][2];
    double temp = So * exp(sigma * step * sqrt(T/step));
    for(int i = 0; i < step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(-2 * sigma * sqrt(T/step));
    }
    double temp1=0.0;
    for(int i = step ; i >= 1; i--)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown;
            pup = (exp((r-q) * T/step) - exp(- sigma * sqrt(T/step)))/(exp(sigma * sqrt(T/step)) - exp(-sigma * sqrt(T/step)));
            pdown =1-pup;
            price[j][1] = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pdown);
            price[j][0] = price[j][0]/exp(sigma * sqrt(T/step));
        }
        if(i == 3)
            temp1 = price[1][1];
    }
    return (temp1 - price[0][1]) * step / T/2;
}


double Option_Pricing::E_ABT_Put_Price(int step)
{
    return (E_BT_Put_Price(step) + E_BT_Put_Price(step + 1))/2;
}

double Option_Pricing::E_ABT_Put_Delta(int step)
{
    return (E_BT_Put_Delta(step) + E_BT_Put_Delta(step + 1))/2;
}

double Option_Pricing::E_ABT_Put_Gamma(int  step)
{
    return (E_BT_Put_Gamma(step) + E_BT_Put_Gamma(step + 1))/2;
}

double Option_Pricing::E_ABT_Put_Theta(int  step)
{
    return (E_BT_Put_Theta(step) + E_BT_Put_Theta(step + 1))/2;
}

double Option_Pricing::A_BT_Put_Price(int step)
{
    double price[step+1][2];
    double temp = So * exp(sigma * step * sqrt(T/step));
    for(int i = 0; i < step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(-2 * sigma * sqrt(T/step));
    }
    for(int i = step ; i >= 1; i--)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown;
            pup = (exp((r-q) * T/step) - exp(- sigma * sqrt(T/step)))/(exp(sigma * sqrt(T/step)) - exp(-sigma * sqrt(T/step)));
            pdown =1-pup;
            price[j][0] = price[j][0]/exp(sigma * sqrt(T/step));
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pdown);
            if(K-price[j][0] < st)
            {
                price[j][1] = st;
            }
            else
            {
                price[j][1] = K - price[j][0];
            }
            
        }
        
        //code fixed part***********************
        if(i == 3)
            cout<<price[0][1]<<","<<price[1][1]<<","<<price[2][1]<<endl;
        if(i == 2)
            cout<<price[0][1]<<","<<price[1][1]<<endl;
        if(i == 1)
            cout<<price[0][1]<<endl;
    }
    return price[0][1];
}

double Option_Pricing::A_BT_Put_Delta(int step)
{
    double price[step+1][2];
    double temp = So * exp(sigma * step * sqrt(T/step));
    for(int i = 0; i < step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(-2 * sigma * sqrt(T/step));
    }
    for(int i = step ; i >= 2; i--)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown;
            pup = (exp((r-q) * T/step) - exp(- sigma * sqrt(T/step)))/(exp(sigma * sqrt(T/step)) - exp(-sigma * sqrt(T/step)));
            pdown =1-pup;
            price[j][0] = price[j][0]/exp(sigma * sqrt(T/step));
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pdown);
            if(K-price[j][0] < st)
            {
                price[j][1] = st;
            }
            else
            {
                price[j][1] = K - price[j][0];
            }
            
        }
    }
    return (price[0][1]-price[1][1])/(price[0][0]-price[1][0]);
}
double Option_Pricing::A_BT_Put_Gamma(int step)
{
    double price[step+1][2];
    double temp = So * exp(sigma * step * sqrt(T/step));
    for(int i = 0; i < step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(-2 * sigma * sqrt(T/step));
    }
    for(int i = step ; i >= 3; i--)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown;
            pup = (exp((r-q) * T/step) - exp(- sigma * sqrt(T/step)))/(exp(sigma * sqrt(T/step)) - exp(-sigma * sqrt(T/step)));
            pdown =1-pup;
            price[j][0] = price[j][0]/exp(sigma * sqrt(T/step));
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pdown);
            if(K-price[j][0] < st)
            {
                price[j][1] = st;
            }
            else
            {
                price[j][1] = K - price[j][0];
            }
            
        }
    }
    return ((price[0][1]-price[1][1])/(price[0][0]-price[1][0])-(price[1][1]-price[2][1])/(price[1][0]-price[2][0]))/(price[0][0] - price[2][0] )* 2;
}
double Option_Pricing::A_BT_Put_Theta(int step)
{
    
    double price[step+1][2];
    double temp = So * exp(sigma * step * sqrt(T/step));
    for(int i = 0; i < step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(-2 * sigma * sqrt(T/step));
    }
    double temp1=0.0;
    for(int i = step ; i >= 1; i--)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown;
            pup = (exp((r-q) * T/step) - exp(- sigma * sqrt(T/step)))/(exp(sigma * sqrt(T/step)) - exp(-sigma * sqrt(T/step)));
            pdown =1-pup;
            price[j][0] = price[j][0]/exp(sigma * sqrt(T/step));
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pdown);
            if(K-price[j][0] < st)
            {
                price[j][1] = st;
            }
            else
            {
                price[j][1] = K - price[j][0];
            }
            
        }
        if(i == 3)
            temp1 = price[1][1];
    }
    return (temp1 - price[0][1]) * step / T / 2;
}

double Option_Pricing::A_ABT_Put_Price(int step)
{
    return (A_BT_Put_Price(step) + A_BT_Put_Price(step + 1))/2;
}

double Option_Pricing::A_ABT_Put_Delta(int step)
{
    return (A_BT_Put_Delta(step) + A_BT_Put_Delta(step + 1))/2;
}
double Option_Pricing::A_ABT_Put_Gamma(int step)
{
    return (A_BT_Put_Gamma(step) + A_BT_Put_Gamma(step + 1))/2;
}

double Option_Pricing::A_ABT_Put_Theta(int step)
{
    return (A_BT_Put_Theta(step) + A_BT_Put_Theta(step + 1))/2;
}

double Option_Pricing::A_TT_Put_Price(int step)
{
    double price[2*step+1][2];
    double temp = So * exp(sigma * step * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    for(int i = 2 * step - 1 ; i >= 1; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            if(K-price[j][0] < st)
            {
                price[j][1] = st;
            }
            else
            {
                price[j][1] = K - price[j][0];
            }
            
        }
    }
    return price[0][1];
}
double Option_Pricing::A_TT_Put_Delta(int step)
{
    double price[2*step+1][2];
    double temp = So * exp(sigma * step * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    for(int i = 2 * step - 1 ; i >= 3; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            if(K-price[j][0] < st)
            {
                price[j][1] = st;
            }
            else
            {
                price[j][1] = K - price[j][0];
            }
            
        }
    }
    return (price[0][1]-price[2][1])/(price[0][0]-price[2][0]);
}
double Option_Pricing::A_TT_Put_Gamma(int step)
{
    double price[2*step+1][2];
    double temp = So * exp(sigma * step * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    for(int i = 2 * step - 1 ; i >= 5; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            if(K-price[j][0] < st)
            {
                price[j][1] = st;
            }
            else
            {
                price[j][1] = K - price[j][0];
            }
            
        }
    }
    return ((price[0][1]-price[2][1])/(price[0][0]-price[2][0])-(price[2][1]-price[4][1])/(price[2][0]-price[4][0]))/(price[1][0]-price[3][0]);
}
double Option_Pricing::A_TT_Put_Theta(int step)
{
    double price[2*step+1][2];
    double temp = So * exp(sigma * step * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    double temp1 = 0.0;
    for(int i = 2 * step - 1 ; i >= 1; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            if(K-price[j][0] < st)
            {
                price[j][1] = st;
            }
            else
            {
                price[j][1] = K - price[j][0];
            }
            
        }
        if(i == 3)
            temp1 = price[1][1];
        
    }
    return (temp1 - price[0][1])/T * step;
}

double Option_Pricing::E_TT_Put_Price(int step)
{
    double price[2*step+1][2];
    double temp = So * exp(sigma * step * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    for(int i = 2 * step - 1 ; i >= 1; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            price[j][1] = st;
            
            
        }
    }
    return price[0][1];
}

double Option_Pricing::E_TT_Put_Delta(int step)
{
    double price[2*step+1][2];
    double temp = So * exp(sigma * step * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    for(int i = 2 * step - 1 ; i >= 3; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            price[j][1] = st;
           
            
        }
    }
    return (price[0][1]-price[2][1])/(price[0][0]-price[2][0]);
}

double Option_Pricing::E_TT_Put_Gamma(int step)
{
    double price[2*step+1][2];
    double temp = So * exp(sigma * step * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    for(int i = 2 * step - 1 ; i >= 5; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            price[j][1] = st;
            
            
        }
    }
    return ((price[0][1]-price[2][1])/(price[0][0]-price[2][0])-(price[2][1]-price[4][1])/(price[2][0]-price[4][0]))/(price[1][0]-price[3][0]);

}

double Option_Pricing::E_TT_Put_Theta(int step)
{
    double price[2*step+1][2];
    double temp = So * exp(sigma * step * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = 0;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    double temp1 = 0.0;
    for(int i = 2 * step - 1 ; i >= 1; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            price[j][1] = st;
                        
        }
        if(i == 3)
            temp1 = price[1][1];
        
    }
    return (temp1 - price[0][1])/T * step;
}

double Option_Pricing::A_TBSM_Put_Price(int step)
{
    double price[2*step-1][2];
    double temp = So * exp(sigma * (step-1) * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step - 1; i++)
    {
        price[i][0] = temp;
        Option_Pricing flag(price[i][0],K,r,q,T/step,sigma);
        double st= flag.E_BSM_Put_Price();
        if(K-price[i][0] < st) price[i][1] = st;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    for(int i = 2 * step - 3 ; i >= 1; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            if(K-price[j][0] < st)
            {
                price[j][1] = st;
            }
            else
            {
                price[j][1] = K - price[j][0];
            }
            
        }
    }
    return price[0][1];

}
double Option_Pricing::A_TBSM_Put_Delta(int step)
{
    double price[2*step-1][2];
    double temp = So * exp(sigma * (step-1) * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step - 1; i++)
    {
        price[i][0] = temp;
        Option_Pricing flag(price[i][0],K,r,q,T/step,sigma);
        double st= flag.E_BSM_Put_Price();
        if(K-price[i][0] < st) price[i][1] = st;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    
    for(int i = 2 * step - 3 ; i >= 3; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            if(K-price[j][0] < st)
            {
                price[j][1] = st;
            }
            else
            {
                price[j][1] = K - price[j][0];
            }
            
        }
    }

    return (price[0][1]-price[2][1])/(price[0][0]-price[2][0]);
}
double Option_Pricing::A_TBSM_Put_Gamma(int step)
{
    double price[2*step-1][2];
    double temp = So * exp(sigma * (step-1) * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step - 1; i++)
    {
        price[i][0] = temp;
        Option_Pricing flag(price[i][0],K,r,q,T/step,sigma);
        double st= flag.E_BSM_Put_Price();
        if(K-price[i][0] < st) price[i][1] = st;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    for(int i = 2 * step - 3 ; i >= 5; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            if(K-price[j][0] < st)
            {
                price[j][1] = st;
            }
            else
            {
                price[j][1] = K - price[j][0];
            }
            
        }
    }
    return ((price[0][1]-price[2][1])/(price[0][0]-price[2][0])-(price[2][1]-price[4][1])/(price[2][0]-price[4][0]))/(price[1][0]-price[3][0]);

}
double Option_Pricing::A_TBSM_Put_Theta(int step)
{
    double price[2*step-1][2];
    double temp = So * exp(sigma * (step-1) * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step - 1; i++)
    {
        price[i][0] = temp;
        Option_Pricing flag(price[i][0],K,r,q,T/step,sigma);
        double st= flag.E_BSM_Put_Price();
        if(K-price[i][0] < st) price[i][1] = st;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    double temp1 = 0.0;
    for(int i = 2 * step - 3 ; i >= 1; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            if(K-price[j][0] < st)
            {
                price[j][1] = st;
            }
            else
            {
                price[j][1] = K - price[j][0];
            }
            
        }
        if(i == 3)
            temp1 = price[1][1];
        
    }
    return (temp1 - price[0][1])/T * step;

}

double Option_Pricing::E_TBSM_Put_Price(int step)
{
    double price[2*step-1][2];
    double temp = So * exp(sigma * (step-1) * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step - 1; i++)
    {
        price[i][0] = temp;
        Option_Pricing flag(price[i][0],K,r,q,T/step,sigma);
        double st= flag.E_BSM_Put_Price();
        price[i][1] = st;
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    for(int i = 2 * step - 3 ; i >= 1; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            price[j][1] = st;
            
            
        }
    }
    return price[0][1];
}
double Option_Pricing::E_TBSM_Put_Delta(int step)
{
    double price[2*step-1][2];
    double temp = So * exp(sigma * (step-1) * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step - 1; i++)
    {
        price[i][0] = temp;
        Option_Pricing flag(price[i][0],K,r,q,T/step,sigma);
        double st= flag.E_BSM_Put_Price();
        price[i][1] = st;
        
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    
    for(int i = 2 * step - 3 ; i >= 3; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            price[j][1] = st;
            
        }
    }
    return (price[0][1]-price[2][1])/(price[0][0]-price[2][0]);

}
double Option_Pricing::E_TBSM_Put_Gamma(int step)
{
    double price[2*step-1][2];
    double temp = So * exp(sigma * (step-1) * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step - 1; i++)
    {
        price[i][0] = temp;
        Option_Pricing flag(price[i][0],K,r,q,T/step,sigma);
        double st= flag.E_BSM_Put_Price();
        price[i][1] = st;
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    for(int i = 2 * step - 3 ; i >= 5; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            price[j][1] = st;
            
            
        }
    }
    return ((price[0][1]-price[2][1])/(price[0][0]-price[2][0])-(price[2][1]-price[4][1])/(price[2][0]-price[4][0]))/(price[1][0]-price[3][0]);

}
double Option_Pricing::E_TBSM_Put_Theta(int step)
{
    double price[2*step-1][2];
    double temp = So * exp(sigma * (step-1) * sqrt(3 * T/step));
    for(int i = 0; i < 2 * step - 1; i++)
    {
        price[i][0] = temp;
        Option_Pricing flag(price[i][0],K,r,q,T/step,sigma);
        double st= flag.E_BSM_Put_Price();
        price[i][1] = st;
        
        temp = temp * exp(- sigma * sqrt(3 * T/step));
    }
    double temp1 = 0.0;
    for(int i = 2 * step - 3 ; i >= 1; i = i - 2)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown,pmid;
            pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
            pmid = 2.0/3.0;
            price[j][0] = price[j+1][0];
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pmid + price[j + 2][1] * pdown);
            price[j][1] = st;
            
            
        }
        if(i == 3)
            temp1 = price[1][1];
        
    }
    return (temp1 - price[0][1])/T * step;
}

double Option_Pricing::A_TBSMRE_Put_Price(int step)
{
    return 2*A_TBSM_Put_Price(step)-A_TBSM_Put_Price(step/2);
}
double Option_Pricing::A_TBSMRE_Put_Delta(int step)
{
    return 2*A_TBSM_Put_Delta(step)-A_TBSM_Put_Delta(step/2);
}
double Option_Pricing::A_TBSMRE_Put_Gamma(int step)
{
    return 2*A_TBSM_Put_Gamma(step)-A_TBSM_Put_Gamma(step/2);
    
}
double Option_Pricing::A_TBSMRE_Put_Theta(int step)
{
    return 2*A_TBSM_Put_Theta(step)-A_TBSM_Put_Theta(step/2);

}

double Option_Pricing::E_TBSMRE_Put_Price(int step)
{
    return 2*E_TBSM_Put_Price(step)-E_TBSM_Put_Price(step/2);
}
double Option_Pricing::E_TBSMRE_Put_Delta(int step)
{
    return 2*E_TBSM_Put_Delta(step)-E_TBSM_Put_Delta(step/2);
}
double Option_Pricing::E_TBSMRE_Put_Gamma(int step)
{
    return 2*E_TBSM_Put_Gamma(step)-E_TBSM_Put_Gamma(step/2);
    
}
double Option_Pricing::E_TBSMRE_Put_Theta(int step)
{
    return 2*E_TBSM_Put_Theta(step)-E_TBSM_Put_Theta(step/2);
    
}

double Option_Pricing::A_BBSM_Put_Price(int step)
{
    double price[step][2];
    double temp = So * exp(sigma * (step-1) * sqrt(T/step));
    for(int i = 0; i < step; i++)
    {
        price[i][0] = temp;
        Option_Pricing flag(price[i][0],K,r,q,T/step,sigma);
        double st= flag.E_BSM_Put_Price();

        if(  st> K - price[i][0]) price[i][1] = st;
        else price[i][1] = K - price[i][0];
        temp = temp * exp(-2 * sigma * sqrt(T/step));
    }
    for(int i = step -1 ; i >= 1; i--)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown;
            pup = (exp((r-q) * T/step) - exp(- sigma * sqrt(T/step)))/(exp(sigma * sqrt(T/step)) - exp(-sigma * sqrt(T/step)));
            pdown =1-pup;
            price[j][0] = price[j][0]/exp(sigma * sqrt(T/step));
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pdown);
            if(K-price[j][0] < st)
            {
                price[j][1] = st;
            }
            else
            {
                price[j][1] = K - price[j][0];
            }
            
        }
    }
    return price[0][1];
}

double Option_Pricing::E_BBSM_Put_Price(int step)
{
    double price[step][2];
    double temp = So * exp(sigma * (step-1) * sqrt(T/step));
    for(int i = 0; i < step; i++)
    {
        price[i][0] = temp;
        Option_Pricing flag(price[i][0],K,r,q,T/step,sigma);
        double st= flag.E_BSM_Put_Price();
        price[i][1] = st;
        temp = temp * exp(-2 * sigma * sqrt(T/step));
    }
    for(int i = step -1 ; i >= 1; i--)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown;
            pup = (exp((r-q) * T/step) - exp(- sigma * sqrt(T/step)))/(exp(sigma * sqrt(T/step)) - exp(-sigma * sqrt(T/step)));
            pdown =1-pup;
            price[j][0] = price[j][0]/exp(sigma * sqrt(T/step));
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pdown);
            price[j][1] = st;
            
        }
    }
    return price[0][1];
}

double Option_Pricing::A_BBSMRE_Put_Price(int step)
{
    return 2*A_BBSM_Put_Price(step)-A_BBSM_Put_Price(step/2);
}

double Option_Pricing::E_BBSMRE_Put_Price(int step)
{
    return 2*E_BBSM_Put_Price(step)-E_BBSM_Put_Price(step/2);
}

vector<double> Option_Pricing::E_BADT_Pricer(int step)
{
   
    vector<double> lambda;
    lambda.reserve(step+1);
    lambda.push_back(1);
    double p_up = (exp((r-q) * T/step)-exp(-sigma * sqrt(T/step)))/(exp(sigma * sqrt(T/step))-exp(-sigma * sqrt(T/step)));
    double p_down = 1- p_up;
    for(int i = 2; i <= step+1; i++)
    {
        double flag = lambda[0];
        for(int j = 0; j < i-1; j++)
        {
            double flag1 = lambda[j];
            
            lambda[j] = lambda[j] * exp(-r * T/step)*p_up + (j != 0)*flag * exp(-r * T/step) * p_down;
            //cout<<setprecision(12)<<lambda[j]<<",";
            flag = flag1;
            
        }
        int j = i-1;
        lambda.push_back((j != 0)*flag * exp(-r * T/step) * p_down);
       //cout<<setprecision(12)<<lambda[j];
       //cout<<endl;
            
    }
    return lambda;
    
}

vector<double> Option_Pricing::E_TADT_Pricer(int step)
{
    double pup = 1.0/6.0 + (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
    double pdown = 1.0/6.0 - (r - q - sigma * sigma /2) * sqrt(T/step/12.0/sigma/sigma);
    double pmid = 2.0/3.0;
    vector<double> lambda1,lambda2;
   
    for(int i = 0; i < 2*step + 1; i ++)
    {lambda1.push_back(1);}
    lambda2 = lambda1;
    for(int i = 3; i <= 2*step+1; i = i + 2)
    {
        for(int j = 0; j < i; j++)
        {
            lambda1[j] = (j<i-2)*lambda2[j] * exp(-r * T/step)*pup + ((j < i-1 )&& (j > 0))*lambda2[j-1] * exp(-r * T/step) * pmid+(j > 1)*lambda2[j-2]*exp(-r * T/step) * pdown;
            //cout<<setprecision(12)<<lambda1[j]<<",";
            
        }
        lambda2 = lambda1;
        //cout<<endl;
    }
    return lambda1;
}

double Option_Pricing::A_BT_Call_Price(int step)
{
    
    ofstream f1("/Users/YuxiLIU/Desktop/out.csv");
    f1.clear();
    
   
    
    
    
    
    
    
    double price[step+1][2];
    double temp = So * exp(sigma * step * sqrt(T/step));
    for(int i = 0; i < step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] < K) price[i][1] = 0;
        else price[i][1] = price[i][0] - K;
        temp = temp * exp(-2 * sigma * sqrt(T/step));
    }
    for(int i = step ; i >= 1; i--)
    {
        double t = 1000000000;
        //double m = 0;
        for(int j = 0; j < i; j++)
        {
            double pup,pdown;
            pup = (exp((r-q) * T/step) - exp(- sigma * sqrt(T/step)))/(exp(sigma * sqrt(T/step)) - exp(-sigma * sqrt(T/step)));
            pdown =1-pup;
            price[j][0] = price[j][0]/exp(sigma * sqrt(T/step));
            double st = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pdown);
            if(price[j][0] - K < st)
            {
                price[j][1] = st;
               // cout<<"2";
            }
            else
            {
                price[j][1] = price[j][0] - K;
                if(price[j][0] < t) t = price[j][0];
                //if(price[j][0] > m) m = price[j][0];
                //cout<<"1";
            }
            
        }
        if(t != 1000000000)
            f1<<T/step * i<<","<<t<<endl;
        //if(m != 0)
            //f1 << T/step * i<<","<<m<<endl;
    }
     f1.close();
    return price[0][1];
}


double Option_Pricing::E_BT_Call_Price(int step)
{
    double price[step+1][2];
    double temp = So * exp(sigma * step * sqrt(T/step));
    for(int i = 0; i < step + 1; i++)
    {
        price[i][0] = temp;
        if(price[i][0] > K) price[i][1] = price[i][0] - K;
        else price[i][1] = 0;
        temp = temp * exp(-2 * sigma * sqrt(T/step));
    }
    for(int i = step ; i >= 1; i--)
    {
        for(int j = 0; j < i; j++)
        {
            double pup,pdown;
            pup = (exp((r-q) * T/step) - exp(- sigma * sqrt(T/step)))/(exp(sigma * sqrt(T/step)) - exp(-sigma * sqrt(T/step)));
            pdown =1-pup;
            price[j][1] = exp(-r * T /step) * (price[j][1] * pup + price[j + 1][1] * pdown);
            price[j][0] = price[j][0]/exp(sigma * sqrt(T/step));
        }
    }
    return price[0][1];

}

double Option_Pricing::E_BSM_Call_Price()
{
    return So * exp( -q * T) * N(d1) - K * exp( - r * T) * N(d2);
}

