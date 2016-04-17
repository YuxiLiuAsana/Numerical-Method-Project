//
//  FD.cpp
//  Finite Different Method
//
//  Created by LiuYuxi on 12/8/15.
//  Copyright © 2015 LiuYuxi. All rights reserved.
//

#include "FD.hpp"

FD :: FD() : x_left( 0 ) , x_right( 0 ) , t_final( 0 ) , M( 0 ) , N( 0 ) , alpha( 0 )
{

}

FD :: FD(double sx_left, double sx_right, double st_final, int sM, int sN, double salpha) : x_left( sx_left ) , x_right( sx_right ) , t_final( st_final ) , M( sM ) , N( sN ) , alpha( salpha )
{
    vector<double> b;
    for(int i = 0; i < N + 1; i ++)
    {
        b.push_back(0);
    }
    for(int i = 0; i < M + 1; i++)
        data.push_back(b);
        
    for(int i = 0; i < N+1; i++)
    {
        data[0][i] = f(i);
    }
    for(int i = 0; i < M +1 ; i++)
    {
        data[i][0] = g_left(i);
        data[i][N] = g_right(i);
    }

}

FD :: ~FD()
{

}

double FD :: f( int i )
{
    return exp( -2.0 + (double)i *4.0 / N ) ;
}

double FD :: g_left( int i )
{
    return exp( (double)i /M - 2 );
}

double FD :: g_right( int i )
{
    return exp( (double)i /M + 2 );
}

double FD :: u_exact(int i, int j)
{
    return exp( -2.0 + (double)i *4.0 / N + (double)j /M );
}

void FD :: ForwardEuler()
{
    
    for(int i = 1; i < M + 1; i++)
    {
        for( int j = 1; j< N; j++)
        {
            data[i][j] = alpha * data[i - 1][j - 1] + ( 1 - 2 * alpha) * data[i-1][j] + alpha * data[i - 1][j + 1];
            
        }
    }
}

void FD :: BackwardEuler(string method)
{
    Matrix<double> a;
    //cout<<alpha<<endl;
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
            result = a.SOR(b,1.2);
        
        for(int j = 1; j < N; j++ )
            data[i][j] = result[j-1];
    }

}

void FD :: CrankNicolson(string method)
{
    Matrix<double> a,b1;
    //cout<<N<<endl;
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
            result = a.linear_solve_Lu_no_pivoting(b);
        else if(method == "SOR")
            result = a.SOR(b,1.2);
        
        for(int j = 1; j < N; j++ )
            data[i][j] = result[j-1];
    }
}

double FD :: MaxPointwiseError()
{
    double max = 0;
    for(int i = 1; i < M+1;i++)
        for(int j = 1; j < N; j++)
        {
            double exact = u_exact(j, i );
            if(abs(exact-data[i][j])>max)
                max = abs(exact-data[i][j]);
        }
    return max;
}

double FD::ErrorRMS()
{
    double sum = 0;
    for(int i = 0; i < N; i ++)
    {
        double exact = u_exact(i,M);
        sum +=(data[M][i]-exact)*(data[M][i]-exact)/exact/exact;
    }
    sum = sum/(1+N);
    sum = sqrt(sum);    
    return sum;
}


ostream& operator<<(ostream &out, FD& source)
{
    for(int i = 0 ; i< source.M + 1; i ++)
    {
        for(int j = 0; j < source.N ; j++)
            out<<setprecision(12)<<source.data[i][j]<<", ";
        out<<setprecision(12)<<source.data[i][source.N]<<endl;
    }
    return out;
}