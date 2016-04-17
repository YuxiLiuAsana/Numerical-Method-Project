#include "Matrix.h"
#include<iomanip>
template <class T>
Matrix<T>::Matrix()
{
	Data.clear();
}

template <class T>
Matrix<T>::Matrix(vector<vector<T>>& source )
{
	Data = source;
}

template <class T>
Matrix<T>::Matrix(Matrix<T>& source)
{
	Data = source.Data;
}

template <class T>
Matrix<T>::~Matrix()
{
	Data.clear();
}

template <class T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>& source)
{
	Data = source.Data;
	return *this;
}


template <class T>
Matrix<T> Matrix<T>::operator+(Matrix<T>& source) const
{
	Matrix<T> original;
	original.Data = Data;
	if(source.Data.size() != Data.size())
	{
		cout<<"Different Row Size"<<endl;
		original.Data.clear();
	}
	else
	{
		for(int i = 0; i < Data.size(); i ++)
		{
			if(source.Data[i].size() != Data[i].size())
			{
				cout<<"Different Line Size!!!"<<endl;
				original.Data.clear();
			}
			else
			{
				for(int j = 0; j < Data[i].size(); j++)
				{
					original.Data[i][j] = Data[i][j] + source.Data[i][j];
				}
			}
		}
	}
	return original;
}

template <class T>
Matrix<T> Matrix<T>::operator-(Matrix<T>& source) const
{
	Matrix<T> original;
	original.Data = Data;
	if(source.Data.size() != Data.size())
	{
		cout<<"Different Row Size"<<endl;
		original.Data.clear();
	}
	else
	{
		for(int i = 0; i < Data.size(); i ++)
		{
			if(source.Data[i].size() != Data[i].size())
			{
				cout<<"Different Line Size!!!"<<endl;
				original.Data.clear();
			}
			else
			{
				for(int j = 0; j < Data[i].size(); j++)
				{
					original.Data[i][j] = Data[i][j] - source.Data[i][j];
				}
			}
		}
	}
	return original;
}

template <class T>
Matrix<T> Matrix<T>::operator*(Matrix<T>& source) const
{
	
	Matrix<T> original;
	vector<double> temp;
	if(Data[0].size() != source.Data.size())
	{
		cout<<" The matrixs can not multiply!!!"<<endl;
		original.Data.clear();
	}
	else
	{
		for(int i = 0; i < Data.size(); i ++)
		{		
			temp.clear();
			for(int j = 0; j < source.Data[0].size(); j ++)
			{
				double sum = 0;
				for(int k = 0; k < Data[0].size(); k ++)
				{
					sum += Data[i][k] * source.Data[k][j];
				}
				temp.push_back(sum);
			}
			original.Data.push_back(temp);
		}
	}
	return original;
}

template <class T>
Matrix<T>& Matrix<T>::operator=(vector<vector<T>>& source)
{
	Data = source;
	return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator+(vector<vector<T>>& source) const
{
	Matrix<T> original ;
	original.Data = Data;
	if(source.size() != Data.size())
	{
		cout<<"Different Row Size"<<endl;
		original.Data.clear();
	}
	else
	{
		for(int i = 0; i < Data.size(); i ++)
		{
			if(source[i].size() != Data[i].size())
			{
				cout<<"Different Line Size!!!"<<endl;
				original.Data.clear();
			}
			else
			{
				for(int j = 0; j < Data[i].size(); j++)
				{
					original.Data[i][j] = Data[i][j] + source[i][j];
				}
			}
		}
	}
	return original;
}

template <class T>
Matrix<T> Matrix<T>::operator-(vector<vector<T>>& source) const
{
	Matrix<T> original;
	original.Data =	Data;
	if(source.size() != Data.size())
	{
		cout<<"Different Row Size"<<endl;
		original.Data.clear();
	}
	else
	{
		for(int i = 0; i < Data.size(); i ++)
		{
			if(source[i].size() != Data[i].size())
			{
				cout<<"Different Line Size!!!"<<endl;
				original.Data.clear();
			}
			else
			{
				for(int j = 0; j < Data[i].size(); j++)
				{
					original.Data[i][j] = Data[i][j] - source[i][j];
				}
			}
		}
	}
	return original;
}

template <class T>
Matrix<T> Matrix<T>::operator*(vector<vector<T>>& source) const
{
	
	Matrix<T> original;
	vector<double> temp;
	if(Data[0].size() != source.size())
	{
		cout<<" The matrixs can not multiply!!!"<<endl;
		original.Data.clear();
	}
	else
	{
		for(int i = 0; i < Data.size(); i ++)
		{		
			temp.clear();
			for(int j = 0; j < source[0].size(); j ++)
			{
				double sum = 0;
				for(int k = 0; k < Data[0].size(); k ++)
				{
					sum += Data[i][k] * source[k][j];
				}
				temp.push_back(sum);
			}
			original.Data.push_back(temp);
		}
	}
	return original;
}

template <class T>
ostream& operator<<(ostream &out, vector<vector<T>>& source) 
{
	for(int i = 0; i < source.size(); i++)
	{
		for(int j = 0; j <source[i].size(); j++)
		{
			out << source[i][j]<<" ";
		}
		out << endl;
	}
	return out;	
}

template<class T>
vector<T> Matrix<T>::operator*(vector<T>& source) const
{
    vector<T> answer;
	if(Data[0].size()!= source.size())
	{
		cout<<"Cannot be multiplied"<<endl;
		answer.clear();
		return answer;
	}
	else
	{
		T sum;
		
		for(int i = 0; i < Data.size();i++)
		{
			sum = 0;
			for(int j = 0; j < Data[i].size(); j++)
			{
				sum += Data[i][j] * source[j];
			}
			answer.push_back(sum);
			
		}
		return answer;
	}
}

template<class T>
bool Matrix<T>::Is_Matrix() const
{
	int i;
	for(i = 0; i < Data.size(); i++)
	{
		if(Data[i].size() != Data[0].size())
			break;
	}
	if(i != Data.size())		
		return false;
	else
		return true;
}

template<class T>
bool Matrix<T>::Is_Square_Matrix() const
{
    bool a = Is_Matrix();
    if(a == false|| Data[0].size() != Data.size())
		return false;
	else
	return true;
}

template<class T>
bool Matrix<T>::Is_Singular() const
{

}
	
template<class T>
bool Matrix<T>::Is_Symatric() const
{

}
	
template<class T>
bool Matrix<T>::Is_SPD() const
{
    return true;
}

template<class T>
bool Matrix<T>::Is_SPSD() const
{

}


template <class T>
void  Matrix<T>::set(vector<vector<T>>& source)
{
	Data.clear();
	Data = source;
}

template <class T>
vector<vector<T>> Matrix<T>::get() const
{
	return Data;
}

template <class T>
void Matrix<T>::Transpose()
{
    TData.clear();
    for(int i = 0 ; i < Data[0].size(); i ++)
    {
        vector<T> b;
        b.clear();
        for(int j = 0; j < Data.size(); j ++)
            b.push_back(Data[j][i]);
        TData.push_back(b);
    }
}

template <class T>	
void Matrix<T>::Inverse()
{
    vector<vector<T>> ind=Data;
    IData = Data;
    for(int i = 0; i < Data.size();i++)
    {
        for(int j = 0; j < Data.size(); j++)
        {
            ind[i][j] = 0;
            if(i == j)
                ind[i][j] = 1;
        }
        
    }
    for(int i = 0; i < Data.size();i++)
    {
        IData[i] = linear_solve_Lu_row_pivoting(ind[i]);
    }
    Matrix<T> temp(IData);
    temp.Transpose();
    IData = temp.TData;
    
}

template <class T>
T Matrix<T>::Determinant() const
{
	return 1;
}

template <class T>	
vector<T> Matrix<T>::eigenvalue() const
{

}

template <class T>	
vector<vector<T>> Matrix<T>::eigenvector() const
{

}

template <class T>	
vector<T> Matrix<T>::eigenvector(T& eigenvalue) const
{

}

template<class T>
vector<T> Matrix<T>::forward_subst(const vector<T>& b) const
{
	vector<T> x;
	x.clear();
	if(Is_Matrix() == false)
	{
		cout<<"It is not a matrix!!!"<<endl;
		return x;
	}
	else if(Data.size() != Data[0].size())
	{
		cout<<"It is not a square matrix!!!"<<endl;
		return x;
	}
	else
	{
		int flag = 0;
		for(int i = 0; i < Data.size(); i++)
		{
			for(int j = i + 1; j < Data[i].size(); j++)
				if( Data[i][j] != 0)
					flag ++;
		}
		if(flag != 0)
		{
			cout<<"It is not a lower matrix!!!"<<endl;
			return x;
		}
		else if(b.size() != Data[0].size())
		{
			cout<<"Dimension of L and b is not fit!!!"<<endl;
			return x;
		}
		else
		{
			x.push_back(b[0]/Data[0][0]);
			for(int i = 1; i < b.size(); i++)
			{
				T sum = 0;
				for(int j = 0; j < i; j++)
					sum +=Data[i][j]*x[j];
				x.push_back((b[i]-sum)/Data[i][i]);
			}
			return x;
		}
	}
}

template <class T>
vector<T> Matrix<T>::backward_subst(const vector<T>& b) const
{
	vector<T> x = b;
	
	if(Is_Matrix() == false)
	{
		cout<<"It is not a matrix!!!"<<endl;
		x.clear();
		return x;
	}
	else if(Data.size() != Data[0].size())
	{
		cout<<"It is not a square matrix!!!"<<endl;
		x.clear();
		return x;
	}
	else
	{
		int flag = 0;
		for(int i = 0; i < Data.size(); i++)
		{
			for(int j = i + 1; j < Data[i].size(); j++)
				if( Data[j][i] != 0)
					flag ++;
		}
		if(flag != 0)
		{
			x.clear();
			cout<<"It is not a upper matrix!!!"<<endl;
			return x;
		}
		else if(b.size() != Data[0].size())
		{
			x.clear();
			cout<<"Dimension of L and b is not fit!!!"<<endl;
			return x;
		}
		else
		{
			x[b.size()-1]=b[b.size()-1]/Data[b.size()-1][b.size()-1];
			for(int i = b.size()-2; i >= 0; i--)
			{
				T sum = 0;
				for(int j = i+1; j < b.size(); j++)
				sum +=Data[i][j]*x[j];				
				x[i]=(b[i]-sum)/Data[i][i];
			}
			return x;
		}
	}
}

template<class T>
void Matrix<T>::lu_no_pivoting()
{
    P.clear();
    L.clear();
    U.clear();
    if(Is_Square_Matrix() == false)
		cout<<"It is not a square matrix !!!"<<endl;
	else 
	{
		int i;
		for(i = 0; i < Data.size(); i++)
		{
			vector<vector<T>> a;
			for(int j = 0; j<=i; j++)
			{
				vector<T> b;
				for(int k =0; k <= i; k++)
					b.push_back(Data[j][k]);
				a.push_back(b);
			}
			Matrix<T> c(a);
			if(c.Determinant()==0)
				break;
		}
		if(i != Data.size())
			cout<<"The matrix cannot be LU decomposited without row pivoting!!!"<<endl;
		else
		{
			L = Data;
			U = Data;
			vector<vector<T>> Temp(Data);
			for(int i = 0; i < Data.size(); i ++)
			{
				for(int j = i+1; j < Data.size(); j++)
				{
					L[i][j] = 0;
					U[j][i] = 0;
				}
				L[i][i] = 1;
			}
			for(int i = 0; i < Data.size()-1; i++)
			{
				U[i][i]=Temp[i][i];
				for(int j = i+1; j < Data.size(); j ++)
				{
					U[i][j]=Temp[i][j];
					L[j][i]=Temp[j][i]/U[i][i];
				}
				for(int j = i + 1; j < Data.size(); j++)
				{
					for(int k = i+1; k < Data.size(); k++)
					{
						Temp[j][k]=Temp[j][k]-L[j][i]*U[i][k];
					}
				}
			}
			U[Data.size()-1][Data.size()-1]=Temp[Data.size()-1][Data.size()-1];
		}
	}
}

template<class T>
void Matrix<T>::lu_row_pivoting()
{
    P.clear();
    L.clear();
    U.clear();
    if(Is_Square_Matrix() == false)
		cout<<"It is not a square matrix !!!"<<endl;
	else if(Determinant() == 0)
		cout<<"The matrix is singular!!!"<<endl;
	else
	{
		L = Data;
		U = Data;
		vector<vector<T>> Temp(Data);
		for(int i = 0; i < Data.size(); i ++)
		{
			P.push_back(i + 1);
		}
		for(int i = 0; i < Data.size(); i ++)
		{
			for(int j = i+1; j < Data.size(); j++)
			{
				L[i][j] = 0;
				U[j][i] = 0;
			}
			L[i][i] = 1;
		}
		for(int i = 0; i < Data.size()-1; i++)
		{
			T value = abs(Temp[i][i]);
			int flag = i;
			for(int j = i+1; j < Data.size();j++)
			{
				if(value < abs(Temp[j][i]))
				{
					flag = j;
					value = abs(Temp[j][i]);
				}
			}
			if(flag != i)
			{
				vector<T> p;
				p = Temp[i];
				Temp[i] = Temp[flag];
				Temp[flag] = p;
				T m;
				for(int j = 0; j < i; j++)
				{
					m = L[i][j];
					L[i][j] = L[flag][j];
					L[flag][j] = m;
				}
				int t;
				t = P[i];
				P[i] = P[flag];
				P[flag] = t;
			}
			U[i][i] = Temp[i][i];
			for(int j = i + 1; j < Data.size(); j++)
			{
				U[i][j]=Temp[i][j];
				L[j][i]=Temp[j][i]/U[i][i];
			}
			for(int j = i + 1; j < Data.size(); j ++)
			{
				for(int k = i + 1; k < Data.size(); k ++)
				{
					Temp[j][k] = Temp[j][k]-L[j][i]*U[i][k];
				}
			}
		}
		U[Data.size()-1][Data.size()-1]=Temp[Data.size()-1][Data.size()-1];

	}
}

template<class T>
void Matrix< T >::cholesky_decomposition()
{
   
    if(Is_SPD() == false)
    {
        cout<< "The matrix is not symetric positive definite!!!"<<endl;
    }
    else
    {
        vector<vector<T>> temp(Data);
        U.clear();
        U = Data;
        for(int i = 1 ; i < Data.size();i++)
            for(int j = 0; j < i ; j++)
                U[i][j] = 0;
        for(int i = 0; i < Data.size()-1 ; i++)
        {
            U[i][i] = sqrt(temp[i][i]);
            for(int j = i+1; j < Data.size(); j++)
                U[i][j] = temp[i][j]/U[i][i];
            for(int j = i + 1; j < Data.size(); j++)
            {
                for(int k = j ; k < Data.size(); k ++)
                {
                    temp[j][k] = temp[j][k] - U[i][j]*U[i][k];
                }
            }
        }
        U[Data.size()-1][Data.size()-1] = sqrt(temp[Data.size()-1][Data.size()-1]);
    }
 
}

template<class T>
vector<T> Matrix<T>::linear_solve_cholesky(const vector<T>& source) 
{
    cholesky_decomposition();
    Matrix<T> a(U),b;
    a.Transpose();
    b = a.TData;
    vector<T> y, x;
    y = b.forward_subst(source);
    x = a.backward_subst(y);
    return x;
}

template<class T>
vector<T> Matrix<T>::linear_solve_Lu_no_pivoting(const vector<T>& source)
{
    lu_no_pivoting();
    Matrix<T> a(U),b(L);
   
    vector<T> y, x;
    y = b.forward_subst(source);
    x = a.backward_subst(y);
    return x;
}

template<class T>
vector<T> Matrix<T>::linear_solve_Lu_row_pivoting(const vector<T>& source1)
{
    lu_row_pivoting();
    
    Matrix<T> a(U),b(L);
    vector<T> source = source1;
    for(int i = 0; i < source1.size(); i++)
    {
        source[i]=source1[P[i]-1];
        //cout<<source[i]<<","<<P[i]-1<<","<<source1[i]<<endl;
    }
    
    vector<T> y, x;
    y = b.forward_subst(source);
    x = a.backward_subst(y);
    return x;
}

template<class T>
vector<T> Matrix<T>::Jacobi(vector<T>& b) const
{
    double tol = 1e-6;
    vector<T> x0;
    vector<vector<T>> temp = Data;
    for(int i = 0; i< Data.size(); i++)
    {
        temp[i][i] = 0;
        x0.push_back(0);
    }
    vector<T> r0 = b;
    vector<T> A_x0 = *this * x0;
    for(int i = 0 ; i < Data.size(); i++)
    {
        r0[i] = b[i] - A_x0[i];
    }
    double rnorm = 0;
    for( int i = 0; i < Data.size(); i++)
    {
        rnorm += r0[i] * r0[i];
    }
    rnorm = sqrt(rnorm);
    tol = tol * rnorm;
    Matrix<T> flag(temp);
    int iter = 0;
    int first = 0;
     cout<<"Jacobi"<<endl;
    while (rnorm > tol)
    {
        iter ++;
        vector<T> sum = flag * x0;
        for(int i = 0; i < Data.size(); i++)
        {
            x0[i] = -1/Data[i][i] * sum[i] + 1/Data[i][i] *b[i];
        }
        A_x0 = *this * x0;
        for(int i = 0; i < Data.size(); i++)
        {
            r0[i] = b[i] - A_x0[i];
        }
        rnorm = 0;
        for(int i = 0; i < Data.size(); i++)
        {
            rnorm += r0[i] * r0[i];
        }
        rnorm = sqrt(rnorm);
       
        if(first < 3)
        {
            for(int i = 0; i < x0.size();i++)
                cout<<setprecision(11)<<x0[i]<<',';
            cout<<endl;
        }
        
        first ++;
        
    }
    cout<<iter<<", "<<setprecision(12)<<rnorm<<endl;
    return x0;
}

template<class T>
vector<T> Matrix<T>::Gauss_Siedel(vector<T>& b) const
{
    double tol = 1e-6;
    vector<T> x0;
    vector<vector<T>> LD = Data, UD = Data ;
    for(int i = 0; i< Data.size(); i++)
    {
        x0.push_back(0);
        for( int j = i + 1; j < Data.size(); j++)
        {
            LD[i][j] = 0;
            UD[j][i] = 0;
        }
        UD[i][i] = 0;
    }
    Matrix<T> DL(LD),DU(UD);
    vector<double> b_new = DL.forward_subst(b);
    vector<T> r0 = b;
    vector<T> A_x0 = *this * x0;
    for(int i = 0 ; i < Data.size(); i++)
    {
        r0[i] = b[i] - A_x0[i];
    }
    double rnorm = 0;
    for( int i = 0; i < Data.size(); i++)
    {
        rnorm += r0[i] * r0[i];
    }
    tol = tol * sqrt(rnorm);
    int iter = 0;
    int first = 0;
    cout<<"GS"<<endl;
    while(rnorm > tol)
    {
        iter++;
        vector<T> m = DL.forward_subst(DU * x0);
        for(int i = 0; i< Data.size(); i++)
        {
            x0[i] = -m[i] + b_new[i];
        }
        vector<T> A_x0 = *this * x0;
        for(int i = 0; i < Data.size(); i++)
        {
            r0[i] = b[i] - A_x0[i];
        }
        rnorm = 0;
        for(int i = 0; i < Data.size(); i++)
        {
            rnorm += r0[i] * r0[i];
        }
        rnorm = sqrt(rnorm);
        
        
        first ++;
    }
    cout<<iter<<", "<<setprecision(12)<<rnorm<<endl;
    return x0;

}

template<class T>
vector<T> Matrix<T>::SOR(vector<T>& b, double w) const
{
    double tol = 1e-6;
    vector<T> x0;
    vector<vector<T>> LD = Data, UD = Data ;
    for(int i = 0; i< Data.size(); i++)
    {
        x0.push_back(0);
        for( int j = i + 1; j < Data.size(); j++)
        {
            LD[i][j] = 0;
            UD[j][i] = 0;
            LD[j][i] = LD[j][i]*w;
            UD[i][j] = -UD[i][j]*w;
        }
        UD[i][i] = (1-w)*UD[i][i];
        
    }
    Matrix<T> DL(LD),DU(UD);
    vector<double> b_new = DL.forward_subst(b);
    vector<T> r0 = b;
    double rnorm = 0;
    for( int i = 0; i < Data.size(); i++)
    {
        rnorm += r0[i] * r0[i];
    }
    tol = tol * sqrt(rnorm);
    int iter = 0;
  
    
    while(true )
    {
        vector<T> x_t=x0;
        iter++;
        vector<T> m = DL.forward_subst(DU * x0);
        for(int i = 0; i< Data.size(); i++)
        {
            x0[i] = m[i] + w * b_new[i];
        }
        vector<T> A_x0 = *this * x0;
        for(int i = 0; i < Data.size(); i++)
        {
            r0[i] = b[i] - A_x0[i];
        }
        rnorm = 0;
        for(int i = 0; i < Data.size(); i++)
        {
            rnorm += (x0[i]-x_t[i]) * (x0[i]-x_t[i]);
        }
        rnorm = sqrt(rnorm);
        if(rnorm < 1e-6)
            break;      
        
       
    }
    
    return x0;
}

template<class T>
vector<T> Matrix<T>::SORP(vector<T>& b, double w, vector<T>& x0) const
{
    //fix the tol
    double tol = 1e-8;
    
    vector<T> x_old = x0;
    vector<T> x_new = x0;
    
    while(true)
    {
        x_old = x_new;
        for(int i = 0; i < x0.size(); i++)
        {
            double sum1 = 0, sum2 = 0;
            for(int j = 0; j < i ; j++)
            {
                sum1 += Data[i][j] * x_new[j];
            }
            for(int j = i+1; j < x0.size(); j++)
            {
                sum2 += Data[i][j] * x_old[j];
            }
            x_new[i] = (1-w) * x_old[i] - w/Data[i][i] * (sum1 + sum2 ) + w * b[i]/Data[i][i];
        }
    
        double sum = 0;
        for(int i = 0; i < x_old.size(); i ++)
        {
            sum += (x_new[i] - x_old[i]) * (x_new[i] - x_old[i]);
        }
        sum = sqrt(sum);
        if(sum < tol)
            break;
    }
    return x_new;
}