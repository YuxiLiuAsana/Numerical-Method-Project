#ifndef MATRIX_H
#define MATRIX_H
#include<iostream>
#include<vector>
#include<cmath>
using namespace std;

template <class T>
class Matrix
{
public:
	vector<vector<T>> Data;
    vector<vector<T>> TData;
	vector<vector<T>> L;
	vector<vector<T>> U;
	vector<int> P;
    vector<vector<T>> IData;
    

public:
	Matrix();
	Matrix(vector<vector<T>>&);
	Matrix(Matrix<T>&);
	~Matrix();

	Matrix<T>& operator=(Matrix<T>&);
	Matrix<T> operator+(Matrix<T>&) const;
	Matrix<T> operator-(Matrix<T>&) const;
	Matrix<T> operator*(Matrix<T>&) const;
	Matrix<T>& operator=(vector<vector<T>>&);
	Matrix<T> operator+(vector<vector<T>>&) const;
	Matrix<T> operator-(vector<vector<T>>&) const;
	Matrix<T> operator*(vector<vector<T>>&) const;
	vector<T> operator*(vector<T>&) const;
	friend  ostream& operator<<(ostream &out, Matrix<T>& source) 
	{
		for(int i = 0; i < source.Data.size(); i++)
		{
			for(int j = 0; j < source.Data[0].size(); j++)
			{
				out << source.Data[i][j]<<" ";
			}
			out << endl;
		}
		return out;	
	}
    
    friend ostream& operator<<(ostream & out, vector<T> & source)
    {
        for(int i = 0; i < source.size()-1;i++)
        {
            out<<source[i] <<",";
        }
        out<<source[source.size()-1]<<endl;
        return out;
    }

	bool Is_Matrix() const;
	bool Is_Square_Matrix() const;
	bool Is_Singular() const;
	bool Is_Symatric() const;
	bool Is_SPD() const;
	bool Is_SPSD() const;

	vector<vector<T>> get() const;
	void set(vector<vector<T>>&);
	void Transpose();
	void Inverse();
	T Determinant() const;
	vector<T> eigenvalue() const;
	vector<vector<T>> eigenvector() const;
	vector<T> eigenvector(T&) const;
	

	vector<T> forward_subst(const vector<T>&) const;
	vector<T> backward_subst(const vector<T>&) const;
	void lu_no_pivoting();
	void lu_row_pivoting();
    void cholesky_decomposition();
    vector<T> linear_solve_cholesky(const vector<T>&);
    vector<T> linear_solve_Lu_no_pivoting(const vector<T>&);
    vector<T> linear_solve_Lu_row_pivoting(const vector<T>&);
    vector<T> Jacobi(vector<T>&) const;
    vector<T> Gauss_Siedel(vector<T>&) const;
    vector<T> SOR(vector<T>&, double ) const;
    vector<T> SORP(vector<T>&, double, vector<T>& ) const;

	
};



#endif