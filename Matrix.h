#pragma once
// This function is used for check the information of the matrix

#ifndef MATRIX_H
#define MATRIX_H
#include<stdlib.h>
#include<string>
#include<iostream>
using namespace std;

// define abstract class
class FloatMatrix; 
class Intarry;   
//super class

// sub-class
class  Matrix
{
protected:
	size_t nRows;     // avoid use int, int has the minus number
	size_t nColumns;

public:
	 Matrix() {};     // default constructor
	 Matrix(int n, int m) { nRows = n; nColumns = m; }
	 ~Matrix() {};	   // destructor

	void checkBounds(int, int);								  //check the boundary of matrix
	int giveNumberOfRows()   const { return nRows; }         // check the row of matrix
	int giveNumberOfColumns()  const { return nColumns; }	 // check the column of matrix
	bool isSquare() const { return (!(nRows - nColumns)); }  // check square matrix
};

Matrix::Matrix()
{
}

Matrix::~Matrix()
{
}

 void Matrix::checkBounds(int i, int j)
 {
	 if (i <= 0) 
	 {
		 cout << i << endl;
		 throw logic_error("matrix error on row <= 0 \n");
	}
	 if (j <= 0)
	 {
		 cout << j << endl;
		 throw logic_error("matrix error on column <= 0 \n");
	 }
	 if (i > nRows)
	 {
		 cout << i << endl;
		 throw logic_error("matrix exceeds the row > nRow");
	 }
	 if (j > nColumns)
		 cout << j << endl;
		throw logic_error("matrix exceeds the column > nColumn");
 }


#endif // !MATRIX_H
