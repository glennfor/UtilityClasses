
// program to work on matrices as defined in maths
//math: rows X columns
//prototype of the Matrix class using C++

/**********************************************************************************************
//implement container with vector or list DS
//implement numbers as fraction objects
//ADD column and row marginals to the matrix to string representation
//through exceptiions instead of logging messages                                                     
***********************************************************************************************/

#if 0
----------------------------------
  + addtion
  - subtraction
  * Matrix multiplication
	 multiply_by(multiplies element-wise)

  == equality testing
     magnitude
     inverse
     transpose
     is_identity
     is_square, is_rectangular, is_diagonal, is_mult_identity
     get_determinant or get_magnitude
     toString
     cofactor_matrix
     get_adjoint//also get_adjugate
     multiply by scalar
 ^  raise to the power
#endif

//for a 2-D array eg  Mat[3][3] to access elements by means of pointers
//Mat[i][j] == *(Mat + j + ncols*i) ncols here being 3


#include <iostream>
#include <conio.h>
#include <cstdio>
#include <cmath>
#include <fstream>
#include "..\\..\\MACROS.h"
#include <iomanip>
#include <sstream>

#include <exception>
#include <stdexcept>
#include <ctime>
#include <cstdlib>

#include <vector>

#include <string>
#include <cstring>

#define EPSILON    0.00001 // e-5


using std::cout;
using std::cin;
using std::cerr;

using std::string;
using std::endl;
using std::vector;

using std::exception;

//function to allocate a pointer to a pinter
//ie a 2-d array


void allocate(double** &pointer, int rows, int columns){
	pointer = (double**)new double[rows];
	for(int i=0; i < columns; ++i)
	    pointer[i] = (double*)new double [columns];
}

	//pointer = (double**)malloc(sizeof(double*)*rows);
	//for(int i=0; i < columns; ++i)
	    //pointer[i] = (double*)malloc(sizeof(double)*rows);

void deallocate(double** &pointer, int rows, int columns){
	for(int i=0; i < columns; ++i)
	    delete[] pointer[i];
	delete[] pointer;
	pointer = 0;
}



//exception class to throw exceptions
class MathException: public exception{
	private:
		string message;
	public:

		~MathException()throw(){
		};

		MathException(std::string msg):message(msg){
		}
		MathException(const char* msg):message(msg){
		}
		virtual const char* what() const throw();
		virtual std::string what_is() {//const throw()
			return message;
		}
		void setMessage(std::string msg){
            this->message=msg;
		}
 		void setMessage(const char* msg){
            this->message=msg;
		}

};

const char* MathException::what()const throw(){
	return this->message.c_str();
}


class Matrix{
    friend Matrix operator*(const Matrix&, const Matrix&)throw (std::string);
	friend Matrix operator*(double, const Matrix&);
	friend bool operator==(const Matrix&, const Matrix&);
	friend std::ostream& operator<<(std::ostream&, const Matrix&);
	friend Matrix makeIdentityMatrix(int );

	protected:
        static int count;//just for fun; counts number of Matrix objects in memory
       //check out later
	private:
		int nrows,//number of rows
		ncols;//number of columns
		bool size_set ;
		double** M;//will contain the elements of the matrix
		bool is_square(void){
			return ncols==nrows;
		}
		bool is_identity(void);

	public:
		Matrix(int, int, double**);//rows by colum
		Matrix();//for cases when itit should be set later
		Matrix(int r, int c)
		{
			int nrows(r), ncols(c);
			allocate(M, nrows,ncols);
		}
		void set_size(int, int);
		void set_matrix(double [][10]);
		void set_values(std::istream& from = std::cin);
		Matrix operator+(Matrix);
		Matrix operator-(Matrix);
		Matrix operator/(double);
		
		double at(int, int);
		bool is_singular(void){
			return fabs(this->determinant()) < EPSILON;
		}
		Matrix multiply_by(Matrix);//multiply elememt_wise
		Matrix transpose(void);
		std::string toString(void);//essentially converts to string for printinging to the console
	
		Matrix minor_submatrix(int , int);//returns a matrix object representing the
		//minor submatrix
		
		//can use magnitude()  or determinant() to get |A| where A is a matrix
		double determinant(void);
		double magnitude();
		Matrix inverse(void);
		Matrix adjoint(void);
		Matrix adjugate(void);
		Matrix cofactor_matrix(void);
		Matrix operator^(int pow);
		Matrix power(int pow);
		Matrix exp(int pow);
		double operator()(int ,int);
		~Matrix();//destructor

};


Matrix makeIdentityMatrix(int);

//maybe this constructor should deallocate the array
Matrix::Matrix(int r, int c, double** Arr){
	//note:: Arr is a 2-D matrix
	this->nrows = r;
	this->ncols = c;
	this->size_set = true;
	allocate(this->M, this->nrows, this->ncols);
	for(int i = 0; i < this->nrows; ++i){
		for(int j = 0; j < this->ncols; j++)
		    M[i][j] = Arr[i][j];
	}
}


Matrix::Matrix(){
	nrows = ncols = 2;
	this->M = 0;
	size_set = false;
}

void Matrix::set_size(int m, int n){
	allocate(this->M, this->nrows, this->ncols);
	size_set = true;
}


Matrix::~Matrix(){
	for(int i = 0; i < this->ncols; ++i)
	    delete[] this->M[i];
	delete[] this->M;
	this->M = 0;
}


void Matrix::set_matrix(double a[][10]){//really not necessary
	if(size_set)
	    for(int i = 0; i < nrows; ++i)
			for(int j = 0; j < ncols; j++)
			    M[i][j] = a[i][j];
   else std::cout << "Matrix size not set , DEFAULT";
}


//modify here for the matrix to state its size in the prompt
void Matrix::set_values(std::istream& from /*= std::cin*/){
	//stringstream or file
	//set the matrix size
	if(!this->size_set){
		std::cout << "Unset Size";
		nrows=ncols=2;
		std::cout << nrows << "x" << ncols << std::endl;
	}
	std::cout<<"Getting Values ..."<<std::endl;
	for(int i = 0; i < nrows; ++i)
	{
		for(int j = 0; j < ncols; ++j)
		    from >> M[i][j];
	}
}


Matrix Matrix::operator+(Matrix B){
   	if(nrows != B.nrows or ncols != B.ncols){
		//THOW AN EXCEPTION HERE : MATRIX ADDITION COMPATIBILITY ERROR
	}
	
	double** ret;
	allocate(ret, nrows, ncols);
	for(int i = 0; i < nrows; ++i){
		for(int j = 0; j < ncols; j++)
		    ret[i][j] = M[i][j] + B.M[i][j];
	}
	return Matrix(nrows, ncols, ret);
}


Matrix Matrix::operator/(double div){
	if(fabs(div) < EPSILON)
	{
		throw MathException("Division By Zero");
	}
	
	double** ret;
	allocate(ret, this->nrows, this->ncols);
	for(int i = 0; i < this->nrows; ++i)
		for(int j = 0; j < this->ncols; j++)
		    ret[i][j] = this->M[i][j]/div;
	
	return Matrix(nrows, ncols, ret);
}


//inverse calculation

Matrix Matrix::inverse(void){
	if(this->is_singular())
	{
        //throw new MathException(std::string("Math Error: Singular matrix cannot have and inverse"));
		//raise exception here
	}

	return this->adjoint()/this->determinant();
}



//calculating the determinant
// done recursively
double Matrix::determinant(void){//same as the magnitude
	if(not this->is_square()){

		//THOW AN EXCEPTION HERE : MATRIX  COMPATIBILITY ERROR
		//"MATH ERROR: Determinant of a non square matrix can not be calculated"
	}
	
	if(ncols==nrows and nrows==1)
	    return this->M[0][0];
	else if(ncols==2){
		return this->M[0][0]*this->M[1][1] - this->M[1][0]*this->M[0][1];
	}
	else {
		double det = 0;
		for(int j = 0; j < ncols; ++j){
			//use the first row to do the math::the row remains the same the column xhnages
			//just put the 0 explicitely even though it representds the first row
		    det += pow(-1, (j+0))*this->M[0][j]*(this->minor_submatrix(0,j)).determinant();
		    //det += ((j+0)%2 ? 1 : -1 )*this->M[0][j]*(this->minor_submatrix(0,j)).determinant();
		  }
		return det;//returning the determinant
	}
}



Matrix Matrix::minor_submatrix(int i , int j){//starting from(0,0)
	//minor submatrix is the matrix gotten after deleting the
	//the row and column specified
	//Return the submatrix obtained by removing the `i`th row
     //and `j`th column 

	if(i<0 or j<0 ){
		//THOW AN EXCEPTION HERE : MATRIX INDEXING ERROR
		// "MATH ERROR : Matrix Indices Cannot be Negative\n\n";

		}
	if(i >= this->nrows or j >= this->ncols){
        //THOW AN EXCEPTION HERE : MATRIX INDEXING  ERROR
		//"MATH ERROR : Matrix Indices Cannot be greater than size\n\n";

		}
	if(not size_set){
        //THOW AN EXCEPTION HERE : MATRIX SIZE error in index ERROR
	}
	//skips the row and column specified(ie deletes it fro the new matrix)
	double** temp;
	allocate(temp, nrows-1, ncols-1);
	int newcol, newrow;
	newcol = newrow = 0;
	for(int row = 0; row < nrows; ++row)
	{
		if (row == i)
		    continue;
		for(int col = 0, newcol=0; col < ncols; col++)
		{
		    if (col == j)
				continue;
		    temp[newrow][newcol] = this->M[row][col];
		    newcol++;
		}
		newrow++;
	}

	return Matrix(nrows-1, ncols-1, temp);
}



//calculaing the cofactor matrix
// done recursively
Matrix Matrix::cofactor_matrix(void){
	if(not this->is_square() and not size_set){
        //THOW AN EXCEPTION HERE : MATRIX COMPATIBILITY ERROR
		// "MATH ERROR: Cofactor matrix of a non square matrix can not be calculated" << "\n\tEXITING";
	}
	if (this->nrows == 1)
	  return (*this);
	else if(this->nrows == 2){//just do for two since it's easy
		double** ret;
		allocate(ret, 2, 2);
		for (int i =0; i < this -> nrows ; ++i)
			for(int j =0; j < this -> ncols; j++)
			    ret[i][j]= (i==j? M[1-i][1-j] : -M[i][j]);
		return Matrix(2,2,ret);
	}
	else{
		double** ret;
		allocate(ret, this->nrows, this->ncols);
		for (int i =0; i < this->nrows ; ++i)
			for(int j = 0; j < this->ncols; j++)
			    ret[i][j]= pow(-1, (i+j))*(this->minor_submatrix(i ,j)).determinant();

		return Matrix(nrows, ncols, ret);
	}

}



Matrix Matrix::adjoint(void){
	return this->cofactor_matrix().transpose();
}


Matrix Matrix::adjugate(void){
	return this->adjoint();
}


//check again
Matrix Matrix::transpose(void){
	double** ret;
	allocate(ret,ncols, nrows);
	for(int i = 0; i < ncols; ++i){
		for(int j = 0; j < nrows; j++)
		    ret[i][j] = this->M[j][i];
	}
	return Matrix(ncols, nrows, ret);
}


//get the item at a particular index(0-based indexing)
double Matrix::at(int i, int j){
	
	bool unqualified_indices = (i >= this->nrows or j >= this->ncols
								 or i < 0  or j < 0 );

	if(unqualified_indices){
		
		char error[100],*msg=
		"Matrix Error:"
		"Matrices indices out of range: Matrix Size is(%dx%d).";
		std::sprintf(error,msg, nrows, ncols);
		throw (MathException(error));
	}
	return this->M[i][j];
}


Matrix operator*(const Matrix& A, const Matrix& B)throw (std::string){
	if(A.ncols != B.nrows){
		char *error,  msg[] =
		"Matrix Error: Matrices of sizes %dx%d and %dx%d cannot be multiplied. ";
		std::sprintf(error, msg, A.nrows, A.ncols, B.nrows, B.ncols);
		//throw  MathException(error);
		throw new std::string(error);
 	}
	//final matrix
	double** mat;
	allocate(mat, A.nrows, B.ncols);
    for(int i = 0; i < A.nrows; ++i)
		for(int j = 0; j < B.ncols; j++)
		{
			int sum = 0;
			
			for(int k = 0; k < A.ncols/*or B.nrows*/ ; ++k )
				sum += A.M[i][k]*B.M[k][j];
				
		    mat[i][j] = sum;
		}

	return Matrix(A.nrows, B.ncols, mat);
}



Matrix operator*(double scalar, const Matrix& A){
	double** ret;
	int nrows = A.nrows,
	ncols = A.ncols;
	allocate(ret, nrows, ncols);
	for(int i = 0; i < nrows; ++i){
		for(int j = 0; j < ncols; j++)
		    ret[i][j] = scalar * A.M[i][j];
	}
	return Matrix(nrows, ncols, ret);
 }



//global function to generate a mathematical identity matrix

Matrix makeIdentityMatrix(int size){
	double** temp;
	allocate(temp, size,size);
	for(int i = 0; i < size; ++i){
		for(int j = 0; j < size; j++)
		    temp[i][j] =(i==j) ? 1 : 0;
	}
	return Matrix(size, size, temp);
	//deallocate(temp, size, size);

}


///CORIGINAL LAM
#if 0 //ss
std::ostream& operator<< (std::ostream& os, const Matrix&  mat){
//matrix will be printed for now in trhe space of 3 numbers and the
//order will be lost if numbers are greater than 3 digits~~????
	os << "\n\n    ";
	TEXTCOLOR(9);
	for(int in=0; in < mat.ncols; ++in)
		os << " [" << in << "] ";
	os << std::endl;
	for(int i = 0; i < mat.nrows; ++i){
        TEXTCOLOR(9);
	    os << "[" << i << "] ";
	    TEXTCOLOR(10);
		for(int j = 0; j < mat.ncols; ++j)
			os << std::setw(3) <<std::right << mat.M[i][j] <<"  ";
		os << std::endl;
	}
	os << std::endl;
	TEXTCOLOR(DEFAULT);
}

#endif //if 0

std::ostream& operator<< (std::ostream& os, const Matrix&  mat){
//matrix will be printed for now in trhe space of 3 numbers and the
//order will be lost if numbers are greater than 3 digits~~????
	os << "\n\n     ";
	TEXTCOLOR(9);
	for(int in=0; in < mat.ncols; ++in)
		os << "  [" << in << "]  ";
	os << std::endl << std::endl;
	for(int i = 0; i < mat.nrows; ++i){
        TEXTCOLOR(9);
	    os << "[" << i << "]  ";
	    TEXTCOLOR(10);
		for(int j = 0; j < mat.ncols; ++j)
			os << std::setw(4) <<std::right
			   << (fabs(mat.M[i][j]) < EPSILON ? 0.0 : mat.M[i][j]) <<"   ";
		os << std::endl <<std::endl;
	}
	os << std::endl;
	TEXTCOLOR(DEFAULT);
}






///testing phase
//
//
//
//



int main(){
	cout << "MATRIX MANIPULATION ~~\n\n";
	srand(time(NULL));
	double** arr, **arr2;
	allocate(arr2, 2, 2);
	allocate(arr, 3,3);

	for(int i = 0; i < 3; ++i)
	    for(int j = 0; j < 3; ++j)
	        arr[i][j] = rand()%10;
	        
	Matrix mat(3,3,arr), mx(2, 2, arr2),
			mat4 = makeIdentityMatrix(4);
			
	cout << mat;
	cout << "\nEnter elemnt of a 3x3 matrix: \n";
	mat.set_values();


	cout << std::endl << "New Matrix is :" << mat;
	cout << std::endl << "Determinat is :" << mat.determinant();
	cout << std::endl << "Transpose is" << mat.transpose();
	cout << std::endl << "Matrix*3/2 " << 3*mat/2;
	cout << std::endl << "Matrix+Matrix " << mat+mat;
	cout << endl << "Matrix Adjoint is " << mat.adjoint();
	cout << "sub matrix of mat at [0,0] " << mat.minor_submatrix(0,0).determinant();
	cout << std::endl << "Cofactor matrix is" << mat.cofactor_matrix();
    cout << std::endl << "Matrix4 " << mat4;


    
	cout << endl << endl;
	
	cout << std::endl << "Enter 2X2 matrix :\n";
	//mx.set_values();
	cout << "\nMatrix is: " << mx;
	cout << "Cofact is " << mx.cofactor_matrix();
	//cout << std::endl << "Adjoint  is : " << mat.adjugate();
	//cout << std::endl << "Inverse is : " << mat.inverse();
	try{
		//cout << "\n\n->" ;
		cout << "Element at (" << 1 << ", " << 1 << ") = " << mat.at(2,1);
	}catch(MathException& ext){
        TEXTCOLOR(FOREGROUND_RED);//|FOREGROUND_INTENSITY);
		cerr << "  Exception Ocurred:\n  " << ext.what();
		TEXTCOLOR(FOREGROUND_RED|FOREGROUND_GREEN|FOREGROUND_BLUE);
	}catch(...){
		cerr << "elixis!!";
		exit(EXIT_FAILURE);
	}
	
	deallocate(arr2, 2, 2);
	deallocate(arr, 3,3);
	getch();
	return 0;
}


/*private members of a class are accessible only from within
 other members of the same class or from
their friends*/

//interface
