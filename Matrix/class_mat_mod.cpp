// program to work on matrices as defined in maths
//math: rows X columns
//prototype of the Matrix class using C++


//make 2 definitions of the matrix class
//one which the constructor initialises and one initialised by methods

/**********************************************************************************************
//This matrix class and some others will be later implemented using STLwhen well versed*******
//eg. std::vector

//ADD column and row marginals to the matrix to string once finished                                                                  ***********
***********************************************************************************************/

//modify programs to throw exception rather than logging error messages
//and assume in trivial cases

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



//all constructs using multidimensional arrays will be changed to
//pointers!~~be careful!!!!!!!!!

//or preferably vectors


/**********************************
**this version uses vectors or two d arrays to initialise the matrix
**************************************/

#include <iostream>
#include <conio.h>
#include <cstdio>
#include <cmath>
#include <fstream>
#include "..\\..\\MACROS.h"
#include <iomanip>
#include <sstream>
#include <exception>

#define EPSILON    0.00001

#include<vector>
using std::vector;
using std::cout;
using std::string;
class Matrix{
	protected:
        static int count;//just for fun; counts number of Matrix objects in memory
       //check out later
	private:
		int nrows,//number of rows
		ncols;//number of columns
		bool size_set ;
		double* M;//will contain the elements of the matrix
		bool is_square(void){
			return ncols==nrows;
		}
		bool is_identity(void);

	public:
		Matrix(double[][10], int, int);//rows by colum
		Matrix();
		Matrix(vector<vector<double> >, int r, int c);//n_rows and cols to be specified

		void set_size(int, int);
		void set_matrix(double [][10]);
		void set_values(std::istream from = std::cin);
		std::string toString(void);//essentially converts to string for printinging to the console


		//can use magnitude()  or determinant() to get |A| where A is a matrix
		double determinant(void);
		double magnitude();
		Matrix inverse(void);
		Matrix adjoint(void);
		Matrix cofactor_matrix(void);
 		Matrix minor_submatrix(int , int);
		
		//double operator[](int, int);
		double operator()(int ,int);
		
		
		Matrix multiply_by(Matrix);//multiply elememt_wise
		
		Matrix transpose(void);
		
		Matrix operator^(int pow);
		
		Matrix operator/(double);
		
		friend Matrix operator+(Matrix&, vector<vector<double> >);
		friend Matrix operator+(Matrix&, Matrix&);
		
		friend Matrix operator-(Matrix&, vector<vector<double> >);
		friend Matrix operator-(Matrix&, Matrix&);
		
		friend Matrix operator*(Matrix&, vector<vector<double> >);
		friend Matrix operator*(Matrix&, Matrix&);
		friend Matrix operator*(double, Matrix&);
		
		friend bool operator==(const Matrix&, const Matrix&);
		friend std::ostream& operator<<(std::ostream&, const Matrix);
		
		
		~Matrix();//destructor
};
Matrix::Matrix(int m, int n, double Arr[][10]){
	this->nrows = m;
	this->ncols = n;
	this->size_set = true;
	this->M = new double[m][n];
	for(int i = 0; i < this->nrows; ++i){
		for(int j = 0; j < this->ncols; j++)
		    M[i][j] = Arr[i][j];
	}
	size_set = true;
}
Matrix::Matrix(){
	nrows = ncols = 2;
	size_set = false;
}
void Matrix::set_size(int m, int n){
	nrows = m;
	ncols = n;
	M = new double[m][n];
	size_set = true;
}
void Matrix::set_matrix(double a[][10]){//really not necessary
	if(size_set)
	    for(int i = 0; i < nrows; ++i)
			for(int j = 0; j < ncols; j++)
			    M[i][j] = a[i][j];
   else std::cout << "Matrix size not set , DEFAULT";
}
void Matrix::reset_values(std::ifstream from = std::cin){//stringstream or file
	//set the matrix size
	if(!size_set){
		std::cout<<"The size of the matrix has not been set:: DEFAULT =2  X 2"<<std::endl;
	}
	std::cout<<"Getting Values ..."<<std::endl;
	for(int i = 0; i < nrows; ++i){
		for(int j = 0; j < ncols; ++i)
		    from >> M[i][j]
	}
	std::cout<<"Values have been set";
}

Matrix Matrix::operator+(Matrix B){
	if(nrows != B.nrows or ncols != B.ncols){
		//THOW AN EXCEPTION HERE : MATRIX ADDITION COMPATIBILITY ERROR
	}
	double ret[nrows][ncols];
	for(int i = 0; i < nrows; ++i){
		for(int j = 0; j < ncols; j++)
		    ret[i][j] = M[i][j] + B.M[i][j];
	}
	return Matrix(nrows, ncols, ret);
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
	else {
		//return sum of deterinanst of the sub matrices
		//which may in turn have sub matrices
		double det = 0;
		for(int j = 0; j < ncols; ++j){
			//use the first row to do the. math::the row remains the same the column xhnages
		//just put the 0 explicitely even though it representds the first row
		    det += int(pow(-1, (j+0)))*this->M[0][j]*(this->minor_submatrix(0,j)).determinant();
		    det += ((j+0)%2 ? 1 : -1 )*this->M[0][j]*(this->minor_submatrix(0,j)).determinant();
		  }
		return det;//returning the determinant
	}
}

Matrix Matrix::minor_submatrix(int i , int j){//starting from(0,0)
	//minor submatrix is the matrix gotten after deleting the
	//the row and column specified
	//Return the submatrix obtained by removing the `i`th row
     //and `j`th column from ``self``.

	if(i<0 or j<0 ){
		//THOW AN EXCEPTION HERE : MATRIX INDEXING ERROR
		// "MATH ERROR : Matrix Indices Cannot be Negative\n\n";
	
		}
	if(i >= this->nrows or j >= this->ncols){
        //THOW AN EXCEPTION HERE : MATRIX INDEXING  ERROR
		//"MATH ERROR : Matrix Indices Cannot be greater than size\n\n";
		

	}
	if( not this->is_square()){
        std::cerr << "MATH ERROR : Minor of a non square matrix is questionable ??\n\n";
		getch();
		exit(EXIT_FAILURE);

	}
	if(not size_set){
        std::cerr << "MATH ERROR : Matrix Size not given\n\n";
		getch();
		exit(EXIT_FAILURE);

	}
	//skips the row and column specified(ie deletes it fro the new matrix)
	double temp[nrows-1][ncols-1];
	int newcol, newrow;
	newcol=newrow = 0;
	for(int row = 0; row < nrows-1; ++row){
		if (row == i)
		    continue;
		for(int col = 0; col < ncols-1; col++){
		    if (col==j)continue;
		    temp[newrow][newcol]= this->M[row][col];
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
	if (nrows == 1)
	  return (*this);
	else if(nrows == 2){//just do for two since it's easy
		double ret[2][2];
		for (int i =0; i <2 ; ++i)
			for(int j =0; j < 2; j++)
			    ret[i][j]= (i==j? M[1-i][1-j] : -M[i][j]);
		return Matrix(2,2,ret);
	}
	else{
    	double ret[nrows][ncols];
		for (int i =0; i < nrows ; ++i)
			for(int j =0; j < ncols; j++)
			    ret[i][j]= int(pow(-1, (i+j)))*this->M[i][j]*
			    (this->minor_submatrix(i ,j)).determinant();
		return Matrix(nrows, ncols, ret);
	}

}


//check again
Matrix Matrix::transpose(void){
	double ret[ncols][nrows];
	for(int i = 0; i < ncols; ++i){
		for(int j = 0; j < nrows; j++)
		    ret[i][j] = M[j][i];
	}
	return Matrix(ncols, nrows, ret);
}

double Matrix::operator[](int i, int j){
	if(i>=nrows or j >= ncols or i < 0  or j < 0){
		TEXTCOLOR(FOREGROUNG_RED);
		std::cerr << "MATH ERROR: Indices are invalid or out of range" << "\n\tEXITING";
		TEXTCOLOR(DEFAULT);
		getch();
		exit(EXIT_FAILURE);
	}

	return M[i][j];
}

Matrix operator*(const Matrix& A, const Matrix& B){
	if(A.ncols != A.nrows){
		return make_identity(A.nrows, A.ncols);
	}
	//final matrix
	double mat[A.nrows][B.ncols];
    for(int i = 0; i < A.nrows; ++i){
		for(int j = 0; j < B.ncols; j++){
			int sum = 0;
			for(int k = 0; k < A.ncols/*or B.nrows*/ ; ++k )
				sum += A.M[i][k]*B.M[k][j];
		    mat[i][j] = sum;
		}
	}
	return Matrix(A.nrows, B.ncols, mat);
}

 Matrix operator*(double scalar, const Matrix& A){
	double ret[A.nrows][A.ncols];
	for(int i = 0; i < nrows; ++i){
		for(int j = 0; j < ncols; j++)
		    ret[i][j] = scalar * A.M[i][j];
	}
	return (A.nrows, A.ncols, ret)
 }



//global function to generate a mathematical identity matrix

Matrix make_identity(int size){ //dimension of the matrix
	double temp[rows][cols];
	//vector<vector<double>> temp;
	for(int i = 0; i < size; ++i){
		for(int j = 0; j < size; j++)
		    temp[i][j] =(i==j) ? 1 : 0;
	}
	return Matix(rows, cols, temp);
}



std::ostream& operator<<(std::ostream& os, const Matrix  mat){
//matrix will be printed for now in trhe space of 3 numbers and the
//order will be lost if numbers are greater than 3 digits~~????
	os << "\n\t";
	for(int in=0; in<mat.ncols; ++in)
		os << " [" << in << "] ";
	os<<std::endl;
	for(int i=0; i<mat.ncols; ++i){
	    os << "[" << i << "]  ";
		for (int j=0; j<mat.rows; ++j){
			os << std::setw(3) << " " << mat.M[j][i] << " ";
		}
		os << std::endl;
	}
}



///testing phase
//
//
//
//


int main(){
	double arr[3][3]={
		{1, 2,-1},
		{3, 0, 1},
		{4, 2, 1}
	};
	Matrix mat(3,3,arr);
	std::cout<<mat;
	getch();
}


