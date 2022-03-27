// program to ionsn fraction prototype and work o

//math:  a/b, b != 0 and a and b are integers

//prototype of the Fraction class using C++



#if 0
----------------------------------
  + addtion
  - subtraction
  * multiplication by another and by a number
  /
  == equality testing
  +=
  -=
  ++
  ++
 ^  raise to the power for integers only.
#endif

#include <iostream>
#include <fstream>
#include "MACROS.h"
#include <cmath>

#define EPSILON    0.001



// ------------------------------------------------------
// Fraction.h
// A numerical class to represent fractions
// ------------------------------------------------------
#ifndef _FRACTION_
#define _FRACTION_
#include <cstdlib>
using namespace std;
class Fraction{
    protected:
       // static int count;//just for fun; counts number of Matrix objects in memory
       //check out later

	private:
		long numerator, denominator;
	public:
		Fraction(long n = 0, long d = 1);
		Fraction reciprocal(){
			return Fraction(denominator, numerator);
		}
		Fraction conjugate(){//same as reciprocal
			return Fraction(denominator, numerator);
		}
		Fraction operator-() const{
			return Fraction(-numerator, denominator);
		}
		Fraction& operator+=(const Fraction& a){
			numerator = a.numerator * denominator
					+ numerator * a.denominator;
			denominator *= a.denominator;
			return *this;
		}
		Fraction& operator-=(const Fraction& a){
			*this += (-a);
			return *this;
		}
		Fraction& operator++(){
			numerator += denominator;
			return *this;
		}
		Fraction& operator--(){
			numerator -= denominator;
			return *this;
		}
		Fraction operator^(int index){
			return Fraction(int(pow(numerator, index)), int(pow(denominator, index)));
		}
		double  value(){
			//the constructor handles problem of division by 0
			return  static_cast<double>(numerator)/denominator;
		}


		friend Fraction operator+(const Fraction&, const Fraction&);
		friend Fraction operator-(const Fraction&, const Fraction&);
		friend Fraction operator*(const Fraction&, const Fraction&);
		friend Fraction operator/(const Fraction&, const Fraction&);
		friend ostream& operator<< (ostream& os, const Fraction& a);
		friend istream& operator>> (istream& is, Fraction& a);
};
#endif
Fraction::Fraction(long n, long d){
	if(d == 0){
		cerr << "\nError: Division by zero!\n";
		exit(1);
	}
	if( n < 0 ) n = -n, d = -d;
	numerator = n; denominator = d;
}

Fraction operator+(const Fraction& a, const Fraction& b){
	Fraction temp;
	temp.denominator = a.denominator * b.denominator;
	temp.numerator = a.numerator*b.denominator + b.numerator * a.denominator;
	return temp;
}
Fraction operator-(const Fraction& a, const Fraction& b ){
	Fraction temp = a; temp += (-b);
	return temp;
}
Fraction operator*(const Fraction& a, const Fraction& b ){
	Fraction temp;
	temp.numerator = a.numerator * b.numerator;
	temp.denominator = a.denominator * b.denominator;
	return temp;
}
Fraction operator/(const Fraction& a, const Fraction& b ){
	if( b.numerator == 0){
		cerr << "\nError: Division by zero!\n";
		exit(1);
	}
	// To multiply a by the inverse of b:
	Fraction temp;
	temp.numerator = a.numerator * b.denominator;
	temp.denominator = a.denominator * b.numerator;
	if( temp.denominator < 0 )
		temp.numerator = -temp.numerator,
		temp.denominator = -temp.denominator;
	return temp;
}

ostream& operator<<(ostream& os, const Fraction& a){
	os << a.numerator << "/" << a.denominator;
	return os;
}
istream& operator>>(istream& is, Fraction& a){
	cout << "Enter a fraction:\n"
	" Numerator: "; is >> a.numerator;
	cout << " Denominator != 0: "; is >> a.denominator;
	if( !is) return is;
	if( a.denominator == 0){
		cout << "\nError: The denominator is 0\n"
		" New denominator != 0: ";
		is >> a.denominator;
		if( a.denominator == 0){
			cerr << "\nError: Division by zero!\n"; exit(1);
		}
	}
	if( a.denominator < 0 )
		a.numerator = -a.numerator,
		a.denominator= -a.denominator;
	return is;
}

int gcd(int a, int b){//use later to reduce factions
	if(b!=0)
	return abs(gcd(b, (a%b)));
	else
	return a;// was return 1; just to check
}
int main(){
	Fraction a(1,3), b(4);
	cout << "\nSome test results:\n\n";
	cout << " a = " << a << endl;
	cout << " b = " << b << endl;
	cout << " a + b = " << (a + b) << endl;
	cout << " a - b = " << (a - b) << endl;
	cout << " a * b = " << (a * b) << endl;
	cout << " a / b = " << (a / b) << endl;
	cout << " --a = " << --a << endl;
	cout << " ++a = " << ++a << endl;
	cout << " a ^ 4 = " << (a ^ 4)<< endl;
	a += Fraction(1,2);
	cout << " a+= 1/2; a = " << a << endl;
	a -= Fraction(1,2);
	cout << " a-= 1/2; a = " << a << endl;
	cout << "-b = " << -b << endl;
	cout << "\nAnd now an input\n";
	cin >> a;
	cout << "\nYour input: " << a << endl << endl;
	system("pause");
	return 0;
}
