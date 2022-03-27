// program to work on the Complex datatype defined in maths
//math:  a + bi  ---  engineering  a + bj, a and b are real numbers
//prototype of the Complex class using C++
#if 0
----------------------------------
  + addtion
  - subtraction
  * multiplication
  / division
  ^ raised to the power(integers only)
	conjugate
	is real for getting the real number from it
 == equality testing
#endif
#include <iostream>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <conio.h>
#define EPSILON    0.001

class Complex{
	protected:
		//count+=1;
        //static int count;//just for fun; counts number of Complex objects in memory
       //check out later
	private:
		double real,//real part of the Complex number
		imag;// imaginary part of the Complex number
		
	public:
		Complex(double, double);
		Complex();//for cases when itit should be set later
		Complex(double re):real(re), imag(0){
		}//when a single argument is supplied
		//Complex(std::string);//for later
		~Complex(){
			//do nothing
		}
		Complex(const Complex& c):real(c.real), imag(c.imag){
		}
		Complex& operator= (const Complex& c){//assignment operator
			real = c.real;
			imag = c.imag;
			return *this;
		}
		Complex& operator() (const Complex& c){
			this -> real = c.real;
			this -> imag = c.imag;
			return *this;
		}
		void set_values(double, double);
		void set_values(const Complex z);
		
		double get_real(void) const;
		double get_imag(void) const;
		
		std::string to_string(void);//essentially converts to string for printing to the console
		Complex conjugate(void)const;
		double abs(void)const;
		double magnitude(void)const;
		
		bool is_real(void){return std::fabs(imag) < EPSILON;}
		
		Complex operator^(int);// exponentiation
		Complex& operator^=(int);// exponentiation in place
		
		Complex operator-() const{
			return Complex(-real, -imag);
		}
		
		friend Complex operator+(const Complex, const Complex);
		friend Complex& operator+=(Complex&, const Complex);
		
		friend Complex operator-(const Complex, const Complex);
		friend Complex& operator-=(Complex&, const Complex);
		
		friend Complex operator*(const Complex, const Complex);
		friend Complex& operator*=(Complex&, const Complex);
		
		friend Complex operator/(const Complex, const double);
		friend Complex& operator/=(Complex&, const double );

		friend Complex operator/(const Complex, const Complex);
		friend Complex& operator/=(Complex&, const Complex);
		
		friend bool operator!=(const Complex&, const Complex&);
		friend bool operator==(const Complex&, const Complex&);
		friend std::ostream& operator<<(std::ostream& out, const Complex z);
		//an input definition for the complex class will not be performed directly
		//using the folllowing definition for reasons of originality and naturality
		//friend istream& operator>>(istream& in, const Complex z)
};

//Complex::count=0;

Complex::Complex(){
    real = imag = 0;
	//sets values to zero if not given
}

Complex::Complex(double x, double y):real(x), imag(y){
	//do nothing
}

void Complex::set_values(double a, double b){
	real = a;
	imag = b;
}
void Complex::set_values(const Complex z){
	real = z.get_real();
	imag = z.get_imag();
}
double Complex::get_real(void)const{
	return real;
}
double Complex::get_imag(void)const{
	return imag;
}

double Complex::abs(void) const{
	//return std::sqrt(std::pow(real, 2)+std::pow(imag, 2));
	return hypot(real, imag);
}

double Complex::magnitude() const {
	return abs();
}
Complex Complex::conjugate(void)const{
	return Complex(real, -imag);
}

Complex Complex::operator^(int power){//still skeptical of this implementation
	Complex z = *this;
	for(int i = 1; i <= power; ++i)
		z = z**this;
	return z;
}

Complex& Complex::operator^=(int power){//still skeptical of this implementation
	Complex z = *this;
	for(int i = 1; i <= power; ++i)
		z = z**this;
	return (*this = z);
}

//makes a suitable string to be printed to the console reprenting the Complex number
//can still use the ostringstream buf will use in other class

std::string Complex::to_string(void){
//    std::ostringstream ostr;
//	//output can be written directly to a string-->> implement completely later
//	ostr << setprecision(3) <<real <<" - " << imag <<" i ";
//	return (ostr.str());//converts it in to a class string object

	char temp[20];
	//sprintf takes a char array not a pointer to char
	if(imag > 0)
		std::sprintf(temp, "%.2g + %.3g i",real, imag);
	else
		std::sprintf(temp, "%.2g - %.3g i",real, std::fabs(imag));
	return std::string(temp);
}


Complex operator+(const Complex c, const Complex z){
	return Complex(c.real + z.real, c.imag + z.imag);
}
Complex& operator+=(Complex& z, const Complex c){
	z.real += c.real;
	z.imag += c.imag;
	return z;
}



Complex operator-(const Complex c, const Complex z){
	return Complex(c.real - z.real, c.imag - z.imag);
}
Complex& operator-=(Complex& z, const Complex c){
	z = z - c;
	return z;
}



Complex operator*(const Complex c, const Complex z){
    double re = (c.real*z.real) - (c.imag * z.imag),
	im = (c.real*z.imag) + (c.imag * z.real );
	return Complex(re, im);

}
Complex& operator*=(Complex& c1, const Complex c2){
	c1 = c1 * c2;
	return c1;
}



Complex operator/(const Complex z, const double div){
	return Complex(z.real/div, z.imag/div);
}

Complex& operator/=(Complex& c, const double div){
	c.real /= div;
	c.imag /= div;
	return c;
}



Complex operator/(const Complex c, const Complex z){
	double denom = (z*z.conjugate()).get_real();
	return ((c * z.conjugate()) / denom);
}
Complex& operator/=(Complex& c, const Complex z){
	double denom = (z*z.conjugate()).get_real();
	c = c * z.conjugate();
	c /= denom;
	return c;
}

bool operator==(const Complex& c, const Complex& z){
	return (std::fabs(c.real-z.real)<EPSILON and std::fabs(c.imag-z.imag)<EPSILON);
}

bool operator!=(const Complex& a, const Complex& b){
	return !(a==b);
}


std::ostream& operator<<(std::ostream& out, const  Complex z){
	const char* sign = z.imag < 0 ? " - " : " + ";
	out << std::setprecision(3) << z.real << sign << std::fabs(z.imag) <<" i ";
	return out;
}


//Time to test the hard work

main(int argc, char *argv[]){
	using std::cout;
	using std::endl;
	
	Complex z(2, 1), u(1, 1);
	Complex a = 2.3, b = 1/a;
	Complex c = a+b*Complex(1, 2.3),
	p = (a/b)+2;
	cout<< z <<endl;
	cout<< z.conjugate() <<endl;
	cout << z / u << endl;
	cout << (u^4) <<endl;
	cout<<std::boolalpha<<(z==u);
	cout << a <<endl;
	cout << b <<endl;
	cout << c <<endl;
	cout << p <<endl;
	cout << (1+2*p)/c <<endl;
	
	Complex z1(10,5),z3(2,2) ,z2=z3/Complex(3,-1);
	cout<<"[-----------------------------------]";
	cout<<z1<<endl;
	cout<<z2<<endl;
	cout<<z1*z2<<endl;
	cout<<z1*z3;
	getch();
}

//improvement
//.. Make constructor to take a string representing a complex number as in math::
	//extend this concept to other classws like matrix and vectors
	//also to take direct input from an istream
// .. check divisin by 0
// .. check division by a 0 complex number(0 +0i)








