// program to work on a 3 dimensional vector defined in maths
//math:  ai + bj + ck or xi + yj + zk
//prototype of the vector class using C++

#if 0
----------------------------------
  + addtion
  - subtraction
  * cross or vector product
     dot or scalar product
  / division by a number is scalar multiplication by its reciprocal
  == equality testing
     direction// ie its angle to the horizontal
     magnitude
     unit vector
     are_parallel, are_perpendicular
     angle bt two vectors
     get_component('x' or 'y' or 'z')
     multiply by scalar
     print
#endif


#include <iostream>
#include <cstdio>
#include <cmath>
#include "..\..\MACROS.h"
#include <iomanip>
#include <fstream>
#include <iomanip>
#include <exception>

#define EPSILON    0.001

class Vector{
	protected:
       // static int count;//just for fun; counts number of Vector objects in memory
       //check out later
	private:
		double x,//x - component
		y,// y - component
		z;// z - component

	public:
		Vector(double, double, double);
		//Vector(std::string)
		Vector();//for cases when it should be set later
		void set(double, double, double);
		double get_component(char);
		double getX();
		double getY();
		double getZ();
		Vector operator+(Vector);
		Vector operator-(Vector);
		Vector operator-()const {
			return Vector(-x, -y, -z );
		}//negate avector
		Vector operator*(Vector);//cross product
		double dotProduct(Vector);
		std::string toString(void);//essentially converts to string for printing to the console
		double magnitude(void);
		Vector operator*(double);//scalar multipliaction
		double direction(void);
		Vector unitVector(void);
		double distanceFrom(Vector);
		double angleTo(Vector);
		bool operator==(Vector);
		friend std::ostream& operator<<(std::ostream&, Vector);

};
Vector::Vector(double x, double y, double z = 0){
	this->x = x;
	this->y = y;
	this->z = z;
}
Vector::Vector(){
	x = y = z = 0;//at origin for unspecified vector
}
void Vector::set(double x, double y, double z){
    this->x = x;
	this->y = y;
	this->z = z;

}


double Vector::get_component(char c){
	switch(tolower(c)){//note there is no need for a break
		case 'x' : return x;
		case 'y' : return y;
		case 'z' : return z;
		default : return 0;//throw logic_error
	}
}

double Vector::getX(){
	return x;
}

double Vector::getY(){
	return y;
}

double Vector::getZ(){
	return z;
}
Vector Vector::operator+(Vector v){
	return Vector(x + v.x,  y + v.y,  z + v.z);
}
Vector Vector::operator-(Vector v){
    return Vector(x - v.x,  y - v.y,  z - v.z);
}

Vector Vector::operator*(Vector v){//cross product
	return Vector(((this->y * v.z)-(this->z*v.y)),  ((this->x * v.z)-(this->z*v.x)),
	((this->x * v.y)-(this->y*v.x))
	);
}

double Vector::dotProduct(Vector v){//dot product
	return (x*v.x + y*v.y + z*v.z);
}
std::string Vector::toString(void){//essentially converts to string for printing to the console
	char tmp[40];
	std::sprintf(tmp, "%.3fi %+.3fj %+.3fk", x, y, z);
	return std::string(tmp);
}
double Vector::magnitude(void){
	return double(std::sqrt(x*x + y*y + z*z));
}
Vector Vector::operator*(double a){
	return Vector(x*a, y*a, z*a);
}
Vector Vector::unitVector(void){//vector/magnitude
	double det = this->magnitude();
	return Vector(x/det, y/det, z/det);
}
double Vector::direction(void){//also angl
	//for 2d vectors
	return (std::atan2(y, x)*180)/M_PI;
}
double Vector::distanceFrom(Vector v){//check out later
	return (*this - v).magnitude();
}
double Vector::angleTo(Vector v){//also same as direction from
	return Vector(x-v.x, y-v.y, z-v.z).direction();
}
bool Vector::operator==(Vector v){
	return (std::fabs((x - v.x) < EPSILON) and std::fabs((x - v.x) < EPSILON) and
				std::fabs((x - v.x) < EPSILON));
}

std::ostream& operator << (std::ostream& os, Vector v){
	os << std::setprecision(3) << v.x <<"i " <<
	 std::showpos << v.y <<"j "<< v.z <<"k ";
	return os;
}
main(){
	Vector A(5,4,1), B(2,2, 0), C;
	std::cout << (-A + B)<< std::endl << std::endl;
	system("pause");
}











