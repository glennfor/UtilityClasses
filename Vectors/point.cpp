
#include <iostream>
#include <cstdio>
#include <cmath>
#include "..\..\MACROS.h"
#include <iomanip>
#include <fstream>
#include <iomanip>
#include <exception>

#define EPSILON    0.001
using std::ostream ;
using std::istream;
using std::cout;
using std::cin;

class Point
{
	private:
	    double m_dX, m_dY, m_dZ;

	public:
	    Point(double dX=0.0, double dY=0.0, double dZ=0.0)
	    {
	    m_dX = dX;
	    m_dY = dY;
	    m_dZ = dZ;
	    }

	    friend ostream& operator<< (ostream &out, Point &cPoint);
	    friend istream& operator>> (istream &in, Point &cPoint);

	    double GetX() { return m_dX; }
	    double GetY() { return m_dY; }
	    double GetZ() { return m_dZ; }
};

ostream& operator<< (ostream &out, Point &cPoint)
{
    out << "(" << cPoint.m_dX << ", " <<
        cPoint.m_dY << ", " <<
        cPoint.m_dZ << ")";
    return out;
}

istream& operator>> (istream &in, Point &cPoint)
{
    in >> cPoint.m_dX;
    in >> cPoint.m_dY;
    in >> cPoint.m_dZ;
    return in;
}




main(){
	cout << "Space point calculations: \n\n";
	cout << "Entre a point :";
	Point p;
	cin >> p;
	cout << "\nPoint is  " << p <<"\n\n" ;
	system("pause");
}
