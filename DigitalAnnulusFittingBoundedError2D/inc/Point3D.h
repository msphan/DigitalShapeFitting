#ifndef POINT3D_H_
#define POINT3D_H_

#include <iostream>
#include <vector>
#include <math.h>
#include <cstdlib>

using namespace std;

class Point3D{

private:
	float x,y,z;
public:
	Point3D(){}
	Point3D(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
	float getX();
	float getY();
	float getZ();
	void setX(float _x);
	void setY(float _y);
	void setZ(float _z);
	void translation(float a,float b,float c);
	void translation(Point3D t);
	void inverse();
	double innerProduct(Point3D t);
	double outerProduct(Point3D t);
	double squareDis();
	double squareDis(Point3D t);
	void permulation(int type);
	bool isInside(vector<Point3D> ps); //TODO : not implement

};

#endif