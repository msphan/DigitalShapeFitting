#ifndef PLANE_H_
#define PLANE_H_

#include "Point3D.h"

class Plane{
protected:
	// ax + by + cz +d = 0
	Point3D normal;
	float d;
public:
	Plane(){}
	Plane(Point3D _normal, float _d) : normal(_normal), d(_d){}
	Plane(float _a, float _b, float _c, float _d) : normal(_a,_b,_c), d(_d){}
	float getNormalX();
	float getNormalY();
	float getNormalZ();
	Point3D getNormal();
	float getIntercept();
	void setNormalX(float _x);
	void setNormalY(float _y);
	void setNormalZ(float _z);
	void setNormal(Point3D _normal);
	void setNormal(float _x, float _y, float _z);
	void setIntercept(float _d);
	double signedDis(Point3D p); //TODO : What is the signed distance
	int principalProjectionPlane(); //TODO : What is the principal projection plane
	// count the number of points is on and below the distance d
	int count(vector<Point3D> ps, float d); 
	// count the number of points is on and below the distance d
	// the results is save in set of vector in
	int count_report(vector<Point3D> ps, float d, vector<int>& in); 
private:	
	Point3D pointProjection(Point3D p); //TODO : How to obtain the projection point of the point given

};

#endif