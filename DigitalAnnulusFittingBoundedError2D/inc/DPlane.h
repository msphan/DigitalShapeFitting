#ifndef DPLANE_H_
#define DPLANE_H_

#include "Plane.h"

class DPlane : public Plane{
private:
	// 0 <= ax + by + cz + d <= w
	float w;
public:
	DPlane(){}
	DPlane(Point3D _normal, float _d, float _w) : Plane(_normal,_d), w(_w) {}
	DPlane(float _a, float _b, float _c, float _d, float _w) : Plane(_a, _b, _c, _d), w(_w) {}
	float getWidth();
	void setWidth(float _w);
	bool contains(Point3D p);
	bool containsWithBoundaryTreatment(Point3D p);
	vector<bool> getInliers(vector<Point3D> pset);
	vector<bool> getInliersWithBoundaryTreatment(vector<Point3D> pset);
	int getInliersNum(vector<Point3D> pset);
	int getInliersNumWithBoundaryTreatment(vector<Point3D> pset);
	vector<Point3D> pointProjection(Point3D p);

};

#endif