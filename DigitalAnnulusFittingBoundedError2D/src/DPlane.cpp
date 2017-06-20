#include "../inc/DPlane.h"

float DPlane::getWidth(){
	return this->w;
}

void DPlane::setWidth(float _w){
	this->w = _w;
}

bool DPlane::contains(Point3D p){
	double distance = p.innerProduct(this->normal)+this->d;
	if( distance >= 0 && distance < this->w )	return true;
	else										return false;
}

bool DPlane::containsWithBoundaryTreatment(Point3D p){
	double distance = p.innerProduct(this->normal)+this->d;
	if( distance >= 0 && distance <= this->w )	return true;
	else										return false;
}

vector<bool> DPlane::getInliers(vector<Point3D> pset){
	int pnum = pset.size();
	vector<bool> inliers (pnum); // initial with all 0 elements
	for(int i = 0; i< pnum; i++)
		if(this->contains(pset[i])) inliers[i] = true;
	return inliers;
}

vector<bool> DPlane::getInliersWithBoundaryTreatment(std::vector<Point3D> pset){
	int pnum = pset.size();
	vector<bool> inliers (pnum); // initial with all 0 elements
	for(int i = 0; i< pnum; i++)
		if(this->containsWithBoundaryTreatment(pset[i])) inliers[i] = true;
	return inliers;
}

int DPlane::getInliersNum(std::vector<Point3D> pset){
	int t = 0;
	for(int i = 0; i < pset.size(); i++)
		if(this->contains(pset[i])) t++;
	return t;

}

int DPlane::getInliersNumWithBoundaryTreatment(std::vector<Point3D> pset){
	int t = 0;
	for(int i = 0; i < pset.size(); i++)
		if(this->containsWithBoundaryTreatment(pset[i])) t++;
	return t;
}

vector<Point3D> DPlane::pointProjection(Point3D p){
	vector<Point3D> proj (2);
	float a = this->normal.innerProduct(p);
	float b = this->normal.squareDis();

	//projection on the first plane
	float t = (-this->d - a)/b;

	float xx = this->normal.getX();
	float x = xx*t + xx;

	float yy = this->normal.getY();
	float y = yy*t + yy;

	float zz = this->normal.getZ();
	float z = zz*t + zz;

	proj[0] = Point3D(x,y,z);

	//projection on the second plane
	t = (this->w - this->d - a)/b;
	x = xx*t + xx;
	y = yy*t + yy;
    z = zz*t+zz;

	proj[1] = Point3D(x,y,z);

	return proj;
}