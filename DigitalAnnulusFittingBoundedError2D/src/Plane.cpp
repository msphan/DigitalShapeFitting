#include "../inc/Plane.h"

float Plane::getNormalX(){
	return this->normal.getX();
}

float Plane::getNormalY(){
	return this->normal.getY();
}

float Plane::getNormalZ(){
	return this->normal.getZ();
}

Point3D Plane::getNormal(){
	return this->normal;
}

float Plane::getIntercept(){
	return this->d;
}

void Plane::setNormalX(float _x){
	this->normal.setX(_x);
}

void Plane::setNormalY(float _y){
	this->normal.setY(_y);
}

void Plane::setNormalZ(float _z){
	this->normal.setZ(_z);
}

void Plane::setNormal(Point3D _normal){
	this->normal = _normal;
}

void Plane::setNormal(float _x, float _y, float _z){
	this->normal = Point3D(_x,_y,_z);
}

void Plane::setIntercept(float _d){
	this->d = _d;
}

double Plane::signedDis(Point3D p){
	return (this->normal.innerProduct(p) + this->d);
}

int Plane::principalProjectionPlane(){
	float a = abs((int)this->normal.getX());
	float b = abs((int)this->normal.getY());
	float c = abs((int)this->normal.getZ());

	if ( c >= a && c >= b )			return 1;	// xy-plane
    else if ( a >= b && a >= c )	return 2; 	// yz-plane
    else				    		return 3;	// zx-plane 
}

Point3D Plane::pointProjection(Point3D p){
	float a = this->normal.innerProduct(p);
	float b = this->normal.squareDis();

	float t = (-this->d - a)/b;
	
	float xx = this->normal.getX();
	float x = xx*t + xx;

	float yy = this->normal.getY();
	float y = yy*t + yy;

	float zz = this->normal.getZ();
	float z = zz*t + zz;

	Point3D proj(x,y,z);

	return proj;

}

int Plane::count(vector<Point3D> ps, float dd){
  
	int t = 0;
	for(int i = 0; i < ps.size();i++){
	  
	  if(fabs(this->signedDis(ps[i])) <= dd) t++;
	  
	}
	  
	  
	return t;
}

int Plane::count_report(vector<Point3D> ps, float dd, vector<int>& in){
  
	int t = 0;
	in.resize(ps.size());
	for(int i = 0; i<ps.size();i++)
	  if(fabs(this->signedDis(ps[i])) <= dd){
	    t++;
	    in[i] = 1;
	    
	  } 
	  else in[i] = 0;
	  
	return t;
}
