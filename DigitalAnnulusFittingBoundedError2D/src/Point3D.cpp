#include "../inc/Point3D.h"

float Point3D::getX(){
	return this->x;
}

float Point3D::getY(){
	return this->y;
}

float Point3D::getZ(){
	return this->z;
}

void Point3D::setX(float _x){
	this->x = _x;
}

void Point3D::setY(float _y){
	this->y = _y;
}

void Point3D::setZ(float _z){
	this->z = _z;
}

void Point3D::translation(float a, float b, float c){
	this->x+=a; this->y+=b; this->z+=c;
}

void Point3D::translation(Point3D t){
	this->x+= t.getX();
	this->y+= t.getY();
	this->z+= t.getZ();
}

void Point3D::inverse(){
	this->x*= -1; this->y*= -1; this->z*= -1;
}

double Point3D::innerProduct(Point3D t){
	return this->x*t.getX()+this->y*t.getY()+this->z*t.getZ();
}

double Point3D::outerProduct(Point3D t){
	return 1;
}

double Point3D::squareDis(){
	double p = this->getX();
	double q = this->getY();
	double r = this->getZ();
	return p*p + q*q + r*r;
}

double Point3D::squareDis(Point3D t){
	double p = this->getX() - t.getX();
	double q = this->getY() - t.getY();
	double r = this->getZ() - t.getZ();
	return p*p + q*q + r*r;
}

void Point3D::permulation(int type){
	if(type < 0 || type > 2) cout<<"error in permulation type";
	else if(type == 1){
		float tem = this->x;
		this->x = this->y;
		this->y = this->z;
		this->z = tem;
	}
	else if(type ==2){
		float tem = this->x;
		this->x = this->z;
		this->z = this->y;
		this->y = tem;
	}

}

bool Point3D::isInside(std::vector<Point3D> ps){
	if ( this->getX()<ps[0].getX() || this->getY()<ps[0].getY() || this->getZ()<ps[0].getZ() 
				|| this->getX()>ps[1].getX() || this->getY()>ps[1].getY() || this->getZ()>ps[1].getZ() ) {
			return false; 
		}
		else {
			return true;
		}
}