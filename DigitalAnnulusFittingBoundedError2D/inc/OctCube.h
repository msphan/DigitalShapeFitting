/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 200x by <your name/ organization here>
 */
#ifndef OCTCUBE_H_
#define OCTCUBE_H_

#include "Plane.h"

class OctCube{
private:
	float min_cube_size;		// the minimum cube side length 
	vector<Point3D> vertices;	// the minimum and maximum vertices of the current cube 
	vector<Point3D> points;		// the set of points in the current cube  
	vector<OctCube*> subcube;	// the childrens' octree cube (maximum 8 children)
	
	float eps_query;		// the parameter value interval for making the query planes, which are defined by   
					// ax+ by+ z+ c =0 where -1 <= a,b <=1, and min_c(a,b)<= c <= max_c(a,b) (this is different 
					// from Dror's setting: -2 <= c <= 1) for the normalized image size (the image size is 
					// considered to be 1).  
					// The parameter values are multiples of "eps_query". 
					// Note that this value is different for each level of the octree cube. 
	
	int*** query_table;			// For each query plane with the value variations of a, b, c, 
	int slope_num;	// a,b interval
	int** intercept_num; // intercept interval 
	

public:
	static int round(float n);
	OctCube(){}
	OctCube(vector<Point3D> v, float _min_cube_size, vector<Point3D> p);
	
	
	float getEpsilonQuery();
	
	/*!
	 * \brief
	 * Intersection with a plane
	 * 
	 * \param p
	 * Plane object
	 * 
	 * \returns
	 * 0 : the cube is above the plane
	 * 1 : the cube intersect with the plane
	 * 2 : the cube is below the plane
	 */
	int intersection(Plane p);

	/*!
	 * \brief
	 * Build the data structure for the halfspace range counting of the query planes (with approximation).
	 */
	void buildQueryTable();
	
	/*!
	 * \brief
	 * Obtain the query plane parameters (a,b,c) from an index (i,j,k).
	 * 
	 * \param i
	 * Index of plane.
	 * 
	 * \param j
	 * Index of plane.
	 * 
	 * \param k
	 * Index of plane.
	 * 
	 * \returns
	 * Plane parameters (a,b,c).
	 */
	vector<float> getQueryPlane(int i, int j, int k);
	
	/*!
	 * \brief
	 * Obtain the intercept index of the most approximated plane after translation ax + by + z + c + t = 0.
	 * 
	 * \param i
	 * Description of parameter i.
	 * 
	 * \param j
	 * Description of parameter j.
	 * 
	 * \param k
	 * Description of parameter k.
	 * 
	 * \returns
	 * Write description of return value here.
	 *
	 * \remarks
	 * Write remarks for translationQueryPlane here.
	 */
	int translationQueryPlane(int i, int j, int k, float t);

	float getCubeSize();
	int getPointNum();

	/*!
	 * \brief
	 * Get the eight cube vertices.
	 * 
	 * \returns
	 * Write description of return value here.
	 * 
	 */
	vector<Point3D> getVertices();

	/*!
	 * \brief
	 * Query a plane.
	 * 
	 * \param p
	 * Description of parameter p.
	 * 
	 * \returns
	 * Write description of return value here.
	 * 
	 * \throws <exception class>
	 * Description of criteria for throwing this exception.
	 * 
	 * Write detailed description for queryPlane here.
	 * 
	 * \remarks
	 * Write remarks for queryPlane here.
	 * 
	 * \see
	 * Separate items with the '|' character.
	 */
	int queryPlane(Plane p);
	
	/*!
	 * \brief
	 * Write brief comment for queryPlane here.
	 * 
	 * \param i
	 * Description of parameter i.
	 * 
	 * \param j
	 * Description of parameter j.
	 * 
	 * \param k
	 * Description of parameter k.
	 * 
	 * \returns
	 * Write description of return value here.
	 * 
	 * \throws <exception class>
	 * Description of criteria for throwing this exception.
	 * 
	 * Write detailed description for queryPlane here.
	 * 
	 * \remarks
	 * Write remarks for queryPlane here.
	 * 
	 * \see
	 * Separate items with the '|' character.
	 */
	int queryPlane(int i, int j, int k);

	float getQueryPlaneEsp();
	int getQueryPlaneSlopeNum();

	/*!
	 * \brief
	 * Write brief comment for getQueryPlaneInterceptNum here.
	 * 
	 * \param a
	 * Description of parameter a.
	 * 
	 * \param b
	 * Description of parameter b.
	 * 
	 * \returns
	 * Write description of return value here.
	 * 
	 * \throws <exception class>
	 * Description of criteria for throwing this exception.
	 * 
	 * Write detailed description for getQueryPlaneInterceptNum here.
	 * 
	 * \remarks
	 * Write remarks for getQueryPlaneInterceptNum here.
	 * 
	 * \see
	 * Separate items with the '|' character.
	 */
	int getQueryPlaneInterceptNum(int a, int b);

};

#endif