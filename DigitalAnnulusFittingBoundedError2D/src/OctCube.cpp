#include "../inc/OctCube.h"

float OctCube::getEpsilonQuery(){
    
    return this->eps_query;
}


int OctCube::round(float n){
	return (n < 0.0) ? ceil(n - 0.5) : floor(n + 0.5);
}

OctCube::OctCube(vector<Point3D> v, float _min_cube_size, vector<Point3D> p){
	
	this->query_table = NULL;
	this->min_cube_size = _min_cube_size;
	this->vertices.push_back(v[0]);
	this->vertices.push_back(v[1]);

	for(int i = 0; i<p.size(); i++) this->points.push_back(p[i]);
	float cube_size = this->getCubeSize();

	// setting epsilon query, it is from Fonsena's paper
	this->eps_query = (this->min_cube_size) /(cube_size);
	
	//initialize query table 
	if(cube_size >= this->min_cube_size){
		
		int slope = this->getQueryPlaneSlopeNum();
		this->slope_num = slope;
		this->query_table = new int** [slope];
		this->intercept_num = new int* [slope];
		for(int i = 0; i<slope; i++){
			this->query_table[i] = new int* [slope];
			this->intercept_num[i] = new int [slope];
			for(int j = 0; j<slope; j++){
				this->query_table[i][j] = NULL;
				int num = this->getQueryPlaneInterceptNum(i,j);
				this->intercept_num[i][j] = num;
				this->query_table[i][j] = new int[num];
				for(int k = 0; k < num; k++)
					this->query_table[i][j][k] = 0;

			}
		}		
	}else{// we have only one element if cube_size > min_cube_size
		this->query_table = new int** [1];
		this->query_table[0] = new int* [1];
		this->query_table[0][0] = new int[1];
		this->query_table[0][0][0] = 0;
		
	}

	float half_cube_size = cube_size/2;

	// make octree subcubes
	if(cube_size >= this->min_cube_size){		
	  
		// get the min and max vertices of 8 subcubes
		vector<Point3D> subvertices[8];
		int u = 0;
		for(int i = 0; i<2; i++)
			for(int j = 0; j<2; j++)
				for(int k = 0; k<2; k++){
					for(int l = 0; l<2; l++){
						subvertices[u].push_back(this->vertices[0]);
						subvertices[u][l].translation(half_cube_size*(i+l),half_cube_size*(j+l),half_cube_size*(k+l));
					}
					u++;
				}
		
		// for storing the set of points contained in each subcube 
		vector<Point3D> subset[8];
		for(int i = 0; i<this->points.size(); i++){
			for(int j = 0; j<8; j++){
				if(this->points[i].isInside(subvertices[j])){
					subset[j].push_back(this->points[i]);
					break;
				}
				
			}
		}
		// make the suboctreecubes if the subcube contains at least one point
		for(int i = 0; i<8; i++)
			if(subset[i].size() > 0){
				this->subcube.push_back(new OctCube(subvertices[i],this->min_cube_size,subset[i]));	
			}
	}

}

// check the intersection between current cube and plane
int OctCube::intersection(Plane p){
	
	vector<Point3D> v = this->getVertices();
	int above = 0;
	int below = 0;
	for(int i = 0; i<8 && above*below == 0; i++){
		float dis = p.signedDis(v[i]);
		if(dis > 0) above++;
		else if(dis < 0) below++;
	}

	if(above==8) return 0; // all the vertices are above the plane
	else if(below==8) return 2; // all the vertices are below the plane 
	else return 1; // some vertices are above the plane while some others are below the plane 	
}

void OctCube::buildQueryTable(){ // using halfspace range searching (see article)
  
	if(this->getCubeSize() < this->min_cube_size){
		this->query_table[0][0][0] = this->getPointNum();
	}
	else{

		// build the query table recursively for the octree subcubes
		for(int i = 0; i < this->subcube.size(); i++)
			this->subcube[i]->buildQueryTable();
		
		// generate all the query planes ax+by+z+c=0
		for(int i = 0; i < this->slope_num; i++)
			for(int j = 0; j < this->slope_num; j++){
			  
				for(int k = 0; k < this->intercept_num[i][j]; k++){
					vector<float> par = this->getQueryPlane(i,j,k);
					Plane p(par[0],par[1],1,par[2]);
					
					for(int l = 0; l < this->subcube.size(); l++){
						int h = this->subcube[l]->queryPlane(p);
						
						this->query_table[i][j][k] += h;
					}
					
			
				}
			}
		
	
	}
}

vector<float> OctCube::getQueryPlane(int i, int j, int k){
	
	vector<float> par (3);
	int pos;
	
	par[0] = -this->eps_query * i;
	par[1] = -this->eps_query * j;

	float _max = 0; 
	float _min = -1;
	if ( par[0]<0 )			// modification depends on the value of a
		_max -= par[0];
	else _min -= par[0];
	if ( par[1]<0 )			// modification depends on the value of b 
		_max -= par[1];
	else _min -= par[1];

	pos = (int)((_max/this->eps_query)+(-_min/this->eps_query)+1)/2; 

	if(k<=pos) par[2] = this->eps_query * k;
	else par[2] = -this->eps_query * (k - pos);

	Point3D t(par[0],par[1],1);
	par[2] = par[2]*this->getCubeSize() - t.innerProduct(this->vertices[0]);

	return par;
}
// Translate plane(i,j,k) to another parallel plane
int OctCube::translationQueryPlane(int i, int j, int k, float t){
	vector<float> par = this->getQueryPlane(i,j,k);
	Point3D p(par[0],par[1],1);
	// translation and scale change 
	float intercept = (par[2] + t + p.innerProduct(this->vertices[0]))/this->getCubeSize(); 
	
	//float intercept = par[2] + t;
	
	// calculate the maximum and minimum values of the intercept
	float _max = 0;
	float _min = -1;
	if(par[0]<0) _max -= par[0];
	else		 _min -= par[0];
	if(par[1]<0) _max -= par[1];
	else		 _min -= par[1];

	/*
	float a = - par[0] * this->getCubeSize();
	float b = - par[1] * this->getCubeSize();
	float temp1 = (a * a + b * b) / ( 2 * this ->getCubeSize());

	if(_max > temp1) _max = temp1; 
	*/
	
	// verify if the new intercept is in the range 
	if ( intercept > _max || intercept < _min ) {
		return -1; 
	}else{
		// obtain the positive part of indexes
		int pos = (int)((_max/this->eps_query)+(-_min/this->eps_query)+1)/2; 
		if(intercept==0) return 0;
		else if(intercept > 0) return this->round(intercept/this->eps_query);
		else return this->round(-intercept/this->eps_query)+pos;
	}
}

float OctCube::getCubeSize(){
	return (this->vertices[1].getX() - this->vertices[0].getX());
}

int OctCube::getPointNum(){
	return this->points.size();
}

vector<Point3D> OctCube::getVertices(){
	vector<Point3D> v (8);
	float _cube_size = this->getCubeSize();
	for(int i = 0; i < 2; i++)
	for(int j = 0; j < 2; j++)
	for(int k = 0; k < 2; k++){
		v[i+j*2+k*4] = Point3D(this->vertices[0]);
		v[i+j*2+k*4].translation(_cube_size*i, _cube_size*j, _cube_size*k);
	}
	return v;
}

// get number of points that is below of the plane
int OctCube::queryPlane(Plane p){
	if(this->getCubeSize() < this->min_cube_size){
		//int d = this->intersection(p); 	
				
		//if ( d==0 ) return 0;
		//else return this->getPointNum();
		return this->getPointNum();
		//return 0;
	}else{
		
		// verify the normal vectors (a,b,1) are in the range -1<=a<=1, -1<=b<=1
		vector<float> par (3);
		par[0] = p.getNormalX();
		par[1] = p.getNormalY();
		if ( par[0]<-1 || par[0]>1 || par[1]<-1 || par[1]>1 )
			return 0;
		// obtain the positive part of indexes for a and b 
		int pos = this->slope_num/2;
		// calculate the two indexes for a and b  
		vector<int> index (3);
		for(int i = 0; i<2; i++){
			
			index[i] = this->round(-par[i]/this->eps_query);
			
			if(index[i] >= this->slope_num) index[i] = this->slope_num-1;
			
			// re-calculate the parameter values of a and b
			par[i] = -this->eps_query*index[i];			
	
		}
		
		// get the intercept c 
		// translation and scale change is necessary for the parameter c  
		par[2] = p.getIntercept() + p.getNormal().innerProduct(this->vertices[0]);
		par[2] = par[2]/this->getCubeSize();

		float max_c = 0;	// initialize the maximum value for c=-z-ax-by (z=0)
		float min_c = -1; 	// initialize the minimum value (z=1)
		
		if(par[0]<0)	 				// modification depends on the value of a
			max_c -= par[0];
		else 
			min_c -= par[0];
		if(par[1]<0)					// modification depends on the value of b 
			max_c -= par[1];
		else 
			min_c -= par[1];

		/*
		float a = - par[0] * this->getCubeSize();
		float b = - par[1] * this->getCubeSize();
		float temp1 = (a * a + b * b) / ( 2 * this ->getCubeSize());

		
		
		if(max_c > temp1) max_c = temp1; 
		*/

		// verify if the parameter is in the rangee 
		if(par[2]>max_c) 						// the cube is below the plane 
			return this->getPointNum();			
		else if(par[2]<min_c) 				// the cube is above the plane 
			return 0; 		
		else{
			// obtain the positive part of indexes
			//pos = (int)(((max_c/this->eps_query)+(int)(min_c/this->eps_query)+1)/2); 	
			pos = (int)((max_c/this->eps_query)+(-min_c/this->eps_query)+1)/2; 
					
			// calculate the index for c 
			index[2] = this->round(par[2]/this->eps_query);
				
			// for zero and positive values, we do not need to change the index except for the one 
			// that is out of the range 
			if(index[2]>pos) 
				return this->getPointNum(); 		// the cube is below the plane 
			// for negative values, 
			else if(index[2]<0){
				index[2] = pos - index[2];		// correct the index 
				if (index[2]>=this->intercept_num[index[0]][index[1]]) 	// the cube is above the plane 
					return 0; 
			}
			//cout<<"slope_num ="<<this->slope_num<<" , inter num = "<<this->intercept_num[index[0]][index[1]]<<"index = "<<index[0]<<" , "<<index[1]<<" , "<<index[2]<<" , query_table= "<<this->query_table[index[0]][index[1]][index[2]]<<endl;
			return this->query_table[index[0]][index[1]][index[2]];
		}

	}

}

// get plane from the index(i,j,k)
int OctCube::queryPlane(int i, int j, int k){
	// verify if the index is adequate 
	if(i<0 || i>=this->slope_num) {
		cout<<"Out of range of the query table.";
		exit(0);
	}
	else if(j<0 || j>=this->slope_num) {
		cout<<"Out of range of the query table.";
		exit(0);
	}
	else if( k<0 || k>=this->intercept_num[i][j]) {
		cout<<"Out of range of the query table.";
		exit(0);
	}
	
	return this->query_table[i][j][k];
}

float OctCube::getQueryPlaneEsp(){
	return this->eps_query;
}

int OctCube::getQueryPlaneSlopeNum(){
	
	if(this->query_table == NULL)
		return (int)ceil(1.0/this->eps_query)*1+1; 		// between -1 and 1 for a and b (between -1 and 1)
	else
		return this->slope_num;

}

int OctCube::getQueryPlaneInterceptNum(int a, int b){
	if(this->query_table[a][b]==NULL){// not yet initialized 	
		// calculate the maximum and minimum values for c, which depend on the values of a and b  
		// obtain the positive part of indexes for a and b
		// calculate the parameter values of a and b 
		float par[2];
		
		par[0] = - this->eps_query * a;
		par[1] = - this->eps_query * b;
		
		// obtain the maximum and minimum values for the parameter c, which depends on the values of a and b 
		float _max = 0; 		// initialize the maximum value for c=-z-ax-by (z=0)
		float _min = -1; 		// initialize the minimum value (z=1)
		if(par[0]<0)	 				// modification depends on the value of a
			_max -= par[0];
		else 
			_min -= par[0];
		if(par[1]<0)					// modification depends on the value of b 
			_max -= par[1];
		else 
			_min -= par[1];

		// get the number of value variations of the parameter c for generating the query planes ax+by+z+c=0
		// Note that Dror's setting is -2<=c<=1.

		// check condition : c <= sqr(a) + sqr(b)

		/*
		float a = - par[0] * this->getCubeSize();
		float b = - par[1] * this->getCubeSize();
		float temp1 = (a * a + b * b) / ( 2 * this ->getCubeSize());

		if(_max > temp1) _max = temp1;  
		*/

		//cout<<"num intercept ="<<(int)((_max/this->eps_query)+(-_min/this->eps_query)+20)<<endl;

		return (int)((_max/this->eps_query)+(-_min/this->eps_query)+1);
	}
	else return this->intercept_num[a][b];
}