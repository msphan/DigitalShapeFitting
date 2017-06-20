#include "../inc/OctCube.h"
#include "../inc/DPlane.h"
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <stdio.h>
#include <ctime>
#include <unistd.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>

#include <GL/glut.h>

float w = 3.0; // default width

// view choice
// view = 1 : annulus view, using approx.method + plane adjusitng
// view = 2 : annulus view, using approx.method
// view = 3 : plane view, using approx.method + plane adjusting
// view = 4 : plane view, using approx.method
int view = 1;

float approx_eps = 3; // default input epsilon

char* fname = 0; // filename of datafile

vector<Point3D> points; // set of points
vector<Point3D> copy_points;


vector<int> plane_in2; 	// inliers and outliers of appox.plane using approx.method
vector<int> plane_in; 	// inliers and outliers of appox.plane using approx.method + plane adjusting
vector<int> plane_io; 	// inliers and outliers of exact plane
Plane pl;	 	// 
Plane copy_pl;

vector<int> circle_in2; // inliers and outliers of approx.annulus using approx.method
vector<int> circle_in; // inliers and outliers of approx.annlus using approx.method + plane adjusting
vector<int> circle_io; // inliers and outliers of exact annulus



int in_points = 400;
int out_points = 100;
int num_inliers = 0;
float IA=100,IB=100,R=90.5; // parameters of true circle; I(IA,IB) and Radius R

// parameters for drawing plane
float xc=0,yc=0,zc=0;
bool pressed=false;
int lx=-1,ly=-1;
float rotx=0.0,roty=0,rotz=0;
float scale=0.08;
float tranx = -90.0, trany = -95.0, tranz = 0.0;

float grid_size = 200;

// annulus parameters
float a,b,rr;
float copy_a,copy_b,copy_rr;

// plane parameters
float ww;
float www;
float wwww;
float ddd;
float c_gamma;
float copy_cgamma;
float c_gam;
float epsilon_query;

// set max min values of data points bounded in 3D grid
void set_max_min(float& x1, float& y1, float& z1, float& x2, float& y2, float& z2, float x, float y, float z){
  
	if (x<x1) x1=x;
	if (y<y1) y1=y;
	if (z<z1) z1=z;
	if (x>x2) x2=x;
	if (y>y2) y2=y;
	if (z>z2) z2=z;
  
}

// get a set of points from file or randomize
void get_points(){
	
	
	int x,y;
	float z;
	
	if (fname){
	  
		// define the max min values
		float _x1=100000,_y1=100000,_z1=100000;
		float _x2=-100000,_y2=-100000,_z2=-100000;
	
		fstream fs(fname,ios::in);
		if (!fs || fs.fail()){
			perror("can't read from input file");
			exit(-1);
		}
		
		// read grid size (the first line)
		fs>>grid_size;
		
		// and the next lines is coordinates parameters
		while (!fs.eof()){
			fs>>x>>y;
			if (!fs.eof()){
				
				z = (x * x + y * y) / (float) ( 2 * grid_size );
				points.push_back(Point3D(x,y,z));
				set_max_min(_x1,_y1,_z1,_x2,_y2,_z2,x,y,z);
			}
		}
	
		// translation to original coordinates (0,0,0)
		for (int i=0;i<points.size();i++){
			points[i].setX(points[i].getX()-_x1);
			points[i].setY(points[i].getY()-_y1);
			points[i].setZ(points[i].getZ()-_z1);
		}
	  
	  
	}else{
	 
	
		// get random points
		srand(1283779988);
		points.resize(0);
		
		// generate set of points for circle
		
		// points with radius smaller than R 0.5
		R -= 0.5;
		bool check = true;		

		while (points.size()<(in_points/2)){ 

			// get random x
			x=(IA - R) + rand() % (int) (2 * R + 1);
			
			// check for redundance, only one coordinates 
			for(int i = 0; i< points.size(); i++)
				if(points[i].getX() == x){
					check = false;
					break;
				}

			if (check == false){
				check = true;
				continue;
			}
			
			// calculate y and z
			float temp = R * R - (x - IA)*(x - IA);
			y = (int)((sqrt(temp) + IB));
			
			z = (x * x + y * y) / (float) ( 2 * grid_size );
			points.push_back(Point3D(x,y,z));
			
			num_inliers ++; 

			// for one x we have two results of y = (+- sqrt(r*r - (x-a)*(x-a))) + b
			if(temp != 0 && points.size()<(in_points/2)){
				y = (int)((-sqrt(temp) + IB));
			
				z = (x * x + y * y) / (float) ( 2 * grid_size );
				points.push_back(Point3D(x,y,z));
				
				num_inliers ++; 

			}

			
		}
		
		// points with radius greater than R 0.5
		R += 1;
		check = true;
		while (points.size()<in_points){ 

			// get random x
			x=(IA - R) + rand() % (int)(2 * R + 1);
			
			for(int i = in_points/2; i< points.size(); i++)
				if(points[i].getX() == x){
					check = false;
					break;
				}

			if (check == false){
				check = true;
				continue;
			}
			
			// else calculate y
			float temp = R * R - (x - IA)*(x - IA);
			y = (int)((sqrt(temp) + IB));
			
			z = (x * x + y * y) / (float) ( 2 * grid_size );
			points.push_back(Point3D(x,y,z));
			num_inliers ++;
			
			
			if(temp != 0 && points.size()<in_points){
				y = (int)((-sqrt(temp) + IB));
				
				z = (x * x + y * y) / (float) ( 2 * grid_size );
				points.push_back(Point3D(x,y,z));
				num_inliers ++;
				

			}

		}
		
		R -= 0.5;
		
		check = true;
		
		// generate outliers
		while(points.size() < (in_points + out_points)){
			x = rand() % (int)grid_size;
			y = rand() % (int)grid_size;

			for(int j = 0; j< points.size(); j++)
				if(points[j].getX() == x && points[j].getY() == y) {
					check = false;
					break;
				}	

			if (check == false){
				check = true;
				continue;
			}

			float temp = (x - IA) * (x - IA) + (y - IB) * (y - IB);

			// check if inliers
			if(temp < (R - 3)*(R - 3) || temp > (R + 3)*(R + 3)){
			  
			  z = (x * x + y * y) / (float) ( 2 * grid_size );

			  points.push_back(Point3D(x,y,z));
			  
			}

			
		}
		
		grid_size = 200;
		cout<<"IA, IB, R2="<<IA<<","<<IB<<","<<R<<endl;
	}

}

// get values from input arguments
void get_args(int argc, char **argv){
	int i = 1;
	while(i < argc) {
		if (!strcmp(argv[i], "-h")){
			printf("Usage: %s [-f file -w plane_width -e error_epsilon -v view]\n",argv[0]);
			printf("file: data points file directory\n");
			printf("default plane_width equal 3.0\n");
			printf("default error epsilon equal 2.0\n");
			printf("default view equal 1 (annulus view after using approx.method + plane adjusting)\n");
			printf("--view equal 2 (annulus view after using approx.method)\n");
			printf("--view equal 3 (plane view after using approx.method + plane adjusting)\n");
			printf("--view equal 4 (plane view after using approx.method)\n");
			exit(-1);
		}
		else if (!strcmp(argv[i], "-f")){
			fname=argv[++i];
		}
		else if (!strcmp(argv[i], "-w")){
			w=atof(argv[++i]);
		}
		else if(!strcmp(argv[i], "-e")){
			approx_eps=atof(argv[++i]);
		}
		else if(!strcmp(argv[i], "-v")){
			view=atof(argv[++i]);
		}
		else if(argv[i][0] == '-'){
			cerr << "Unrecognized option.(-" << argv[i] << ")\n";
			exit(-1);
		}
		i++;
	}
}

void draw_circle(int a, int b, float r){
    
    glBegin(GL_LINES);
    float x,y;
    float PI = 3.14159;
    for (int i = 0; i < 180; i++)
    {
      float Angle = i * (2.0*PI/180);
      float X = cos( Angle )*r;
      float Y = sin( Angle )*r;
      X += a;
      Y += b;
      glVertex2f( X, Y );

    }
    glEnd();
  
}

void draw_plane(float A,float B,float C,float D){
  
    glBegin(GL_POLYGON);
    
    float xc=grid_size/2;		
    float yc=grid_size/2;
    float zc=grid_size/2;
    
    float x,y,z;
	
    glColor3f(1.0, 0.0, 1.0);
    x=0;
    y=0;
    z= (-A*x-B*y-D)/C;
    glVertex3f(x,y,z);
    
    x=xc*2;
    y=0;
    z= (-A*x-B*y-D)/C;
    glVertex3f(x,y,z);
    
    x=xc*2;
    y=yc*2;
    z= (-A*x-B*y-D)/C;
    glVertex3f(x,y,z);

    x=0;
    y=yc*2;
    z= (-A*x-B*y-D)/C;
    glVertex3f(x,y,z);

    glEnd();
}

void plane_view(int cases){
  
    glBegin(GL_POINTS);
    
    vector<int> plane_temp;
    Plane p_temp;
    float plane_width;
    
    // display results if using approx.method + plane adjusting
    if(cases == 1){
      
      p_temp = pl;
      plane_width = wwww;
      plane_temp.resize(plane_in.size());
      for(int i = 0; i< plane_in.size(); i++)
	plane_temp[i] = plane_in[i];
    }
    
    // display results if using only approx.method
    else if(cases == 2){
      
      p_temp = copy_pl;
      plane_width = copy_cgamma;
      plane_temp.resize(plane_in2.size());
      for(int i = 0; i< plane_in2.size(); i++)
	plane_temp[i] = plane_in2[i];
    }
	
    for (int i=0;i<copy_points.size();i++)
    {
	    // draw points generated in 3D
	
	    if (!fname && plane_io[i]) glColor3f(0.0,0.0,1.0);
	    else glColor3f(0.0,0.0,0.0);
	    
	    
	    glVertex3f(copy_points[i].getX(), copy_points[i].getY(), copy_points[i].getZ());
	
    }
    
    glEnd();
    
   
    // draw a approximate plane
    draw_plane(p_temp.getNormalX(),p_temp.getNormalY(),p_temp.getNormalZ(),p_temp.getIntercept()-(plane_width/2.0));
    draw_plane(p_temp.getNormalX(),p_temp.getNormalY(),p_temp.getNormalZ(),p_temp.getIntercept()+(plane_width/2.0));
    
    glPopMatrix();
    glutSwapBuffers();

}

void circle_view(int cases){
  
  
    glBegin(GL_POINTS);
    
    vector<int> circle_temp;
    float circle_width;
    float ca,cb,cr;
    
    circle_temp.resize(copy_points.size());
    
    // display results if using only approx.method
    if(cases == 1){
	for(int i = 0; i<circle_in2.size();i++)
	  circle_temp[i] = circle_in2[i];
	circle_width = c_gamma;
	ca = copy_a;
	cb = copy_b;
	cr = copy_rr;
      
    }
    
    // display results if using approx.method + plane adjusting 
    else if(cases == 2){
	for(int i = 0; i<circle_in.size();i++)
	  circle_temp[i] = circle_in[i];
	circle_width = ww;
	ca = a;
	cb = b;
	cr = rr;
    }
      
      
    for (int i=0;i<copy_points.size();i++)
    {
	    // draw original points in 2D (circle)
	    
	    if (!fname && circle_io[i]) glColor3f(0.0,0.0,1.0);
	    else glColor3f(0.0,0.0,0.0);
	    
	    glVertex3f(copy_points[i].getX(), copy_points[i].getY(), 0);
	    
    }
    
    glEnd();
    
    if (!fname){
      
      // draw a true circle
      glLineWidth(1.0);
      glColor3f(0.0, 0.0, 1.0);
      draw_circle(IA,IB,R - w/2.0);
      draw_circle(IA,IB,R + w/2.0);
      
    }
    
    
    // draw a appox.circle
    glLineWidth(2.0);
    glColor3f(1.0, 0.0, 0.0);
    draw_circle(ca,cb,cr - circle_width/2.0);
    draw_circle(ca,cb,cr + circle_width/2.0);    
    
    glPopMatrix();
    glutSwapBuffers();
   
}

void draw(void){
  
	GLdouble red, green, blue;
	int i;

	glClearColor( 1.0f, 1.0f, 1.0f, 0.5f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glClear(GL_DEPTH_BUFFER_BIT);
	glDisable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);

	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_FLAT);

	glPushMatrix();
	glMatrixMode(GL_MODELVIEW);
	glTranslatef(0,0,-100);
	glScalef(scale,scale,scale);
	glRotatef(rotx,0,0,1);
	glRotatef(roty,0,1,0);
	glRotatef(rotz,1,0,0);
	glTranslatef(tranx,trany,tranz);
	glPointSize(3.0);
	
	
	// Draw points inliers and outliers
	
	if(view == 3) plane_view(1);
	else if(view == 4) plane_view(2);
	else if(view == 2) circle_view(1);
	else circle_view(2);
}

void keyboard(unsigned char ch, int x, int y){
	
	switch (ch) {
case 'a':
	rotx+=1.0;
	break;
case 'A':
	rotx-=1.0;
	break;
case 'z':
	roty+=1.0;
	break;
case 'Z':
	roty-=1.0;
	break;
case 'e':
	rotz+=1.0;
	break;
case 'E':
	rotz-=1.0;
	break;
case 'q':
	scale *= 1.2;
	break;
case 'Q':
	scale *= 0.8;
	break;
case 's':
	tranx += 2.0;
	break;	
case 'S':
	tranx -= 2.0;
	break;
case 'd':
	trany += 2.0;
	break;	
case 'D':
	trany -= 2.0;
	break;	
case 'f':
	tranz += 2.0;
	break;	
case 'F':
	tranz -= 2.0;
	break;	
case 27:             /* ESC */
	exit(0);
	break;
	}
	draw();
}

void processMouseActiveMotion(int x, int y){
	if (lx!=-1 && pressed){
		roty += (x-lx)/2;
		rotz += -(y-ly)/2;
	}
	lx=x;
	ly=y;
	glutPostRedisplay();
}

void processMousePassiveMotion(int x, int y){
}

void processMouseEntry(int state){
	if (state == GLUT_LEFT)	{
	}
	else{
	}
}

void processMouse(int button, int state, int x, int y){
	if (state == GLUT_DOWN)	{
	  
		if (button == GLUT_LEFT_BUTTON){
			lx=x;
			ly=y;
			pressed=true;
		}
	} else{
		if (button == GLUT_LEFT_BUTTON) {
			pressed=false;
		}
	}
}

// count the n.of.inliers of cirlce
int circle_count(int a, int b, float r, float width){
  
    int t = 0;
    int x,y;
    float temp;
    for(int i = 0; i< copy_points.size(); i++){
     
      x = copy_points[i].getX();
      y = copy_points[i].getY();
      
      temp = (x - a) * (x - a) + (y - b) * (y - b);
      
      if(temp >= (r - width) * (r - width) &&
	 temp <=(r + width) * (r + width))
	t++;
    }
    return t;
}

// count the n.of.inliers of cirlce, results obtained is saved in vector re
int circle_count_report(int a, int b, float r, float width, vector<int>& re){
  
    int t = 0;
    int x,y;
    float temp;
    re.resize(copy_points.size());
    for(int i = 0; i< copy_points.size(); i++){
     
      x = copy_points[i].getX();
      y = copy_points[i].getY();
      
      temp = (x - a) * (x - a) + (y - b) * (y - b);
      
      if(temp >= (r - width) * (r - width) &&
	 temp <=(r + width) * (r + width)){
	    re[i] = 1;
	    t++;
	}
	else re[i] = 0;
	
    }
    return t;
}

// check condition of C parameter and compute "W" (see article)
float checkC(float A, float B, float C){
  
    float a = - A * grid_size;
    float b = - B * grid_size;
    float r = (a * a + b * b) - ( C * 2 * grid_size);
    if(r < 0) return 0;
    else{
      
      r = r + float(w /2.0);
      float gam = sqrt(r) * (w) / grid_size;
            
      return gam; 
    }
}

int main(int argc,char *argv[]){

	// obtain the arguments from io stream (screen)
	get_args(argc,argv);

	// get data points from file, otherwise generate randomly data points 
	get_points();
	
	cout<<"grid_size = "<<grid_size<<endl;
	
	if(!fname) cout<<"exact.n.of.inliers (counted by randomizing data points) = "<<num_inliers<<endl;
	cout<<"epsilon:"<<approx_eps<<endl;
	
	// save 2D data (annulus) to file
	fstream f2("circle_data_400_100.txt",ios::out);
	// wirte grid size
	f2<<grid_size<<endl;
	
	// write data points
	for (int i=0; i < points.size(); i++)
	{
		f2<<points[i].getX()<<" "<<points[i].getY()<<endl;
	}
	f2.close();
	
	// copy data points before using halfspace range searching
	copy_points.resize(points.size());
	for(int i = 0; i< points.size(); i++)
	  copy_points[i] = points[i];
	
	// count the number of inliers of exact digital annulus by using circle_count_report() 	  
	if(!fname){
	  
	    cout<<"exact.n.of.inliers (counted by circle function) = "<<circle_count_report(IA,IB,R,w/2.0,circle_io)<<endl;
	
	    // convert the exact digital annulus to the exact digital plane
	    float AA = -IA/grid_size;
	    float BB = -IB/grid_size;
	    float CC = ( IA * IA + IB * IB - (R) * (R) ) / ( 2 * grid_size ); // middle intercept
	    www = w * R / grid_size;
	    
	    float CC1 = ( IA * IA + IB * IB - (R - w/2.0) * (R - w/2.0) ) / ( 2 * grid_size ); // lower-bound intercept
	    float CC2 = ( IA * IA + IB * IB - (R + w/2.0) * (R + w/2.0) ) / ( 2 * grid_size ); // upper-bound intercept
	    
	    cout<<"exact plane :"<<endl;
	    cout<<"A, B, C, C'= "<<AA<<" , "<<BB<<" , "<<CC1<<" , "<<CC2<<endl;
	    Plane ep(AA,BB,1,CC);
	    
	    // count the number of inliers of exact 
	    cout<<"exact.n.of.inliers (counted by plane function) = "<<ep.count_report(points,www/2.0,plane_io)<<endl;
	  
	} 
	
	/// APPROXIMATE METHOD FOR DIGITAL ANNULUS FITTING (SEE ARTICLE)
	
	int n_opt = 0;
	DPlane p_opt;

	// initialize time counter
	timeval t1;
	timeval t2;
	gettimeofday(&t1, 0);

	// copy points
	copy_points.resize(points.size());
	for(int i = 0; i< points.size(); i++)
	  copy_points[i] = points[i];

	
	cout<<"number of points :"<<points.size()<<endl;
	vector<Point3D> vertices (2);
	vertices[0] = Point3D(0,0,0);
	vertices[1] = Point3D(grid_size,grid_size,grid_size);

	// make an octree cube 
	OctCube* octcube = new OctCube(vertices, approx_eps, points);	
	
	// build the query table for the octree cube
	octcube->buildQueryTable();
	cout<<"build ok..."<<endl;
	cout<<"-------------------------------------------------------------------"<<endl;
	
	// get all the slope indices 
	int slope = octcube->getQueryPlaneSlopeNum(); 

	epsilon_query = octcube->getEpsilonQuery();
	cout<<"epsilon query = "<<epsilon_query<<endl;

	// width of the approximate digital plane
	int gamma ;
	
	for(int i=0; i<slope; i++)
		for(int j=0; j<slope; j++){
			
			// get the intercept indices 
			int intercept = octcube->getQueryPlaneInterceptNum(i,j);
			
			// for each query plane 
			for(int k=0; k<intercept; k++){
				
				vector<float> par = octcube->getQueryPlane(i,j,k);
				
				// check 'C' conditions and convert "w" to "W(D)". see papers
				float gam = checkC(par[0],par[1],par[2]);
				
				if(gam < 0) continue; 		
				
				// the efficient width for querying
				gamma = gam + 5 * approx_eps;
				
				// get query result from the lower plane with index (i,j,k) 
				int n1 = octcube->queryPlane(i, j, k);
				
				// get next index of the upper plane 		
				int new_index = octcube->translationQueryPlane(i, j, k, gamma);
				
				// get parameters of the upper plane
				vector<float> newpar = octcube->getQueryPlane(i,j,new_index);
				
				if(!checkC(newpar[0],newpar[1],newpar[2])) continue;
				
				if(new_index>=0 && new_index<intercept){
			
					// get query result from the upper plane
					int n2 = octcube->queryPlane(i, j, new_index);
					int n = n2 - n1;
					
					// find the best digital plane
					if(n>n_opt) {
						
						n_opt = n;
						
						DPlane temp(par[0], par[1], 1, par[2], gamma);
						p_opt = temp;	//set the support plane of the best digital plane 
						
						c_gamma = gamma; // mark the width
						c_gam = gam; // mark the converting of "w"
					
					
					}
				}
			}
		}

	cout<<"query ok..."<<endl;
	cout<<"-------------------------------------------------------------------"<<endl;

	
	// get the middle intercept of p_opt
	float dl = p_opt.getIntercept() + (c_gamma / 2.0);
	
	int nl = 0;
	
	float A = p_opt.getNormalX();
	float B = p_opt.getNormalY() ;
	float C = dl ;
	
	// set the middle plane
	Plane pl0(A,B,1,C);
	
	// convert to annulus
	a = -A * grid_size;
	b = -B * grid_size;
	rr = sqrt(a * a + b * b - C * 2 * grid_size);

	// copy parameters (for drawing)
	copy_a = a;
	copy_b = b;
	copy_rr = rr;
	
	cout<<"Approximate method for digital annulus fitting :"<<endl;
	cout<<"Plane :"<<endl;
	cout<<"Width (W + 5 approx_eps) = "<<c_gamma<<endl;
	cout<<"A = "<<p_opt.getNormalX()<<" , B = "<<p_opt.getNormalY()<<" , C' = "<<p_opt.getIntercept()<<" , C = "<<p_opt.getIntercept() + c_gamma<<endl;
	cout<<"n.of.inliers (width = W) ="<<pl0.count(copy_points,c_gam/2.0)<<endl;
	cout<<"n.of.inliers (width = W + 5 approx_eps) = "<<pl0.count_report(copy_points,c_gamma/2.0,plane_in2)<<endl;
	cout<<"n.of.inliers (n_opt) = "<<n_opt<<endl;
	
	// copy for drawing
	copy_cgamma = c_gamma;
	
	// copy plane (for drawing)
	copy_pl = pl0;
	
	// convert width of digital plane to width of digital annulus
	c_gamma = c_gamma * grid_size / rr;
	
	cout<<"Annulus :"<<endl;
	cout<<"Width (W + 5 approx_eps) = "<<c_gamma<<endl;
	cout<<"r = "<<rr<<", a = "<<a<<", b = "<<b<<endl;
	cout<<"n.of.inliers (width = w) ="<<circle_count(a,b,rr,w/2.0)<<endl;
	cout<<"n.of.inliers (width = (W + 5 approx_eps) * grid_size / r ) = "<<circle_count_report(a,b,rr,c_gamma/2.0,circle_in2)<<endl;
	
	// get runtime
	gettimeofday(&t2, 0);
	long seconds  = t2.tv_sec  -t1.tv_sec;
	long useconds = t2.tv_usec - t1.tv_usec;
	cout<<"Time is:"<<((seconds) * 1000 + useconds/1000.0) + 0.5<<endl;
	
	cout<<"-------------------------------------------------------------------"<<endl;
		
	/// OPTIMIZING USING PLANE ADJUSTING
	
	// FIRST OPTIMIZING (rotate slope and shift intercept parameters)
	
	// optimizing slope interval
	float ocenter = 5 * approx_eps;
	
	float aaa = 0;
	float bbb = 0;
		
	for(float h = -ocenter; h <= ocenter; h+=1)
	  for(float k = -ocenter; k<= ocenter; k+=1){
	    
	    for(float i = - c_gamma ; i <= c_gamma; i+= 1){
	  
	      // check condition of "c" parameters
	      if((rr + i - w > 0 && rr + i + w < grid_size) && (a + h > 0 && a + h < grid_size) && (b + k > 0 && b + k < grid_size)){
		
		int ntemp = circle_count(a+h,b+k,rr + i,w/2);
		
		if( ntemp > nl){
		  nl = ntemp;
		  ddd = i;
		  aaa = h;
		  bbb = k;
		}
		
	      }

	    }
	  }
	
	// adjust the parameters
	a += aaa;
	b += bbb;
	rr += ddd;
	
	// make a new width of approximate digital annulus
	ww = w + 3 * approx_eps;
	
	// convert annulus to plane
	float A2 = -a/grid_size;
	float B2 = -b/grid_size;
	float C2 = ( a * a + b * b - rr * rr ) / ( 2 * grid_size );	
	
	www = w * rr / grid_size;
	wwww = ww * rr / grid_size;
	
	Plane pl1(A2,B2,1,C2);
	pl = pl1;
	
	/*
	
	cout<<"The first optimizing :"<<endl;
	cout<<"Plane :"<<endl;
	cout<<" A = "<<A2<<", B = "<<B2<<", C = "<<C2<<endl;
	cout<<"n.of.inliers (width = W) ="<<pl.count(copy_points,www/2.0)<<endl;
	cout<<"n.of.inliers (width = (w + 3 * eps) * r / grid_size) = "<<pl.count(copy_points,wwww/2.0)<<endl;
	
	cout<<"Annulus :"<<endl;
	
	cout<<"r = "<<rr<<", a = "<<a<<", b = "<<b<<endl;
	cout<<"n.of.inliers (width = w) ="<<circle_count(a,b,rr,w/2.0)<<endl;
	cout<<"n.of.inliers (width = w + 3 * eps) = "<<circle_count(a,b,rr,ww/2.0)<<endl;
	
	// get runtime
	gettimeofday(&t2, 0);
	seconds  = t2.tv_sec  -t1.tv_sec;
	useconds = t2.tv_usec - t1.tv_usec;
	cout<<"Time is:"<<((seconds) * 1000 + useconds/1000.0) + 0.5<<endl;
	
	cout<<"-------------------------------------------------------------------"<<endl;
	
	*/

	// SECOND OPTIMIZING (shift r along to width (w + 3 * approx_eps)
	
	float iii = 0;
	float n_circle = circle_count(a,b,rr,w/2.0);
	
	// save to file

	for(float i = w/2.0; i<= ww - w /2.0; i+=0.1){
	  
	  int n_temp = circle_count(a,b,((rr+ww/2.0) - i),w/2.0);
	
	  if(n_temp > n_circle){
	    
	      n_circle = n_temp;
	      iii = i;    
	  }
	  
	}
	
	// adjust radius
	if(iii!=0)
	  rr = ((rr+ww/2.0) - iii);
	
	// get runtime
	gettimeofday(&t2, 0);
	seconds  = t2.tv_sec  -t1.tv_sec;
	useconds = t2.tv_usec - t1.tv_usec;
	
	// recompute plane parameter
	C2 = ( a * a + b * b - rr * rr ) / ( 2 * grid_size );	
	
	www = w * rr / grid_size;
	wwww = ww * rr / grid_size;
	
	Plane pl2(A2,B2,1,C2);
	pl = pl2;
	
	cout<<"optimizing using plane adjusting :"<<endl;
	cout<<"Plane :"<<endl;
	cout<<"A = "<<A2<<", B = "<<B2<<", C = "<<C2<<endl;
	cout<<"n.of.inliers (width = W) ="<<pl.count(copy_points,www/2.0)<<endl;
	cout<<"n.of.inliers (width = (w + 3 * eps) * r / grid_size) = "<<pl.count_report(copy_points,wwww/2.0,plane_in)<<endl;
	cout<<"Annulus :"<<endl;
	cout<<"r = "<<rr<<", a = "<<a<<", b = "<<b<<endl;
	cout<<"n.of.inliers (width = w) ="<<circle_count(a,b,rr,w/2.0)<<endl;
	cout<<"n.of.inliers (width = w + 3 * eps) = "<<circle_count_report(a,b,rr,ww/2.0,circle_in)<<endl;
	cout<<"Time is:"<<((seconds) * 1000 + useconds/1000.0) + 0.5<<endl;
	cout<<"-------------------------------------------------------------------"<<endl;
	
	// display results using openGL .
	
	if(view == 3 || view == 4){
	  
	  roty = -44, rotx = -50, rotz = -95;
	  scale = 0.06;
	  tranx = -165.0, trany = -35.0, tranz = -50.0;
	}
	
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);
	glutInitWindowSize(600,600);
	glutCreateWindow("Digital annulus fitting");
	glutDisplayFunc(draw);
	glClearDepth(1.0);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	glOrtho(-10.0, 10.0, -10.0, 10.0, 1, 1000);	
	glMatrixMode(GL_MODELVIEW);

	glutKeyboardFunc(keyboard);
	glutMouseFunc(processMouse);
	glutMotionFunc(processMouseActiveMotion);
	glutPassiveMotionFunc(processMousePassiveMotion);
	glutEntryFunc(processMouseEntry);

	glutMainLoop();

	
	return 0;
}