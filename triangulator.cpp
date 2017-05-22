#define triREAL double
/*#define radius 1.0
#define edgelength 0.02*/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

#include "triangle.h"
#include "mpi.h"

using namespace std;

double radius = 0.0;
double pii = 2.0 * asin(1.0);


void convert_xl(double x, double y, double z, double* lat, double* lon);

class pnt {
public:
	double x, y, z;

	pnt(double x_, double y_, double z_)
		: x(x_), y(y_), z(z_) { }

	pnt()
		: x(0.0), y(0.0), z(0.0) { }

	pnt& operator=(const pnt &p){
		x = p.x;
		y = p.y;
		z = p.z;
		return *this;
	}

	bool operator==(const pnt &p) const {
		return (x == p.x) && (y == p.y) && (z == p.z);
	}

	pnt operator-(const pnt &p) const {
		double x_, y_, z_;

		x_ = x-p.x;
		y_ = y-p.y;
		z_ = z-p.z;

		return pnt(x_,y_,z_);
	}

	pnt operator+(const pnt &p) const {
		double x_, y_, z_;

		x_ = x+p.x;
		y_ = y+p.y;
		z_ = z+p.z;

		return pnt(x_,y_,z_);
	}

	pnt operator*(double d) const {
		double x_, y_, z_;

		x_ = x*d;
		y_ = y*d;
		z_ = z*d;

		return pnt(x_,y_,z_);
	}

	pnt operator/(double d) const {
		double x_, y_, z_;

		if(d == 0.0){
			
		}

		x_ = x/d;
		y_ = y/d;
		z_ = z/d;
		return pnt(x_,y_,z_);
	}

	pnt& operator/=(double d){
		if(d == 0.0){
			
		}
		x = x/d;
		y = y/d;
		z = z/d;
		return *this;
	}

	pnt& operator+=(const pnt &p){
		x += p.x;
		y += p.y;
		z += p.z;
		return *this;
	}

	void normalize(){
		double norm;

		norm = x*x + y*y + z*z;
		norm = sqrt(norm);

		x = x/norm;
		y = y/norm;
		z = z/norm;
	}

	double dot(const pnt &p) const {
		double junk;
		junk = x*p.x+y*p.y+z*p.z;
		return junk;
	}

	pnt cross(const pnt &p) const {
		double x_, y_, z_;

		x_ = y*p.z - p.y*z;
		y_ = z*p.x - p.z*x;
		z_ = x*p.y - p.x*y;

		return pnt(x_,y_,z_);
	}

	double magnitude() const{
		return sqrt(x*x + y*y + z*z);
	}

	double magnitude2() const{
		return x*x + y*y + z*z;
	}

};
	
	class triangle {
		public:
			int p1, p2, p3;

		triangle()
			:p1(0), p2(0), p3(0) {}

		triangle(int p1_, int p2_, int p3_)
			:p1(p1_), p2(p2_), p3(p3_) {}

		bool operator==(const triangle &t) const {
		return (((p1 == t.p1) && (p2 == t.p2) && (p3 == t.p3)) || ((p2 == t.p1) && (p3 == t.p2) && (p1 == t.p3)) || ((p3 == t.p1) && (p1 == t.p2) && (p2 == t.p3)));
		}
	};

	class region {
		public:
		int num_halos;
		pnt center;
		vector<pnt> localpointlist; 
		vector<pnt> halopoints;
		vector<triangle> triangles;
		vector<triangle> halotriangles;
	};


	double arcdistance(pnt a, pnt b) {
		double angle = a.dot(b) / (radius*radius);
		angle = acos(angle);
		return fabs(angle * radius);
	}

	double areaOfTriangle(pnt p1, pnt p2, pnt p3) {
		double d1, d2, d3, p;
		d1 = sqrt((p2.x-p1.x)*(p2.x-p1.x) + (p2.y-p1.y)*(p2.y-p1.y) + (p2.z-p1.z)*(p2.z-p1.z)); //p1-p2
		d2 = sqrt((p3.x-p2.x)*(p3.x-p2.x) + (p3.y-p2.y)*(p3.y-p2.y) + (p3.z-p2.z)*(p3.z-p2.z));
		d3 = sqrt((p1.x-p3.x)*(p1.x-p3.x) + (p1.y-p3.y)*(p1.y-p3.y) + (p1.z-p3.z)*(p1.z-p3.z));
		p = (d1+d2+d3)/2;
		return sqrt(p*(p-d1)*(p-d2)*(p-d3));
	}

	pnt calculate_midpoint(pnt a, pnt b) {
		pnt m = (b-a)/2;
		m = a+m;
		m.normalize();
		return m*a.magnitude();	
	}

	pnt calculate_circumcenter(pnt a, pnt b, pnt c) {
		pnt ca = c-a;
		pnt ba = b-a;
		pnt baca = ba.cross(ca);
		pnt bacaba = (baca.cross(ba))*(ca.magnitude2());
		pnt cabaca = (ca.cross(baca))*(ba.magnitude2());
		pnt m = a + (bacaba + cabaca) / (2.0*baca.magnitude2());

		return m;
	}

	double density_function(pnt p) {
		/*
		double lat, lon, r, norm, width, trans_center, min_val;
		convert_xl(p.z, p.y, p.z, &lat, &lon);
		
		r = fabs(lat);
		width = .08;
		trans_center = 0.40;
      min_val = pow((1.0/5.0), 4);
      norm = 1.0/(1.0-min_val);
      return (((tanh((trans_center-r)*(1.0/width))+1.0)/2))/norm + min_val;
		*/
		return 1;
	}





/******************************************************************************************************************************************/
/* isHaloOrNatural :: checks to see if the 3 points provided yield a natural, halo, or discard triangle in the provided region.           */
/******************************************************************************************************************************************/

int isNaturalOrHalo(region r, int a, int b, int c) {
	int count = 0;
	if(a < 0 || b < 0 || c < 0) cout << "\n\n\nERROR: isNaturalOrHalo\n\n" << endl;
	if (a < r.localpointlist.size()) {count++; }
	if (b < r.localpointlist.size()) {count++; }
	if (c < r.localpointlist.size()) {count++; }
	
	if(count == 0) return -1;
	else if(count == 1 || count == 2) return 0;
	else if(count == 3) return 1;
	return -1;
}

void convert_xl(double x, double y, double z, double* lat, double* lon) {
	double clat, eps=1.e-10;
    int i;
    double dl = sqrt(x*x + y*y + z*z);
    *lat = asin(z/dl);
   
    //check for being close to either pole
    if (fabs(x) > eps) {
        if (fabs(y) > eps) {
    		*lon = atan(fabs(y/x));
            if ((x <= 0.0) && (y >= 0.0)) *lon = pii-*lon;
            else if ((x <= 0.0) && (y < 0.0)) *lon = *lon+pii;
            else if ((x >= 0.0) && (y <= 0.0)) *lon = 2*pii-*lon;
         }   
   
         else {// we're either on longitude 0 or 180
            if (x > 0.0) *lon = 0.0;
            else *lon = pii;
		 }
	}
   
     else if (fabs(y) > eps) {
         if (y > 0.0) *lon = pii/2.0;
         else *lon = 3.0*pii/2.0;
     }
   
     else *lon = 0.0; // we are at a pole
   
}

void convert_lx(double* x, double* y, double* z, double lat, double lon, double radius) {
	
	*z = radius * sin(lat);
    *x = radius * cos(lon) * cos(lat);
    *y = radius * sin(lon) * cos(lat);
}

double schmidt_transform(double lat, double beta) {
	double A = beta * (1 + sin(lat)) / (1 - sin(lat)) - 1;
	double B = beta * (1 + sin(lat)) / (1 - sin(lat)) + 1;
	return asin(A / B);
}

/******************************************************************************************************************************************/


int main(int argc, char* argv[]) {
	

	int my_id;
	int nTasks;
	MPI_Datatype MPI_Point;
	MPI_Datatype MPI_Tri;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_Point);
	MPI_Type_contiguous(3, MPI_INT, &MPI_Tri);
	MPI_Type_commit(&MPI_Point);
	MPI_Type_commit(&MPI_Tri);

	int i, j, k;
	clock_t t1, t2, tstart, tend;
	tstart = clock();
	string flags_str = "QBPIOYYiz";
	char * flags;
	flags = new char[flags_str.size()+1];
	strcpy(flags,flags_str.c_str());
	double tPartition=0.0, tTriangulation=0.0, tProjection=0.0, tLloyd=0.0, tSendPoints=0.0, tSendTriangles=0.0;
	int nIterations = stoi(argv[2]);
	vector<pnt*> icoVerticesG1;
	vector<region> regions;

	if (my_id == 0) {

			  /*for (j = 0; j<3; j++) {
				  for (i = 0; i<4; i++) {
					  if (j == 0) icoVerticesG1.push_back(new pnt(0.0,  2.0*(0.0-i%2) + 1.0,  (1.0 + (i < 2.0) * (-2.0))*(1.0+sqrt(5.0))/2.0));
					  else if (j == 1) icoVerticesG1.push_back(new pnt(2.0*(0.0-i%2) + 1.0,  (1.0 + (i < 2.0) * (-2.0))*(1.0+sqrt(5.0))/2.0,  0.0));
					  else if (j == 2) icoVerticesG1.push_back(new pnt((1.0 + (i < 2.0) * (-2.0))*(1.0+sqrt(5.0))/2.0,  0.0,  2.0*(0.0-i%2) + 1.0));
				  }
			  }*/
			  if(argc < 1) {
				  cout << "\n:: Must supply command line of form exec filename" << endl; 
				  return -1;
			  }
			  cout << "Number of iterations to perform: " << nIterations << endl;
			  ifstream pointmaster(argv[1]);
			  double x_, y_, z_;
			  i=0;
			  while (pointmaster >> x_ >> y_ >> z_) {
				  icoVerticesG1.push_back(new pnt(x_, y_, z_));
				  i++;
			  }

			  cout << "There are " << icoVerticesG1.size() << " = " << i << " points in " << argv[1] << endl;
			  /*for (i = 0; i<icoVerticesG1.size(); i++) {
				  icoVerticesG1[i]->id = i;
				  printf("Point %d: (%.2f, %.2f, %.2f)\n", icoVerticesG1[i]->id, icoVerticesG1[i]->x, icoVerticesG1[i]->y, icoVerticesG1[i]->z);
				  pointmaster << icoVerticesG1[i]->x << " " << icoVerticesG1[i]->y << " " << icoVerticesG1[i]->z << " " << icoVerticesG1[i]->id << endl;
			  }*/
					  //string edge(argv[2]);

			  radius = sqrt(x_*x_ + y_*y_ + z_*z_);
			  cout << "Radius is " << radius << endl;


		  /******************************************************************************************************************************************/

			  //setupregions(regions); // setup regions must create the appropriate number of regions and decide their center points
			  regions.push_back(region());
			  regions.push_back(region());
			  regions.push_back(region());
			  regions.push_back(region());
			  regions[0].center = pnt(1.0, 0.0, -1/sqrt(2));
			  regions[1].center = pnt(-1.0, 0.0, -1/sqrt(2));
			  regions[2].center = pnt(0.0, 1.0, 1/sqrt(2));
			  regions[3].center = pnt(0.0, -1.0, 1/sqrt(2));
			  ofstream centers("rcenters.txt");	
			  for (i=0; i<regions.size(); i++){
				  regions[i].center.normalize();	
				  printf("Region %d Center: (%.2f, %.2f, %.2f)\n", i, regions[i].center.x, regions[i].center.y, regions[i].center.z);
				  centers << regions[i].center.x << "   " << regions[i].center.y << "   " << regions[i].center.z << endl;
			  }

			}


		  /*
		  pnt test1(0.0, 0.0, -1.0);
		  pnt test2(0.0, 1.0, 0.0);
		  pnt test3(1.0, 0.0, 0.0);

		  double lat, lon;

		  convert_xl(test1.x, test1.y, test1.z, &lat, &lon);
		  printf("(%f, %f, %f) -> %6f, %6f", test1.x, test1.y, test1.z, lat, lon);
		  convert_lx(&test1.x, &test1.y, &test1.z, lat, lon, 1.0);
		  lat = schmidt_transform(lat, beta);
		  printf("->(%f, %f, %f) with transformed l/l %6f, %6f\n", test1.x, test1.y, test1.z, lat, lon);

		  convert_xl(test2.x, test2.y, test2.z, &lat, &lon);
		  printf("(%f, %f, %f) -> %6f, %6f", test2.x, test2.y, test2.z, lat, lon);
		  convert_lx(&test2.x, &test2.y, &test2.z, lat, lon, 1.0);
		  lat = schmidt_transform(lat, beta);
		  printf("->(%f, %f, %f) with transformed l/l %6f, %6f\n", test2.x, test2.y, test2.z, lat, lon);

		  convert_xl(test3.x, test3.y, test3.z, &lat, &lon);
		  printf("(%f, %f, %f) -> %6f, %6f", test3.x, test3.y, test3.z, lat, lon);
		  convert_lx(&test3.x, &test3.y, &test3.z, lat, lon, 1.0);
		  lat = schmidt_transform(lat, beta);
		  printf("->(%f, %f, %f) with transformed l/l %6f, %6f\n", test3.x, test3.y, test3.z, lat, lon);
		  */

			  int iIteration = 0;
			  int ipercent = 0;
			  tstart = clock();
			  if (my_id == 0) cout << "0%%" << flush;
			  while(iIteration < nIterations){
				  iIteration++;
				if (my_id == 0) {
				  if (int(100.0 * float(iIteration) / float(nIterations)) > ipercent) {
					  ipercent = int(100.0 * float(iIteration) / float(nIterations));
					  printf("\b\b\b%02d%%", ipercent);
					  cout << flush;
				  }
				  double d;
				  t2 = clock();
				  //cout << "Partitioning Points... " << flush;
				  for (i=0; i<regions.size(); i++) {
					  if (!regions[i].localpointlist.empty()) regions[i].localpointlist.clear();
					  if (!regions[i].triangles.empty()) regions[i].triangles.clear();
				  }
				  for (i = 0; i < icoVerticesG1.size(); i++) {
					  int closestRegion = 0;
					  double min = arcdistance(*icoVerticesG1[i], regions[0].center);
					  for (j = 1; j < regions.size(); j++){
						  d = arcdistance(*icoVerticesG1[i], regions[j].center);

						  if (d < min) {
							  closestRegion = j;
							  min = d;
						  }
					  }
					  regions[closestRegion].localpointlist.push_back(*icoVerticesG1[i]);
				  }
				  t1 = clock();
				  tPartition += (float)(t1-t2)/CLOCKS_PER_SEC;
				  //cout << (float)(t1-t2)/CLOCKS_PER_SEC << "s" << endl;

			  /*	for (i=0; i<regions.size(); i++) {
					  cout << "Region " << i << "local points:" << endl;
					  for (j=0; j<regions[i].localpointlist.size(); j++) {
						  printf("\t%d: (%.2f, %.2f, %.2f)\n", regions[i].localpointlist[j].id, regions[i].localpointlist[j].x, regions[i].localpointlist[j].y, regions[i].localpointlist[j].z);
					  }
					  cout << "\nRegion " << i << "halo points:" << endl;
					  for (j=0; j<regions[i].halopoints.size(); j++) {
						  printf("\t%d: (%.2f, %.2f, %.2f)\n", regions[i].halopoints[j].id, regions[i].halopoints[j].x, regions[i].halopoints[j].y, regions[i].halopoints[j].z);
					  }
				  }*/
				  //cout << "Checking for errors in point lists... " << endl;
				  int errcount = 0;
				  /*for (i=0; i<regions[0].localpointlist.size(); i++) {
					  for (j = 0; j< regions[1].localpointlist.size(); j++) {
						  if (regions[0].localpointlist[i] == regions[1].localpointlist[j]) {cout << "\n\n\n ERROR: a point is in the local point list of both regions\n\n" << endl; errcount++;}
					  }
				  }*/
				  //cout << "Errcount: " << errcount << endl;

				  //for (i=0; i<regions.size(); i++) {
				  //	cout << "Region " << i << " has " << regions[i].localpointlist.size() << " natural points and " << regions[i].halopoints.size() << " halo points." << endl;
				  //}

			  /******************************************************************************************************************************************/

				  vector<double> sizes;
				  double area;
				  int found = 0;
				  int numtri = 0;
			}
			
		
		MPI_Barrier(MPI_COMM_WORLD);
	printf("Task %d checking in.\n", my_id);

		region myRegion;
		int nPoints;
	
		if (my_id == 0) {
			t2 = clock();
			for(i=1; i<nTasks; i++) {
				double iCenter[3];
				int iNpts = regions[i].localpointlist.size();
				iCenter[0] = regions[i].center.x;
				iCenter[1] = regions[i].center.y;
				iCenter[2] = regions[i].center.z;

				MPI_Send(&iCenter, 1, MPI_Point, i, 1, MPI_COMM_WORLD);
				MPI_Send(&iNpts, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
				
				MPI_Send(regions[i].localpointlist.data(), regions[i].localpointlist.size(), MPI_Point, i, 3, MPI_COMM_WORLD);
			}
			myRegion.center = regions[0].center;
			myRegion.localpointlist.clear();
			myRegion.localpointlist = regions[0].localpointlist;
			nPoints = myRegion.localpointlist.size();
			t1 = clock();
			tSendPoints+= (float)(t1-t2)/CLOCKS_PER_SEC;
		}

		else {
			double myCenter[3];
			myRegion.localpointlist.clear();
			MPI_Recv(&myCenter, 1, MPI_Point, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("Region %d has center (%f, %f, %f)\n", my_id, myCenter[0], myCenter[1], myCenter[2]);
			MPI_Recv(&nPoints, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("Task %d sees there are %d points.\n", my_id, nPoints);
			pnt *myPoints = new pnt[nPoints];
			MPI_Recv(myPoints, nPoints, MPI_Point, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			myRegion.localpointlist.insert(myRegion.localpointlist.begin(), myPoints, myPoints+nPoints);
			delete myPoints;
		}


			//cout << "\nRegion " << i << " Loop" << endl;
			t2 = clock();
			struct triangulateio in, out, vorout;
			pnt center = myRegion.center;
			pnt axis(center.x, center.y, center.z);
			double min_dir = min(fabs(axis.x),min(fabs(axis.y),fabs(axis.z)));

			if(min_dir == fabs(axis.x)){
				axis.x = 1.0;
				axis.y = 0.0;
				axis.z = 0.0;
			} else if (min_dir == fabs(axis.y)){
				axis.x = 0.0;
				axis.y = 1.0;
				axis.z = 0.0;
			} else if (min_dir == fabs(axis.z)){
				axis.x = 0.0;
				axis.y = 0.0;
				axis.z = 1.0;
			}

			pnt x_hat = center.cross(axis);
			x_hat.normalize();
			pnt y_hat = center.cross(x_hat);
			y_hat.normalize();
			double s;
			pnt Q;

			in.numberofpoints = myRegion.localpointlist.size() + myRegion.halopoints.size();
			in.numberofpointattributes = 0;
			in.pointlist = (double *) malloc(in.numberofpoints * 2 * sizeof(double));
			in.numberofsegments = 0;
			in.numberofholes = 0;
			in.numberofregions = 0;
			in.regionlist = (double *) NULL;
			in.pointmarkerlist = (int *) NULL;

			for (j = 0; j < myRegion.localpointlist.size(); j++) {
				
				s = 2.0/center.dot(center + (myRegion.localpointlist[j]));
				Q = myRegion.localpointlist[j] * s + center * (s - 1.0);
				in.pointlist[j*2] = x_hat.dot(Q);
				in.pointlist[j*2+1] = y_hat.dot(Q);			
			}

			out.pointlist = (double *)NULL;
			out.trianglelist = (int *)NULL;
			//t1 = clock();
			//tProjection += (float)(t1-t2)/CLOCKS_PER_SEC;
			//cout << (float)(t1-t2)/CLOCKS_PER_SEC << "s for projection" << endl;
			//cout << "Triangulating ... " << flush;
			triangulate(flags,&in,&out,&vorout);
			//t2 = clock();
			//tTriangulation += (float)(t2-t1)/CLOCKS_PER_SEC;
			//cout << (float)(t2-t1)/CLOCKS_PER_SEC << "s" << endl;

			j = 0;
			int l = 0;
			k = 1;
			
			//cout << "Processing Triangles and performing Lloyd iteration... 0%%" << flush;
			//t1 = clock();
			int percent = 0;
			for (k = 0; k < out.numberoftriangles; k++) {
				//if (100 * k / out.numberoftriangles > percent) {
				//	percent = 100 * k / out.numberoftriangles; 
				//	printf("\b\b\b%02d%%", percent);
				//	cout << flush;
				//}

				//regions[i].triangles.push_back(triangle(l, regions[i].localpointlist[out.trianglelist[3*k]].id, regions[i].localpointlist[out.trianglelist[3*k+1]].id, regions[i].localpointlist[out.trianglelist[3*k+2]].id));
				myRegion.triangles.push_back(triangle(out.trianglelist[3*k], out.trianglelist[3*k+1], out.trianglelist[3*k+2]));
				//numtri++;
			}


			MPI_Barrier(MPI_COMM_WORLD);
			if (my_id == 0) {
				t1 = clock();
				tTriangulation += (float)(t1-t2)/CLOCKS_PER_SEC;
			}

			/* 
			 * Reduce all triangles to task 0
			 */

			if (my_id == 0) {
				t2 = clock();
				int nTris[nTasks];
				regions[0].triangles = myRegion.triangles;
				for (i=1; i<nTasks; i++) {
					MPI_Recv(nTris+i, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					triangle *iTriangles = new triangle[nTris[i]];
					MPI_Recv(iTriangles, nTris[i], MPI_Tri, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					regions[i].triangles.clear();
					regions[i].triangles.insert(regions[i].triangles.begin(), iTriangles, iTriangles+nTris[i]);
					delete iTriangles;
				}
				t1 = clock();
				tSendTriangles+= (float)(t1-t2)/CLOCKS_PER_SEC;
			}

			else {
				int ntris = myRegion.triangles.size();
				MPI_Send(&ntris, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
				MPI_Send(myRegion.triangles.data(), ntris, MPI_Tri, 0, 2, MPI_COMM_WORLD);
			}



		if (my_id == 0) {
			t2 = clock();
			for (i=0; i<nTasks; i++) {
				pnt* newCenter =  new pnt[regions[i].localpointlist.size()];
				double* areaCell = new double[regions[i].localpointlist.size()];
				long npts = regions[i].localpointlist.size();
				// iterate through all triangles computing weights
				for (k = 0; k<regions[i].triangles.size(); k++) {
					triangle tri_k = regions[i].triangles[k];
					if(tri_k.p1 >= npts || tri_k.p2 >= npts || tri_k.p3 >= npts) {
						cout << tri_k.p1 << ", "<< tri_k.p2 << ", "<< tri_k.p3 << " :: " << npts << endl;
					}
	
					pnt p1 = regions[i].localpointlist[tri_k.p1];
					pnt p2 = regions[i].localpointlist[tri_k.p2];
					pnt p3 = regions[i].localpointlist[tri_k.p3];
					
					pnt r = calculate_circumcenter(p1, p2, p3);
	
					pnt w12 = calculate_midpoint(p1, p2);
					pnt w23 = calculate_midpoint(p2, p3);
					pnt w31 = calculate_midpoint(p3, p1);
	
					pnt q1 = calculate_midpoint(p1, r);
					pnt q2 = calculate_midpoint(p2, r);
					pnt q3 = calculate_midpoint(p3, r);
					double area = areaOfTriangle(p1, w31, r) + areaOfTriangle(p1, w12, r);
					newCenter[tri_k.p1] += q1 * area * density_function(q1);	
					areaCell[tri_k.p1] += area;
	
					area = areaOfTriangle(p2, w12, r) + areaOfTriangle(p2, w23, r);
					newCenter[tri_k.p2] += q2 * area * density_function(q2);	
					areaCell[tri_k.p2] += area;
	
					area = areaOfTriangle(p3, w23, r) + areaOfTriangle(p3, w31, r);
					newCenter[tri_k.p3] += q3 * area * density_function(q3);	
					areaCell[tri_k.p3] += area;
				}

				// iterate through all corners and relax according to weights computed

				for (k=0; k<regions[i].localpointlist.size(); k++) {
					newCenter[k] = regions[i].localpointlist[k] + newCenter[k] / areaCell[k];
				}
			}			
			t1 = clock();
			tLloyd += (float)(t1-t2)/CLOCKS_PER_SEC;
		}
	}

	//int numtri2 = numtri;
	//numtri = 0;

/*
	cout << "\nWriting triangles to file ... 0%%" << flush;
	t1 = clock();
	ofstream triangleout("output_triangulation.txt");
	int percent = 0;
	for (i = 0; i < regions.size(); i++) {
		for (j=0; j<regions[i].triangles.size(); j++){
			triangleout << regions[i].triangles[j].p1 << " " << regions[i].triangles[j].p2 << " " << regions[i].triangles[j].p3 << endl;
			numtri++;
			area = areaOfTriangle(*icoVerticesG1[regions[i].triangles[j].p1], *icoVerticesG1[regions[i].triangles[j].p2], *icoVerticesG1[regions[i].triangles[j].p3]);
			found = 0;
			for (k=0; k<sizes.size(); k++) {
				if (area == sizes[k]) found = 1;
			}
			if (!found) sizes.push_back(area);
		if (100 * numtri / numtri2 > percent) {
            percent = 100 * numtri / numtri2;
            printf("\b\b\b%02d%%", percent);
            cout << flush;
        }
		}
		for (j=0; j<regions[i].halotriangles.size(); j++){
			triangleout << regions[i].halotriangles[j].p1 << " " << regions[i].halotriangles[j].p2 << " " << regions[i].halotriangles[j].p3 << endl;
			numtri++;
			area = areaOfTriangle(*icoVerticesG1[regions[i].halotriangles[j].p1], *icoVerticesG1[regions[i].halotriangles[j].p2], *icoVerticesG1[regions[i].halotriangles[j].p3]);
			found = 0;
			for (k=0; k<sizes.size(); k++) {
				if (area == sizes[k]) found = 1;
			}
			if (!found) sizes.push_back(area);
		if (100 * numtri / numtri2 > percent) {
            percent = 100 * numtri / numtri2;
            printf("\b\b\b%02d%%", percent);
            cout << flush;
        }
		}
	}

	t2 = clock();
	cout << "\b\b\b" << (float)(t2-t1)/CLOCKS_PER_SEC << "s" << endl;

	cout << "There are " << numtri << "=" << numtri2 << "=" << regions[0].triangles.size() + regions[0].halotriangles.size() + regions[1].triangles.size() + regions[1].halotriangles.size() << "triangles" << endl;
	cout << "There are " << sizes.size() << " different areas of these triangles:" << endl;
	double biggest = 0.0;
	double smallest = 1000000000000000000.0;

	for (i=0; i<sizes.size(); i++) {
		if (sizes[i] < smallest) smallest = sizes[i];
		else if (sizes[i] > biggest) biggest = sizes[i];
	}

	cout << "Largest Area: " << biggest << "\nSmallest Area: " << smallest << endl;

*/

	if (my_id == 0) {
	tend = clock();
	cout << "Wall clock time in serial: " << (float)(tend-tstart)/CLOCKS_PER_SEC << endl;
	cout << "Total tPartition: " << tPartition << endl;
	cout << "Total tTriangulation (nTasks == " << nTasks << "): " << tTriangulation << endl;
	cout << "Total tLloyd: " << tLloyd << "\n" << endl;
	cout << "nIterations: " << nIterations << endl;
	cout << "Time in ideally parallelized Lloyd over " << nTasks << " processors: " << tPartition + tTriangulation + tLloyd / float(nTasks) << endl;
 	}

	MPI_Finalize();

	return 0;

}
