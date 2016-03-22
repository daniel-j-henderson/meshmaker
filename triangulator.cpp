#define REAL double
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


using namespace std;

double radius = 0.0;
REAL pii = 2.0 * asin(1.0);

class pnt {
public:
	REAL x, y, z;
	int id;

	pnt(double x_, double y_, double z_)
		: x(x_), y(y_), z(z_) { }

	pnt(int id_, double x_, double y_, double z_)
		: id(id_), x(x_), y(y_), z(z_) { }

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
			int id;
			int p1, p2, p3;

		triangle() 
			:id(-1) {}

		triangle(int p1_, int p2_, int p3_)
			:id(-1), p1(p1_), p2(p2_), p3(p3_) {}

		triangle(int id_, int p1_, int p2_, int p3_)
			:id(id_), p1(p1_), p2(p2_), p3(p3_) {}

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


	REAL arcdistance(pnt a, pnt b) {
		REAL angle = a.dot(b) / (radius*radius);
		angle = acos(angle);
		return fabs(angle * radius);
	}

	REAL areaOfTriangle(pnt p1, pnt p2, pnt p3) {
		REAL d1, d2, d3, p;
		d1 = sqrt((p2.x-p1.x)*(p2.x-p1.x) + (p2.y-p1.y)*(p2.y-p1.y) + (p2.z-p1.z)*(p2.z-p1.z)); //p1-p2
		d2 = sqrt((p3.x-p2.x)*(p3.x-p2.x) + (p3.y-p2.y)*(p3.y-p2.y) + (p3.z-p2.z)*(p3.z-p2.z));
		d3 = sqrt((p1.x-p3.x)*(p1.x-p3.x) + (p1.y-p3.y)*(p1.y-p3.y) + (p1.z-p3.z)*(p1.z-p3.z));
		p = (d1+d2+d3)/2;
		return sqrt(p*(p-d1)*(p-d2)*(p-d3));
	}





/******************************************************************************************************************************************/
/* isHaloOrNatural :: checks to see if the 3 points provided yield a natural, halo, or discard triangle in the provided region.
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

void convert_xl(REAL x, REAL y, REAL z, REAL* lat, REAL* lon) {
	REAL clat, eps=1.e-10;
    int i;
    REAL dl = sqrt(x*x + y*y + z*z);
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

void convert_lx(REAL* x, REAL* y, REAL* z, REAL lat, REAL lon, REAL radius) {
	
	*z = radius * sin(lat);
    *x = radius * cos(lon) * cos(lat);
    *y = radius * sin(lon) * cos(lat);
}

double schmidt_transform(REAL lat, REAL beta) {
	REAL A = beta * (1 + sin(lat)) / (1 - sin(lat)) - 1;
	REAL B = beta * (1 + sin(lat)) / (1 - sin(lat)) + 1;
	return asin(A / B);
}

/******************************************************************************************************************************************/


int main(int argc, char* argv[]) {
	int i, j, k;
	clock_t t1, t2;
	string flags_str = "QBPIOYYiz";
	char * flags;
	flags = new char[flags_str.size()+1];
	strcpy(flags,flags_str.c_str());

	vector<pnt*> icoVerticesG1;
	/*for (j = 0; j<3; j++) {
		for (i = 0; i<4; i++) {
			if (j == 0) icoVerticesG1.push_back(new pnt(0.0,  2.0*(0.0-i%2) + 1.0,  (1.0 + (i < 2.0) * (-2.0))*(1.0+sqrt(5.0))/2.0));
			else if (j == 1) icoVerticesG1.push_back(new pnt(2.0*(0.0-i%2) + 1.0,  (1.0 + (i < 2.0) * (-2.0))*(1.0+sqrt(5.0))/2.0,  0.0));
			else if (j == 2) icoVerticesG1.push_back(new pnt((1.0 + (i < 2.0) * (-2.0))*(1.0+sqrt(5.0))/2.0,  0.0,  2.0*(0.0-i%2) + 1.0));
		}
	}*/
	if(argc < 3) {
		cout << "\n:: Must supply command line of form exec filename edgelength betavalue" << endl; 
		return -1;
	}
	ifstream pointmaster(argv[1]);
	REAL x_, y_, z_;
	i=0;
	while (pointmaster >> x_ >> y_ >> z_) {
		icoVerticesG1.push_back(new pnt(i, x_, y_, z_));
		i++;
	}

	cout << "There are " << icoVerticesG1.size() << " = " << i << " points in " << argv[1] << endl;
	/*for (i = 0; i<icoVerticesG1.size(); i++) {
		icoVerticesG1[i]->id = i;
		printf("Point %d: (%.2f, %.2f, %.2f)\n", icoVerticesG1[i]->id, icoVerticesG1[i]->x, icoVerticesG1[i]->y, icoVerticesG1[i]->z);
		pointmaster << icoVerticesG1[i]->x << " " << icoVerticesG1[i]->y << " " << icoVerticesG1[i]->z << " " << icoVerticesG1[i]->id << endl;
	}*/
			//string edge(argv[2]);

	double edgelength = std::stod(argv[2]);
	double beta;
	if(argc == 4) beta = stod(argv[3]);
	else beta = 0.0; // No Shmidt Transform
	
	radius = sqrt(x_*x_ + y_*y_ + z_*z_);



/******************************************************************************************************************************************/

	vector<region> regions;
	//setupregions(regions); // setup regions must create the appropriate number of regions and decide their center points
	regions.push_back(region());
	regions.push_back(region());
	regions[0].center.z = radius;
	regions[1].center.z = -radius;
	printf("Region 0 Center: (%.2f, %.2f, %.2f)\n", regions[0].center.x, regions[0].center.y, regions[0].center.z);
	printf("Region 1 Center: (%.2f, %.2f, %.2f)\n", regions[1].center.x, regions[1].center.y, regions[1].center.z);

	ofstream latstream("original_lats.txt");	

	REAL d = -1;
	if (beta != 0.0) {
		for (i=0; i<icoVerticesG1.size(); i++) {
			REAL lat, lon;
			convert_xl(icoVerticesG1[i]->x, icoVerticesG1[i]->y, icoVerticesG1[i]->z, &lat, &lon);
			latstream << lat << endl;
			if(lat != pii / 2) lat = schmidt_transform(lat, beta);
			convert_lx(&icoVerticesG1[i]->x, &icoVerticesG1[i]->y, &icoVerticesG1[i]->z, lat, lon, radius);
		}
	cout << "Performed Schmidt Transform" << endl;
	}

	ofstream pointsout("output_points.txt");
	for (i=0; i<icoVerticesG1.size(); i++) {
		pointsout << icoVerticesG1[i]->x << "    " << icoVerticesG1[i]->y << "    " << icoVerticesG1[i]->z << endl;
	}
	



/*
pnt test1(0.0, 0.0, -1.0);
pnt test2(0.0, 1.0, 0.0);
pnt test3(1.0, 0.0, 0.0);

REAL lat, lon;

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
	cout << "Generating Pointlists... \r" << endl;
	for (i = 0; i < icoVerticesG1.size(); i++) {
		int closestRegion = 0;
		REAL min = arcdistance(*icoVerticesG1[i], regions[0].center);
		for (j = 1; j < regions.size(); j++){
			d = arcdistance(*icoVerticesG1[i], regions[j].center);

			if (d < min) {
				closestRegion = j;
				min = d;
			}
		}
		regions[closestRegion].localpointlist.push_back(*icoVerticesG1[i]);
		for (j = 0; j < regions.size(); j++) {
			if (j != closestRegion) {
				d = arcdistance(*icoVerticesG1[i], regions[j].center);
				if (fabs(min - d) < 1.5*edgelength) {
					regions[j].halopoints.push_back(*icoVerticesG1[i]);
				}
			}
		}
		
	}

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
	cout << "Checking for errors in point lists... \r" << endl;
	int errcount = 0;
	for (i=0; i<regions[0].localpointlist.size(); i++) {
		for (j = 0; j< regions[1].localpointlist.size(); j++) {
			if (regions[0].localpointlist[i] == regions[1].localpointlist[j]) {cout << "\n\n\n ERROR: a point is in the local point list of both regions\n\n" << endl; errcount++;}
		}
	}

	for (i=0; i<regions.size(); i++) {
		for (j=0; j<regions[i].localpointlist.size(); j++) {
			for (k=0; k < regions[i].halopoints.size(); k++) {
				if (regions[i].localpointlist[j] == regions[i].halopoints[k]) {errcount++; cout << "\n\n\nERROR: a point is in both the halo and local lists for a region\n\n" << endl;}
			}
		}
	}
	cout << "Errcount: " << errcount << endl;

	for (i=0; i<regions.size(); i++) {
		cout << "Region " << i << " has " << regions[i].localpointlist.size() << " natural points and " << regions[i].halopoints.size() << " halo points." << endl;
	}

/******************************************************************************************************************************************/

	vector<REAL> sizes;
	REAL area;
	int found = 0;
	int numtri = 0;

	for (i=0; i<regions.size(); i++) {
		cout << "\nRegion " << i << " Loop" << endl;

		struct triangulateio in, out, vorout;
		pnt center = regions[i].center;
		pnt axis(center.x, center.y, center.z);
		REAL min_dir = min(fabs(axis.x),min(fabs(axis.y),fabs(axis.z)));

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

		in.numberofpoints = regions[i].localpointlist.size() + regions[i].halopoints.size();
		in.numberofpointattributes = 0;
		in.pointlist = (double *) malloc(in.numberofpoints * 2 * sizeof(double));
		in.numberofsegments = 0;
		in.numberofholes = 0;
		in.numberofregions = 0;
		in.regionlist = (double *) NULL;
		in.pointmarkerlist = (int *) NULL;

		for (j = 0; j < regions[i].localpointlist.size(); j++) {
			
			s = 2.0/center.dot(center + (regions[i].localpointlist[j]));
			Q = regions[i].localpointlist[j] * s + center * (s - 1.0);
			in.pointlist[j*2] = x_hat.dot(Q);
			in.pointlist[j*2+1] = y_hat.dot(Q);			
		}

		for (k = 0; k<regions[i].halopoints.size(); k++) {
			s = 2.0/center.dot(center + (regions[i].halopoints[k]));
			Q = regions[i].halopoints[k] * s + center * (s - 1.0);
			in.pointlist[j*2] = x_hat.dot(Q);
			in.pointlist[j*2+1] = y_hat.dot(Q);
			j++;
		}

		out.pointlist = (double *)NULL;
		out.trianglelist = (int *)NULL;
		cout << "Triangulating ... " << flush;
		t1 = clock();
		triangulate(flags,&in,&out,&vorout);
		t2 = clock();
		cout << (float)(t2-t1)/CLOCKS_PER_SEC << "s" << endl;

		j = 0;
		int l = 0;
		k = 1;
		
		cout << "Processing Triangles ... 0%%" << flush;
		t1 = clock();
		int percent = 0;
		for (k = 0; k < out.numberoftriangles; k++) {
		if (100 * k / out.numberoftriangles > percent) {
			percent = 100 * k / out.numberoftriangles; 
			printf("\b\b\b%02d%%", percent);
			cout << flush;
		}
			if ((out.trianglelist[3*k] != out.trianglelist[3*k+1]) && (out.trianglelist[3*k] != out.trianglelist[3*k+2]) && (out.trianglelist[3*k+2] != out.trianglelist[3*k+1])) {
				int test = isNaturalOrHalo(regions[i], out.trianglelist[3*k], out.trianglelist[3*k+1], out.trianglelist[3*k+2]);
				switch (test)
				{
					default :
						break;
					case 0:
						int x, y, z;
						if(out.trianglelist[3*k] < regions[i].localpointlist.size()) x = regions[i].localpointlist[out.trianglelist[3*k]].id;
						else x = regions[i].halopoints[out.trianglelist[3*k] - regions[i].localpointlist.size()].id;

						if(out.trianglelist[3*k+1] < regions[i].localpointlist.size()) y = regions[i].localpointlist[out.trianglelist[3*k+1]].id;
						else y = regions[i].halopoints[out.trianglelist[3*k+1] - regions[i].localpointlist.size()].id;

						if(out.trianglelist[3*k+2] < regions[i].localpointlist.size()) z = regions[i].localpointlist[out.trianglelist[3*k+2]].id;
						else z = regions[i].halopoints[out.trianglelist[3*k+2] - regions[i].localpointlist.size()].id;

						regions[i].halotriangles.push_back(triangle(j, x, y, z));
						j++;
						if(test == 0) numtri++;
						break;
					case 1:
						regions[i].triangles.push_back(triangle(l, regions[i].localpointlist[out.trianglelist[3*k]].id, regions[i].localpointlist[out.trianglelist[3*k+1]].id, regions[i].localpointlist[out.trianglelist[3*k+2]].id));
						l++;
						if(test == 1) numtri++;
						break;

				}
			}
		}


		
		t2 = clock();
		cout << "\b\b\b" << (float)(t2-t1)/CLOCKS_PER_SEC << "s" << endl;
		cout << "And erasing duplicates ..." << flush;
		t1 = clock();
		int numer = 0;
		
		for (k=0; k<i; k++) {
			int cond = 0;
			vector<triangle>::iterator tri_itr = regions[i].halotriangles.begin();
			while (tri_itr !=regions[i].halotriangles.end()) { //for (tri_itr = regions[i].halotriangles.begin(); tri_itr !=regions[i].halotriangles.end(); tri_itr++) {
				for (l = 0; l < regions[k].halotriangles.size(); l++) {
					cond = 0;
					if ((*tri_itr) == regions[k].halotriangles[l]) {
						numer++;
						tri_itr = regions[i].halotriangles.erase(tri_itr);
						l = regions[k].halotriangles.size();
						cond = 1;
					}
				}
				if (!cond) tri_itr++;
			}
		}

		t2 = clock();
		cout << (float)(t2-t1)/CLOCKS_PER_SEC << "s" << endl;
		cout << "Erased " << numer << " duplicate triangles" << endl;
		numtri -= numer;
	}


	int numtri2 = numtri;
	numtri = 0;

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
	REAL biggest = 0.0;
	REAL smallest = 1000000000000000000.0;

	for (i=0; i<sizes.size(); i++) {
		if (sizes[i] < smallest) smallest = sizes[i];
		else if (sizes[i] > biggest) biggest = sizes[i];
	}

	cout << "Largest Area: " << biggest << "\nSmallest Area: " << smallest << endl;

	
	return 0;

}
