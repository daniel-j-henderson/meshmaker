#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <netcdf>

using namespace netCDF;
using namespace netCDF::exceptions;
using namespace std;

double pii = 2.0 * asin(1.0);

int main(int argc, char* argv[]) {

	if(argc != 4) {
        cout << "\n:: Must supply command line of form exec inputfilename outputfilename betavalue" << endl;
        return -1;
    }
	// argv[1] = filename. argv[2] = beta
	double beta = stoi(argv[3]);
		cout << "A" << endl;
	vector<double> lats;
	ifstream latstream(argv[1]);
	double lat;
    while (latstream >> lat) {
        lats.push_back(lat);
    }
	NcFile ncfile(argv[2], NcFile::write);
	NcDim nCells = ncfile.getDim("nCells");
	    cout << "B" << endl;
	NcVar meshDensity = ncfile.getVar("meshDensity");
        cout << "C" << endl;
	//NcVar latCell = ncfile.getVar("latCell");
        cout << "C" << endl;
	size_t num_cells = nCells.getSize();
	//double* lats = new double[num_cells];
	//latCell.getVar(lats);
        cout << "C" << endl;
	if(num_cells != lats.size()) {
		cout << "nCells is different than the number of input lats" << endl;
		return -1;
	}
	double mesh_density[lats.size()];
	double dx = 2 * sqrt(beta) / ((beta+1) + (beta-1) * sin(pii / 2));
        cout << "D" << endl;
	for (int i=0; i<num_cells; i++) {
		double ds = 2 * sqrt(beta) / ((beta+1) + (beta-1) * sin(lats[i]));
		mesh_density[i] = pow(dx / ds, 4);
	}
	vector<size_t> start, count;
	start.push_back(0);
	count.push_back(num_cells);
        cout << "E" << endl;
	meshDensity.putVar(start, count, mesh_density);
	return 0;
}
