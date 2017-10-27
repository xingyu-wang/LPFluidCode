#include <iostream>
#include "vtkfile_processor.h"
#include <cmath>

const double PI = atan(1) * 4.;

using namespace std;

int main(int argc, const char* argv[]) {

// 2d shock tube problem
/*
	double yCen = 0;
	double xLeft = -4.95, xRight = 4.95, spacing = 0.05;

	for(size_t step=0; step<=200; step++) {
		string ifname = "data/out_fluid" + VTKFileProcessor::rightFlush(step,7) + ".vtk";
		string ofname = "output/out_fluid_1d" + VTKFileProcessor::rightFlush(step,7);
	
		VTKFileProcessor p1;
		p1.read(ifname);	
		p1.reduce2dY(yCen,xLeft,xRight,spacing,ofname);
	}
*/	
	
// 3d gas ball expansion into vacuum

	double yCen = 0, zCen = 0;
	double xLeft = -0.4, xRight = 0.4, spacing = 0.01;

	for(size_t step=0; step<=8; step++) {
		string ifname = "data/out_fluid" + VTKFileProcessor::rightFlush(step,7) + ".vtk";
		string ofname = "output/out_fluid_1d" + VTKFileProcessor::rightFlush(step,7);
	
		VTKFileProcessor p1;
		p1.read(ifname);	
		p1.reduce3dYZ(yCen,zCen,xLeft,xRight,spacing,ofname);
	}

	return 0;
}
