#include <iostream>
#include "vtkfile_processor.h"
#include <cmath>

const double PI = atan(1) * 4.;

using namespace std;

int main(int argc, const char* argv[]) {
	
	VTKFileProcessor p1;
	p1.read("out_fluid0000025.vtk");
	p1.changeVelocity(0,300,0);
	p1.rotate(66./180.*PI);
	p1.write("out_fluid_rotate0000025.vtk");
	
	VTKFileProcessor p2;
	p2.read("out_fluid0000025.vtk");
	p2.changeVelocity(0,300,0);
	p2.changeTag(2);
	p2.rotate(114./180.*PI);	
	p2.translate(0,24,0);
	


	p2.combine("out_fluid_rotate0000025.vtk");
	p2.write("out_fluid_combined0000025.vtk");

	return 0;
}
