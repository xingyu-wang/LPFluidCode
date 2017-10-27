#include "geometry_ballexp.h"
#include <iostream>
#include <cmath>


////////////////////////////////////////////////////////////////////////////////////////
// Start of Ballexp3D
////////////////////////////////////////////////////////////////////////////////////////

Ballexp3D::Ballexp3D():radius(0.2), xCen(0), yCen(0), zCen(0) {}

bool Ballexp3D::operator()(double x, double y, double z) const {	
	double xd = x-xCen, yd = y-yCen, zd = z-zCen;

	double dist = sqrt(xd*xd+yd*yd+zd*zd);	

	return dist <= radius; 	
}

void Ballexp3D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = xCen-radius;
	xmax = xCen+radius;
	ymin = yCen-radius;
	ymax = yCen+radius;
	zmin = zCen-radius;
	zmax = zCen+radius;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Ballexp3D
////////////////////////////////////////////////////////////////////////////////////////
