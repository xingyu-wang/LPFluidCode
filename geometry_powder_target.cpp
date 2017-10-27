#include "geometry_powder_target.h"
#include <iostream>
#include <cmath>


////////////////////////////////////////////////////////////////////////////////////////
// Start of PowderTarget2D
////////////////////////////////////////////////////////////////////////////////////////

PowderTarget2D::PowderTarget2D():radius(0.6), height(0.6) {}

bool PowderTarget2D::operator()(double x, double y, double z) const {
  if (y <= 0)
    return (x*x + y*y < radius*radius);
  else
    return (fabs(x) < radius && y < height); 	
}

void PowderTarget2D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = -radius;
	xmax = radius;
	ymin = -radius;
	ymax = height;
	zmin = 0;
	zmax = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of PowderTarget2D
////////////////////////////////////////////////////////////////////////////////////////
