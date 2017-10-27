
#include "geometry_powder_target_3d.h"
#include <iostream>
#include <cmath>


////////////////////////////////////////////////////////////////////////////////////////
// Start of PowderTarget3D
////////////////////////////////////////////////////////////////////////////////////////

PowderTarget3D::PowderTarget3D():radius(0.75), height(0.75), length(4.0) {}

bool PowderTarget3D::operator()(double x, double y, double z) const {
  double ht;
//  double ht, L1 = 0.6, L2 = 0.8;

  //    ht = height + 0.06*sin(2*M_PI*x/L1) * sin(2*M_PI*z/L2);
  ht = height;
  if (y <= 0)
    return (x*x + y*y < radius*radius) && (fabs(z) < length/2.0 );
  else
    return (fabs(x) < radius && y < ht) && (fabs(z) < length/2.0); 	
}

void PowderTarget3D::getBoundingBox(double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax) {
	xmin = -radius;
	xmax = radius;
	ymin = -radius;
	ymax = height + 0.1;
	zmin = -length/2.0;
	zmax = length/2.0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of PowderTarget3D
////////////////////////////////////////////////////////////////////////////////////////
