#include "state_powder_target_3d.h"
#include <iostream>
#include <cmath>


////////////////////////////////////////////////////////////////////////////////////////
// Start of PowderTarget3DState
////////////////////////////////////////////////////////////////////////////////////////

PowderTarget3DState::PowderTarget3DState(): 

  //Thimble parameters
  rho_0(9.625), p0(0.001), p_deposition(1.0), sigma_x(0.2), sigma_y(0.2) {} 

  //Powder target parameters
//  rho_0(13.53), p0(1.0), p_deposition(15000.0), sigma_x(0.3), sigma_y(0.12) {} 


double PowderTarget3DState::pressure(double x, double y, double z) { 
  return p0 + p_deposition*exp(-x*x/(sigma_x*sigma_x) - (y-0.25)*(y-0.25)/(sigma_y*sigma_y) );
}

double PowderTarget3DState::density(double x, double y, double z) {
	return rho_0;
}

void PowderTarget3DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = 0;
	vY = 0;
	vZ = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Gresho3DState
////////////////////////////////////////////////////////////////////////////////////////
