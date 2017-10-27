#include "state_powder_target.h"
#include <iostream>
#include <cmath>


////////////////////////////////////////////////////////////////////////////////////////
// Start of PowderTarget2DState
////////////////////////////////////////////////////////////////////////////////////////

PowderTarget2DState::PowderTarget2DState(): 

  //Thimble parameters
  rho_0(13.53), p0(1.0), p_deposition(15000.0), sigma_x(0.4), sigma_y(0.4) {} 

  //Powder target parameters
//  rho_0(13.53), p0(1.0), p_deposition(15000.0), sigma_x(0.3), sigma_y(0.12) {} 


double PowderTarget2DState::pressure(double x, double y, double z) {
	return p0 + p_deposition*exp(-x*x/(sigma_x*sigma_x) - y*y/(sigma_y*sigma_y));
}

double PowderTarget2DState::density(double x, double y, double z) {
	return rho_0;
}

void PowderTarget2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = 0;
	vY = 0;
	vZ = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Gresho2DState
////////////////////////////////////////////////////////////////////////////////////////
