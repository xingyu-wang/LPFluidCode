#include "state_ballexp.h"
#include <iostream>
#include <cmath>

//static const double PI =  atan(1.0)*4.0;


////////////////////////////////////////////////////////////////////////////////////////
// Start of Ballexp3DState
////////////////////////////////////////////////////////////////////////////////////////

Ballexp3DState::Ballexp3DState(): 
m_fDen(1.327e-6), m_fPressure(4.4873e-2) {} 


double Ballexp3DState::pressure(double x, double y, double z) {
	return m_fPressure;
}

double Ballexp3DState::density(double x, double y, double z) {
	return m_fDen;
}

void Ballexp3DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = 0;
	vY = 0;
	vZ = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Ballexp3DState
////////////////////////////////////////////////////////////////////////////////////////

Ballrotate3DState::Ballrotate3DState():
m_fDen(1), m_fPressure(0) {}


double Ballrotate3DState::pressure(double x, double y, double z) {
        return m_fPressure;
}

double Ballrotate3DState::density(double x, double y, double z) {
        return m_fDen;
}

void Ballrotate3DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        vX = -5.0*y;
        vY = 5.0*x;
        vZ = 0;
}

Ballpressurewave3DState::Ballpressurewave3DState():
m_fDen(1), m_fPressure(100) {}


double Ballpressurewave3DState::pressure(double x, double y, double z) {
	if ((x*x+y*y+z*z)<25.0)
	        return m_fPressure;
	else
		return 0.0;
}

double Ballpressurewave3DState::density(double x, double y, double z) {
        return m_fDen;
}

void Ballpressurewave3DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        vX = 0;
        vY = 0;
        vZ = 0;
}

