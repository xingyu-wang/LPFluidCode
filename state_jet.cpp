#include "state_jet.h"
#include <iostream>
#include <cmath>

static const double PI =  atan(1.0)*4.0;

////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet3DState
////////////////////////////////////////////////////////////////////////////////////////

Jet3DState::Jet3DState(): 
m_fDen(13.534), m_fVelX(0), m_fVelY(0), m_fVelZ(0), m_fPressure(0), m_fPMax(110*1000), m_fMagnitude(0.025), m_fAngle(PI/32.) {} 


double Jet3DState::pressure(double x, double y, double z) {

	double xprime = cos(m_fAngle)*x - sin(m_fAngle)*z;
	double zprime = sin(m_fAngle)*x + cos(m_fAngle)*z; 
	double pTmp = m_fPMax*exp(-(y*y + zprime*zprime)/0.0265 - (xprime*xprime)/20.0);
	 
	double tmp = double(rand()) / double(RAND_MAX);
	double ratio = double(rand()) / double(RAND_MAX);	
	if(tmp<=0.5) ratio*=-1;
	//std::cout<<"tmp="<<tmp<<" ratio="<<ratio<<std::endl;
	return pTmp/(1+0.08*cos(10.*x)+0.06*cos(20.*x)+0.07*sin(15.*y))+pTmp*m_fMagnitude*ratio;
}

double Jet3DState::density(double x, double y, double z) {
	return m_fDen;
}

void Jet3DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = m_fVelX;
	vY = m_fVelY;
	vZ = m_fVelZ;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet3DState
////////////////////////////////////////////////////////////////////////////////////////







////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet1DState
////////////////////////////////////////////////////////////////////////////////////////

Jet1DState::Jet1DState(): 
m_fDen(13.534), m_fVelX(0), m_fPressure(0), m_fPCenX(0), m_fPPeak(110*1000), m_fPCoeff(-1) {} 


double Jet1DState::pressure(double x, double y, double z) {
	return m_fPPeak * exp(m_fPCoeff*( (x-m_fPCenX)*(x-m_fPCenX) ) );
}

double Jet1DState::density(double x, double y, double z) {
	return m_fDen;
}

void Jet1DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = m_fVelX;
	vY = 0;
	vZ = 0;
}


Jet1DLeftState::Jet1DLeftState():
m_fDen(2.1161e-3), m_fVelX(1e4), m_fPressure(1.2778e4) {}


double Jet1DLeftState::pressure(double x, double y, double z) {
        return m_fPressure;
}

double Jet1DLeftState::density(double x, double y, double z) {
        return m_fDen;
}

void Jet1DLeftState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        vX = m_fVelX;
        vY = 0;
        vZ = 0;
}

Jet1DRightState::Jet1DRightState():
m_fDen(2.1161e-3), m_fVelX(-1e4), m_fPressure(1.2778e4) {}


double Jet1DRightState::pressure(double x, double y, double z) {
        return m_fPressure;
}

double Jet1DRightState::density(double x, double y, double z) {
        return m_fDen;
}

void Jet1DRightState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        vX = m_fVelX;
        vY = 0;
        vZ = 0;
}

Jet1DCenterState::Jet1DCenterState():
m_fDen(2.1161e-3), m_fVelX(1e4), m_fPressure(1.2778e4) {} //1e4


double Jet1DCenterState::pressure(double x, double y, double z) {
        return m_fPressure;
}

double Jet1DCenterState::density(double x, double y, double z) {
        return m_fDen;
}

void Jet1DCenterState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	if(x<-0.000001)
	        vX = m_fVelX;
	else
		vX = -m_fVelX;
        vY = 0;
        vZ = 0;
}

Jet1DLaterState::Jet1DLaterState():
m_fDen(2.1161e-3), m_fVelX(1e4), m_fPressure(1.2778e4), m_fDenPost(7.3415e-3), m_fPressurePost(3.1008e+05), m_fShockSpeed(4049.7), m_fTime(1e-6) {}


double Jet1DLaterState::pressure(double x, double y, double z) {
	if(fabs(x)<m_fShockSpeed*m_fTime)
		return m_fPressurePost;
	else
	        return m_fPressure;
}

double Jet1DLaterState::density(double x, double y, double z) {
	if(fabs(x)<m_fShockSpeed*m_fTime)
		return m_fDenPost;
	else
	        return m_fDen;
}

void Jet1DLaterState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX=0;
        if(x<-m_fShockSpeed*m_fTime)
                vX = m_fVelX;
        if(x>m_fShockSpeed*m_fTime)
                vX = -m_fVelX;
        vY = 0;
        vZ = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet1DState
////////////////////////////////////////////////////////////////////////////////////////









////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet2DExpansionState
////////////////////////////////////////////////////////////////////////////////////////

Jet2DExpansionState::Jet2DExpansionState(): 
m_fDen(1.327e-6), m_fPressure(4.4873e-2), m_fVelocity(0) {} 


double Jet2DExpansionState::pressure(double x, double y, double z) {
	return m_fPressure;
}

double Jet2DExpansionState::density(double x, double y, double z) {
	return m_fDen;
}

void Jet2DExpansionState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = 0;
	vY = m_fVelocity;
	vZ = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet2DExpansionState
////////////////////////////////////////////////////////////////////////////////////////






////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet3DExpansionState
////////////////////////////////////////////////////////////////////////////////////////

Jet3DExpansionState::Jet3DExpansionState(): 
m_fDen(1.327e-6), m_fPressure(4.4873e-2), m_fVelocity(0) {} 


double Jet3DExpansionState::pressure(double x, double y, double z) {
	return m_fPressure;
}

double Jet3DExpansionState::density(double x, double y, double z) {
	return m_fDen;
}

void Jet3DExpansionState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = 0;
	vY = m_fVelocity;
	vZ = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet3DExpansionState
////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet2DMergeUpperState
////////////////////////////////////////////////////////////////////////////////////////

Jet2DMergeUpperState::Jet2DMergeUpperState(): 
m_fDen(1.327e-6), m_fPressure(4.4873e-2), m_fVelocity(4.0e+3), alpha(-12.0/180.0*M_PI) {} 


double Jet2DMergeUpperState::pressure(double x, double y, double z) {
	return m_fPressure;
}

double Jet2DMergeUpperState::density(double x, double y, double z) {
	return m_fDen;
}

void Jet2DMergeUpperState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = cos(alpha)*m_fVelocity;
	vY = sin(alpha)*m_fVelocity;
	vZ = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet2DMergeUpperState
////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet2DMergeLowerState
////////////////////////////////////////////////////////////////////////////////////////

Jet2DMergeLowerState::Jet2DMergeLowerState(): 
m_fDen(1.327e-6), m_fPressure(4.4873e-2), m_fVelocity(4.0e+3), alpha(12.0/180.0*M_PI) {} 


double Jet2DMergeLowerState::pressure(double x, double y, double z) {
	return m_fPressure;
}

double Jet2DMergeLowerState::density(double x, double y, double z) {
	return m_fDen;
}

void Jet2DMergeLowerState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = cos(alpha)*m_fVelocity;
	vY = sin(alpha)*m_fVelocity;
	vZ = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet2DMergeLowerState
////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet2DCollisionState
////////////////////////////////////////////////////////////////////////////////////////

Jet2DCollisionState::Jet2DCollisionState(): 
m_fDen(1.327e-6), m_fPressure(4.4873e-2), m_fVelocity(8.0e+2) {} 


double Jet2DCollisionState::pressure(double x, double y, double z) {
	return m_fPressure;
}

double Jet2DCollisionState::density(double x, double y, double z) {
	return m_fDen;
}

void Jet2DCollisionState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	if(x<0.0)	vX=m_fVelocity;
	else vX=-m_fVelocity;
	vY = 0;
	vZ = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet2DCollisionState
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
// Start of Jet2DMergeState
////////////////////////////////////////////////////////////////////////////////////////

Jet2DMergeState::Jet2DMergeState(): 
m_fDen(1.327e-6), m_fPressure(4.4873e-2), m_fVelocity(4.0e+3), alpha(-12.0/180.0*M_PI) {} 


double Jet2DMergeState::pressure(double x, double y, double z) {
	return m_fPressure;
}

double Jet2DMergeState::density(double x, double y, double z) {
	return m_fDen;
}

void Jet2DMergeState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = cos(alpha)*m_fVelocity;
	if(y>0)
		vY = sin(alpha)*m_fVelocity;
	else
		vY = -sin(alpha)*m_fVelocity;
	vZ = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Jet2DMergeState
////////////////////////////////////////////////////////////////////////////////////////

