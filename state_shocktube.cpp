#include "state_shocktube.h"
#include <iostream>
#include <cmath>

//static const double PI =  atan(1.0)*4.0;

////////////////////////////////////////////////////////////////////////////////////////
// Start of Shocktube2DState
////////////////////////////////////////////////////////////////////////////////////////

Shocktube2DState::Shocktube2DState(): 
m_fDenL(1.0), m_fDenR(1.0), m_fPressureL(1.0), m_fPressureR(0.1) {} 


double Shocktube2DState::pressure(double x, double y, double z) {
	if(x<=0) return m_fPressureL;
	else return m_fPressureR;
		
}

double Shocktube2DState::density(double x, double y, double z) {
	if(x<=0) return m_fDenL;
	else return m_fDenR;
}

void Shocktube2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	vX = 0;
	vY = 0;
	vZ = 0;
}

TPShocktube2DState::TPShocktube2DState():
m_fDen1(1.0), m_fDen2(0.125), m_fDen3(1.0), m_fPressure1(1.0), m_fPressure2(0.1), m_fPressure3(0.1) {}


double TPShocktube2DState::pressure(double x, double y, double z) {
        if(x<=1) return m_fPressure1;
        else if(y>1.5) return m_fPressure2;
	else return m_fPressure3;

}

double TPShocktube2DState::density(double x, double y, double z) {
        if(x<=1) return m_fDen1;
        else if(y>1.5) return m_fDen2;
	else return m_fDen3;
}

void TPShocktube2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        vX = 0;
        vY = 0;
        vZ = 0;
}


////////////////////////////////////////////////////////////////////////////////////////
// End of Shocktube2DState
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
// Start of NormalShock2DState
////////////////////////////////////////////////////////////////////////////////////////

NormalShock2DState::NormalShock2DState(): 
m_fDenL(1.327e-6), m_fDenR(1.327e-6), m_fPressureL(4.4873e-2), m_fPressureR(4.4873e-2), m_fVelocityL(8.3165e01), m_fVelocityR(-8.3165e01) {} 


double NormalShock2DState::pressure(double x, double y, double z) {
	if(x<=0) return m_fPressureL;
	else return m_fPressureR;
		
}

double NormalShock2DState::density(double x, double y, double z) {
	if(x<=0) return m_fDenL;
	else return m_fDenR;
}

void NormalShock2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	if(x<=0) vX=m_fVelocityL;
	else	vX=m_fVelocityR;
	vY = 0;
	vZ = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of NormalShock2DState
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
// Start of SimpleWave2DState
////////////////////////////////////////////////////////////////////////////////////////

SimpleWave2DState::SimpleWave2DState():
m_fDenL(1.0), m_fPressureL(1.0), m_fVPost(0.927452506899452), m_fGamma(1.4), m_fITime(3.0) {}


double SimpleWave2DState::pressure(double x, double y, double z) {

        double gamma=m_fGamma;
        double rho_left=m_fDenL;
        double p_left=m_fPressureL;
        double v_post=m_fVPost;
	double time=m_fITime;

//        double u_left=0.0;
        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);
        double x1=-c_left*time;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time;

        if (x<x1)
		return p_left;
        else if (x<x2)
        {
                double c=mu*mu*(-x/time)+(1-mu*mu)*c_left;
                double rho=rho_left*pow((c/c_left),(2/(gamma-1)));
		return p_left*pow(rho/rho_left,gamma);
        }
        else
	{
                double c=mu*mu*(-x2/time)+(1-mu*mu)*c_left;
                double rho=rho_left*pow((c/c_left),(2/(gamma-1)));
                return p_left*pow(rho/rho_left,gamma);
	}

}

double SimpleWave2DState::density(double x, double y, double z) {

        double gamma=m_fGamma;
        double rho_left=m_fDenL;
        double p_left=m_fPressureL;
        double v_post=m_fVPost;
        double time=m_fITime;

//        double u_left=0.0;
        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);
        double x1=-c_left*time;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time;

        if (x<x1)
                return rho_left;
        else if (x<x2)
        {
                double c=mu*mu*(-x/time)+(1-mu*mu)*c_left;
                return rho_left*pow((c/c_left),(2/(gamma-1)));
        }
        else 
        {
                double c=mu*mu*(-x2/time)+(1-mu*mu)*c_left;
                return rho_left*pow((c/c_left),(2/(gamma-1)));
        }

}

void SimpleWave2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {

        double gamma=m_fGamma;
        double rho_left=m_fDenL;
        double p_left=m_fPressureL;
        double v_post=m_fVPost;
        double time=m_fITime;

        double u_left=0.0;
        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);
        double x1=-c_left*time;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time;

        if (x<x1)
                vX =  u_left;
        else if (x<x2)
        {
		vX = (1-mu*mu)*(x/time+c_left);
        }
        else
        {
		vX = (1-mu*mu)*(x2/time+c_left);
        }

        vY = 0;
        vZ = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of SimpleWave2DState
////////////////////////////////////////////////////////////////////////////////////////


SodShocktube2DState::SodShocktube2DState():
m_fDenL(1.0), m_fDenR(0.125), m_fPressureL(1.0), m_fPressureR(0.1) {}


double SodShocktube2DState::pressure(double x, double y, double z) {
        if(x<=-0.00001) return m_fPressureL;
        else return m_fPressureR;
//	else return (m_fPressureL+m_fPressureR*exp(x/0.0565))/(1+exp(x/0.0565));
}

double SodShocktube2DState::density(double x, double y, double z) {
        if(x<=-0.00001) return m_fDenL;
        else return m_fDenR;
//        else return (m_fDenL+m_fDenR*exp(x/0.0565))/(1+exp(x/0.0565));
}

void SodShocktube2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        vX = 0;
        vY = 0;
        vZ = 0;
}

SodShocktube2DLaterState::SodShocktube2DLaterState():
m_fDenL(1.0), m_fDenR(0.125), m_fPressureL(1.0), m_fPressureR(0.1), time(0.5) {}


double SodShocktube2DLaterState::pressure(double x, double y, double z) {
//      x=x-0.005;
        double gamma=1.4;
        double rho_left=1.0;
        double p_left=1.0;
//        double u_left=0.0;

//        double rho_right=0.125;
        double p_right=0.1;
//        double u_right=0.0;

        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);

//        double rho_middle=0.42632;
//        double rho_post=0.26557;
        double p_post=0.30313;
        double v_post=0.92745;
        double v_shock=1.75216;

        double x1=-c_left*time;
        double x3=v_post*time;
        double x4=v_shock*time;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time;

        if(x<x1)
                return p_left;
        if(x<x2)
        {
                double c=mu*mu*(-x/time)+(1-mu*mu)*c_left;
                double rho=rho_left*pow((c/c_left),(2/(gamma-1)));
                return p_left*pow(rho/rho_left,gamma);
        }
        if(x<x3)
                return p_post;
        if(x<x4)
                return p_post;
        return p_right;

}

double SodShocktube2DLaterState::density(double x, double y, double z) {

//      x=x-0.005;
        double gamma=1.4;
        double rho_left=1.0;
        double p_left=1.0;
//        double u_left=0.0;

        double rho_right=0.125;
//        double p_right=0.1;
//        double u_right=0.0;

        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);

        double rho_middle=0.42632;
        double rho_post=0.26557;
//        double p_post=0.30313;
        double v_post=0.92745;
        double v_shock=1.75216;

        double x1=-c_left*time;
        double x3=v_post*time;
        double x4=v_shock*time;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time;
        if(x<x1)
                return rho_left;
        if(x<x2)
        {
                double c=mu*mu*(-x/time)+(1-mu*mu)*c_left;
                return rho_left*pow((c/c_left),(2/(gamma-1)));
        }
        if(x<x3)
                return rho_middle;
        if(x<x4)
                return rho_post;
        return rho_right;
}

void SodShocktube2DLaterState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
//      x=x-0.005;
        double gamma=1.4;
        double rho_left=1.0;
        double p_left=1.0;
        double u_left=0.0;

//        double rho_right=0.125;
//        double p_right=0.1;
        double u_right=0.0;

        double mu=sqrt((gamma-1)/(gamma+1));
        double c_left=sqrt(gamma*p_left/rho_left);

//        double rho_middle=0.42632;
//        double rho_post=0.26557;
//        double p_post=0.30313;
        double v_post=0.92745;
        double v_shock=1.75216;

        double x1=-c_left*time;
        double x3=v_post*time;
        double x4=v_shock*time;
        double c2=c_left-(gamma-1)/2*v_post;
        double x2=(v_post-c2)*time;
        if(x<x1)
                vX=u_left;
        else if(x<x2)
                vX=(1-mu*mu)*(x/time+c_left);
        else if(x<x3)
                vX=v_post;
        else if(x<x4)
                vX=v_post;
        else vX=u_right;

	vY=0.0;
	vZ=0.0;
}


RayleighTaylor2DState::RayleighTaylor2DState():
m_fDenU(2.0), m_fDenL(1.0), m_fPressure(1.0), m_fWaveLength(1.0/3.0), m_fMagnitude(0.01), m_fG(-0.1) {}


double RayleighTaylor2DState::pressure(double x, double y, double z) {
	if (y<0) return m_fPressure+m_fG*m_fDenL*y;
        return m_fPressure+m_fG*m_fDenU*y;
}

double RayleighTaylor2DState::density(double x, double y, double z) {
        if(y<=m_fMagnitude*cos(2*M_PI*x/m_fWaveLength)) return m_fDenL;
        else return m_fDenU;
}

void RayleighTaylor2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        vX = 0;
        vY = 0;
        vZ = 0;
}

RayleighTaylor3DState::RayleighTaylor3DState():
m_fDenU(2.0), m_fDenL(1.0), m_fPressure(1.0), m_fWaveLengthX(1.0/3.0), m_fWaveLengthY(1.0/3.0), m_fMagnitude(0.01), m_fG(-0.1) {}


double RayleighTaylor3DState::pressure(double x, double y, double z) {
        if (z<0) return m_fPressure+m_fG*m_fDenL*z;
        return m_fPressure+m_fG*m_fDenU*z;
}

double RayleighTaylor3DState::density(double x, double y, double z) {
        if(z<=m_fMagnitude*cos(2*M_PI*x/m_fWaveLengthX)*cos(2*M_PI*y/m_fWaveLengthY)) return m_fDenL;
        else return m_fDenU;
}

void RayleighTaylor3DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        vX = 0;
        vY = 0;
        vZ = 0;
}

DamBreak2DState::DamBreak2DState():
m_fDen(1000.), m_fPressure(0),  m_fYLength(0.5), m_fG(-10.0) {}


double DamBreak2DState::pressure(double x, double y, double z) {
        return m_fPressure+m_fG*m_fDen*(y-m_fYLength);
}

double DamBreak2DState::density(double x, double y, double z) {
        return m_fDen;
}

void DamBreak2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        vX = 0;
        vY = 0;
        vZ = 0;
}

BoundaryTest2DState::BoundaryTest2DState():
m_fDen(1000.), m_fPressure(0),  m_fYLength(0.1), m_fG(-10.0) {}


double BoundaryTest2DState::pressure(double x, double y, double z) {
        return m_fPressure+m_fG*m_fDen*(y-m_fYLength);
}

double BoundaryTest2DState::density(double x, double y, double z) {
        return m_fDen;
}

void BoundaryTest2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        vX = 1.;
        vY = 0;
        vZ = 0;
}

KelvinHelmholtz2DState::KelvinHelmholtz2DState():
m_fDen1(1.0), m_fDen2(2.0), m_fVel1(0.5), m_fVel2(-0.5), m_fDelta(0.025), m_fDeltay(0.01), m_fPressure(2.5), m_fWaveLength(0.5){}


double KelvinHelmholtz2DState::pressure(double x, double y, double z) {
        return m_fPressure;
}

double KelvinHelmholtz2DState::density(double x, double y, double z) {
	double my=y;
	double m_fDenm=0.5*(m_fDen1-m_fDen2);
	if(my<0.25) return m_fDen1-m_fDenm*exp((my-0.25)/m_fDelta);
        else if(my<0.5) return m_fDen2+m_fDenm*exp((0.25-my)/m_fDelta);
        else if(my<0.75) return m_fDen2+m_fDenm*exp((my-0.75)/m_fDelta);
        else return m_fDen1-m_fDenm*exp((0.75-my)/m_fDelta);
}

void KelvinHelmholtz2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	double my=y;
	double mx=x;
        double m_fVelm=0.5*(m_fVel1-m_fVel2);
        if(my<0.25) vX = m_fVel1-m_fVelm*exp((my-0.25)/m_fDelta);
        else if(my<0.5) vX = m_fVel2+m_fVelm*exp((0.25-my)/m_fDelta);
        else if(my<0.75) vX = m_fVel2+m_fVelm*exp((my-0.75)/m_fDelta);
        else vX = m_fVel1-m_fVelm*exp((0.75-my)/m_fDelta);
	
	vY = m_fDeltay*sin(2.0*M_PI*mx/m_fWaveLength);
	vZ = 0.0;
}

