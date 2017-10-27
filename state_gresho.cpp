#include "state_gresho.h"
#include <iostream>
#include <cmath>

//static const double PI =  atan(1.0)*4.0;

////////////////////////////////////////////////////////////////////////////////////////
// Start of Gresho2DState
////////////////////////////////////////////////////////////////////////////////////////

Gresho2DState::Gresho2DState(): 
m_fDensity(1) {} 


double Gresho2DState::pressure(double x, double y, double z) {
	double val;
	double r = sqrt(x*x + y*y + z*z);

	if (r <= 0.2) val = 5.0 + 12.5*r*r;
	else if  (r > 0.2 && r <= 0.4) val = 9.0 - 4*log(0.2) + 12.5*r*r - 20*r + 4*log(r);
	else val = 3.0 + 4*log(2.0);
	return val;
}

double Gresho2DState::density(double x, double y, double z) {
	return m_fDensity;
}

void Gresho2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	double r = sqrt(x*x + y*y + z*z);
	double th = atan2(y,x);
	double u;

	if (r <= 0.2) 
	  {
	    u = 5*r;
	    vX = -u*sin(th);
	    vY = u*cos(th);
	  }
	else if  (r > 0.2 && r <= 0.4) 
	  {
	    u = 2.0 - 5*r;
	    vX = -u*sin(th);
	    vY = u*cos(th);
	  }	    
	else 
	  {
	    vX = 0;
	    vY = 0;
	  }
	vZ = 0;
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Gresho2DState
////////////////////////////////////////////////////////////////////////////////////////
Yee2DState::Yee2DState(){}

double Yee2DState::pressure(double x, double y, double z) {
	double val=0;
	double r2=x*x+y*y;
	double gamma=1.4;
	double beta=5.0;
	double T=1-(gamma-1)*beta*beta/8.0/gamma/M_PI/M_PI*exp(1.0-r2);
	val=pow(T,gamma/(gamma-1.0));

        return val;
}

double Yee2DState::density(double x, double y, double z) {
        double val=0;
        double r2=x*x+y*y;
        double gamma=1.4;
        double beta=5.0;
        double T=1-(gamma-1)*beta*beta/8.0/gamma/M_PI/M_PI*exp(1.0-r2);
        val=pow(T,1.0/(gamma-1.0));

        return val;
}

void Yee2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        double r2=x*x+y*y;
        double beta=5.0;
	vX=beta/2.0/M_PI*exp((1.0-r2)/2.0)*(-y);
        vY=beta/2.0/M_PI*exp((1.0-r2)/2.0)*x;
	vZ=0.0;

}

Sedov2DState::Sedov2DState(){}

double Sedov2DState::pressure(double x, double y, double z) {
/*        double val=0;
        double r2=x*x+y*y;
        double gamma=1.666667;
        double e_0=1.0;
        double dx=0.01;
        if(r2<(0.5*dx*dx))
                val=val+(gamma-1)*e_0/(sqrt(3.0)*0.5*dx*dx);

        return val;*/

        double val=1.0;
        double r=sqrt(x*x+y*y);
        double gamma=1.666667;
        double e_0=1.0;
	double dx=0.01;
	double h=10*dx;
	double u=r/h;
	double c=0.212224/17.1539;
//	double c=0.212224/4.28849;
	double p_0=(gamma-1)*e_0/(M_PI*sqrt(3.0)*0.5*dx*dx);
	if(u<1.0)
		val=val+c*(1-1.5*u*u+0.75*u*u*u)*p_0;
	else if(u<2.0)
		val=val+c*0.25*(2-u)*(2-u)*(2-u)*p_0;

        return val;

/*	double val=0;
	double r=sqrt(x*x+y*y);
	double c=0.212224/0.212224;
//	double sigma=0.25;
	val=c/(r*r+1e-3);
//	val=c*exp(-r*r/sigma/sigma/2);	
	return val;*/
}

double Sedov2DState::density(double x, double y, double z) {
        return 1.0;
}

void Sedov2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        vX=0.0;
        vY=0.0;
        vZ=0.0;

}

Noh2DState::Noh2DState(){}

double Noh2DState::pressure(double x, double y, double z) {
        return 1e-6;
}

double Noh2DState::density(double x, double y, double z) {
        return 1.0;
}

void Noh2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
        double r=sqrt(x*x+y*y);
	vX=-x/r;
        vY=-y/r;
	vZ=0.0;

}

////////////////////////////////////////////////////////////////////////////////////////
//// Start of ConvergentShock2DState
//////////////////////////////////////////////////////////////////////////////////////////

ConvergentShock2DState::ConvergentShock2DState():
m_fDensity(1), m_fPressure(1e-6), m_fVelocity(1){}


double ConvergentShock2DState::pressure(double x, double y, double z) {
        return m_fPressure;
}

double ConvergentShock2DState::density(double x, double y, double z) {
        return m_fDensity;
}

void ConvergentShock2DState::velocity(double x, double y, double z, double& vX, double& vY, double& vZ) {
	double r=sqrt(x*x+y*y);
	if (r>0)
	{
	        vX = -x/r*m_fVelocity;
        	vY = -y/r*m_fVelocity;
        	vZ = 0;
	}
	else
	{
		vX = 0;
		vY = 0;
		vZ = 0;
	}
}

////////////////////////////////////////////////////////////////////////////////////////
//// End of ConvergentShock2DState
//////////////////////////////////////////////////////////////////////////////////////////
