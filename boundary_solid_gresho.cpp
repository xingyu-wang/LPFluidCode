#include "boundary_solid_gresho.h"
#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;


////////////////////////////////////////////////////////////////////////////////////////
// Start of Gresho2DSolidBoundary
////////////////////////////////////////////////////////////////////////////////////////

Gresho2DSolidBoundary::Gresho2DSolidBoundary():radius(1.), thickness(0.3) {
	bo = radius-thickness;
}

int Gresho2DSolidBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,  
	vector<double>& xb, vector<double>& yb, vector<double>& zb, 
	vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb, vector<double>& vzb) {	
	
	double dist = sqrt(x*x+y*y);

	if(dist < bo) return 0; // inside
	
	if(dist > radius) return 0; // outside	
	
	if(dist==0) return 0; // origin (center of circle)

	double factor = (2.*radius-dist)/dist;
//	double normal_vx = x/dist;
//	double normal_vy = y/dist;

	xb.push_back(factor*x);
	yb.push_back(factor*y);
//	pressureb.push_back(pressure);
	pressureb.push_back(3+4*log(2));
	vxb.push_back(0);
	vyb.push_back(0);
//	double ip = ip2d(vx,vy,normal_vx,normal_vy);
/*	if(ip<=0) { // leaving solid boundary
		vxb.push_back(0);
		vyb.push_back(0);	
	}
	else {
		vxb.push_back(vx-2.*ip*normal_vx);
		vyb.push_back(vy-2.*ip*normal_vy);
	}
*/	
	return 1;	

}

Sedov2DSolidBoundary::Sedov2DSolidBoundary():radius(1.), thickness(0.3) {
        bo = radius-thickness;
}

int Sedov2DSolidBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,  
	vector<double>& xb, vector<double>& yb, vector<double>& zb, 
	vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb, vector<double>& vzb) {	
	
	double dist = sqrt(x*x+y*y);

	if(dist < bo) return 0; // inside
	
	if(dist > radius) return 0; // outside	
	
	if(dist==0) return 0; // origin (center of circle)

	double factor = (2.*radius-dist)/dist;
//	double normal_vx = x/dist;
//	double normal_vy = y/dist;

	xb.push_back(factor*x);
	yb.push_back(factor*y);
	pressureb.push_back(pressure);
//	pressureb.push_back(3+4*log(2));
	vxb.push_back(0);
	vyb.push_back(0);
//	double ip = ip2d(vx,vy,normal_vx,normal_vy);
/*	if(ip<=0) { // leaving solid boundary
		vxb.push_back(0);
		vyb.push_back(0);	
	}
	else {
		vxb.push_back(vx-2.*ip*normal_vx);
		vyb.push_back(vy-2.*ip*normal_vy);
	}
*/	
	return 1;	
}

////////////////////////////////////////////////////////////////////////////////////////
// End of Gresho2DSolidBoundary
Yee2DSolidBoundary::Yee2DSolidBoundary():radius(5.), thickness(3.0) {
        bo = radius-thickness;
}

int Yee2DSolidBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,
        vector<double>& xb, vector<double>& yb, vector<double>& zb,
        vector<double>& pressureb, vector<double>& vxb, vector<double>& vyb, vector<double>& vzb) {

        double dist = sqrt(x*x+y*y);

        if(dist < bo) return 0; // inside

        if(dist > radius) return 0; // outside  

        if(dist==0) return 0; // origin (center of circle)

        double factor = (2.*radius-dist)/dist;
//        double normal_vx = x/dist;
//        double normal_vy = y/dist;

        xb.push_back(factor*x);
        yb.push_back(factor*y);

	double r=2.*radius-dist;
	
//      pressureb.push_back(pressure);
        pressureb.push_back(1.0);
        vxb.push_back(2.5/M_PI*exp(0.5-0.5*r*r)*(-factor*y));
        vyb.push_back(2.5/M_PI*exp(0.5-0.5*r*r)*(factor*x));
//      double ip = ip2d(vx,vy,normal_vx,normal_vy);
/*      if(ip<=0) { // leaving solid boundary
                vxb.push_back(0);
                vyb.push_back(0);       
        }
        else {
                vxb.push_back(vx-2.*ip*normal_vx);
                vyb.push_back(vy-2.*ip*normal_vy);
        }
*/
        return 1;

}
///////////////////////////////////////////////////////////////////////////////////////
