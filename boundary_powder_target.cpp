#include "boundary_powder_target.h"
#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;


////////////////////////////////////////////////////////////////////////////////////////
// Start of PowderTarget2DSolidBoundary
////////////////////////////////////////////////////////////////////////////////////////

PowderTarget2DSolidBoundary::PowderTarget2DSolidBoundary():R(0.6), L(0.2), H(0.6) {}

int PowderTarget2DSolidBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,  
	vector<double>& xb, vector<double>& yb, vector<double>& zb, 
	vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb, vector<double>& vzb) {	

  if (x >= R-L && x <= R && y >= 0 && y <= H) {// right wall
    reflectRight(x, y, pressure, vx, vy, xb, yb, pressureb, vxb, vyb);
    return 1;
  }

  if (x <= -R+L && x >= -R && y >= 0 && y <= H) {// left wall
    reflectLeft(x, y, pressure, vx, vy, xb, yb, pressureb, vxb, vyb);
    return 1;
  }
	
  double dist = sqrt(x*x+y*y);

  if (dist >= R-L && dist <= R && y < 0) {//low semicircle

    double factor = (2.*R-dist)/dist;
    double normal_vx = x/dist;
    double normal_vy = y/dist;
    
    xb.push_back(factor*x);
    yb.push_back(factor*y);
    pressureb.push_back(pressure);
    double ip = ip2d(vx,vy,normal_vx,normal_vy);
    if(ip<=0) { // leaving solid boundary
      vxb.push_back(0);
      vyb.push_back(0);	
    }
    else {
      vxb.push_back(vx-2.*ip*normal_vx);
      vyb.push_back(vy-2.*ip*normal_vy);
    }
    return 1;	
  }
  
  return 0; //outside the container
}


void PowderTarget2DSolidBoundary::reflectRight(double x, double y, double pressure, double vx, double vy,  
	vector<double>& xb, vector<double>& yb, vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb) {
	
        xb.push_back(2*R-x);
	yb.push_back(y);
	pressureb.push_back(pressure);
	if(ip2d(vx,vy,1,0)<=0) { // approaching solid boundary
		vxb.push_back(0);
		vyb.push_back(0);	
	}
	else {
		vxb.push_back(-vx);
		vyb.push_back(vy);
	}

}

void PowderTarget2DSolidBoundary::reflectLeft(double x, double y, double pressure, double vx, double vy,  
	vector<double>& xb, vector<double>& yb, vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb) {
	
        xb.push_back(-2*R-x);
	yb.push_back(y);
	pressureb.push_back(pressure);
	if(ip2d(vx,vy,-1,0)<=0) { // approaching solid boundary
		vxb.push_back(0);
		vyb.push_back(0);	
	}
	else {
		vxb.push_back(-vx);
		vyb.push_back(vy);
	}

}


////////////////////////////////////////////////////////////////////////////////////////
// End of PowderTarget2DSolidBoundary
////////////////////////////////////////////////////////////////////////////////////////
