#include "boundary_solid_tpshocktube.h"
#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;
static const double onesq2 = 1.0/std::sqrt(2.0);

////////////////////////////////////////////////////////////////////////////////////////
// Start of Shocktube2DSolidBoundary
////////////////////////////////////////////////////////////////////////////////////////

TPShocktube2DSolidBoundary::TPShocktube2DSolidBoundary():lengthX(7), lengthY(3), thickness(0.8) {
	rb = lengthX-thickness;
	lb = thickness;
	nb = lengthY-thickness;
	sb = thickness;
	
	rbo = lengthX;
	lbo = 0;
	nbo = lengthY;
	sbo = 0;
	
	epsilon=1e-3;
}

int TPShocktube2DSolidBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,  
	vector<double>& xb, vector<double>& yb, vector<double>& zb, 
	vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb, vector<double>& vzb) {	
	
	if(x<rb && x>lb && y<nb && y>sb) {
		//cout<<"x="<<x<<" y="<<y<<" z="<<z<<endl;
		//cout<<"rb="<<rb<<" lb="<<lb<<" nb="<<nb<<" sb="<<sb<<endl;
		return 0; // inside
	}
	if(x>rbo || x<lbo || y>nbo || y<sbo) return 0; // outside	

	
	if(x>=rb && y>=nb) { // right north
		reflectRight(x, y, pressure, vx, vy,  
					 xb, yb, pressureb, vxb, vyb);
		reflectNorth(x, y, pressure, vx, vy,  
					 xb, yb, pressureb, vxb, vyb);
		reflectRightNorth(x, y, pressure, vx, vy,  
					      xb, yb, pressureb, vxb, vyb);
		return 3;
	}
	else if(x>=rb && y<=sb) { // right south
		reflectRight(x, y, pressure, vx, vy,  
					 xb, yb, pressureb, vxb, vyb);
		reflectSouth(x, y, pressure, vx, vy,  
					 xb, yb, pressureb, vxb, vyb);
		reflectRightSouth(x, y, pressure, vx, vy,  
					      xb, yb, pressureb, vxb, vyb);
		return 3;
	}
	else if(x<=lb && y>=nb) { // left north
		reflectLeft(x, y, pressure, vx, vy,  
					 xb, yb, pressureb, vxb, vyb);
		reflectNorth(x, y, pressure, vx, vy,  
					 xb, yb, pressureb, vxb, vyb);
		reflectLeftNorth(x, y, pressure, vx, vy,  
					      xb, yb, pressureb, vxb, vyb);
		return 3;
	}
	else if(x<=lb && y<=sb) { // left south
		reflectLeft(x, y, pressure, vx, vy,  
					 xb, yb, pressureb, vxb, vyb);
		reflectSouth(x, y, pressure, vx, vy,  
					 xb, yb, pressureb, vxb, vyb);
		reflectLeftSouth(x, y, pressure, vx, vy,  
					      xb, yb, pressureb, vxb, vyb);
		return 3;
	}
	else if(x>=rb) { // right
		reflectRight(x, y, pressure, vx, vy,  
					 xb, yb, pressureb, vxb, vyb);
		return 1;
	}
	else if(x<=lb) { // left
		reflectLeft(x, y, pressure, vx, vy,  
					 xb, yb, pressureb, vxb, vyb);
		return 1;
	}
	else if(y>=nb) { // north
		reflectNorth(x, y, pressure, vx, vy,  
					 xb, yb, pressureb, vxb, vyb);
		return 1;
	}
	else if(y<=sb) { // south
		reflectSouth(x, y, pressure, vx, vy,  
					 xb, yb, pressureb, vxb, vyb);
		return 1;
	}
	else {
		assert(false);
	}

}

 
void TPShocktube2DSolidBoundary::reflectRight(double x, double y, double pressure, double vx, double vy,  
	vector<double>& xb, vector<double>& yb, vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb) {
	
	xb.push_back(2*rbo-x);
	yb.push_back(y);
	pressureb.push_back(pressure);
	if((ip2d(vx,vy,1,0)+epsilon)<=0) { // approaching solid boundary
		vxb.push_back(0);
		vyb.push_back(0);	
	}
	else {
		vxb.push_back(-vx);
		vyb.push_back(vy);
	}

}

void TPShocktube2DSolidBoundary::reflectLeft(double x, double y, double pressure, double vx, double vy,  
	vector<double>& xb, vector<double>& yb, vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb) {
	
	xb.push_back(2*lbo-x);
	yb.push_back(y);
	pressureb.push_back(pressure);
	if((ip2d(vx,vy,-1,0)+epsilon)<=0) { // approaching solid boundary
		vxb.push_back(0);
		vyb.push_back(0);	
	}
	else {
		vxb.push_back(-vx);
		vyb.push_back(vy);
	}

}


void TPShocktube2DSolidBoundary::reflectNorth(double x, double y, double pressure, double vx, double vy,  
	vector<double>& xb, vector<double>& yb, vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb) {
	
	xb.push_back(x);
	yb.push_back(2*nbo-y);
	pressureb.push_back(pressure);
	if((ip2d(vx,vy,0,1)+epsilon)<=0) { // approaching solid boundary
		vxb.push_back(0);
		vyb.push_back(0);	
	}
	else {
		vxb.push_back(vx);
		vyb.push_back(-vy);
	}

}


void TPShocktube2DSolidBoundary::reflectSouth(double x, double y, double pressure, double vx, double vy,  
	vector<double>& xb, vector<double>& yb, vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb) {
	
	xb.push_back(x);
	yb.push_back(2*sbo-y);
	pressureb.push_back(pressure);
	if((ip2d(vx,vy,0,-1)+epsilon)<=0) { // approaching solid boundary
		vxb.push_back(0);
		vyb.push_back(0);	
	}
	else {
		vxb.push_back(vx);
		vyb.push_back(-vy);
	}

}

void TPShocktube2DSolidBoundary::reflectRightNorth(double x, double y, double pressure, double vx, double vy,  
	vector<double>& xb, vector<double>& yb, vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb) {
		
	xb.push_back(2*rbo-x);
	yb.push_back(2*nbo-y);
	pressureb.push_back(pressure);
	double ip = ip2d(vx,vy,onesq2,onesq2);
	if((ip+epsilon)<=0) { // approaching solid boundary
		vxb.push_back(0);
		vyb.push_back(0);
	}
	else {
		vxb.push_back(vx-2.*ip*onesq2);
		vyb.push_back(vy-2.*ip*onesq2);
	}

}


void TPShocktube2DSolidBoundary::reflectRightSouth(double x, double y, double pressure, double vx, double vy,  
	vector<double>& xb, vector<double>& yb, vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb) {
		
	xb.push_back(2*rbo-x);
	yb.push_back(2*sbo-y);
	pressureb.push_back(pressure);
	double ip = ip2d(vx,vy,onesq2,-onesq2);
	if((ip+epsilon)<=0) { // approaching solid boundary
		vxb.push_back(0);
		vyb.push_back(0);
	}
	else {
		vxb.push_back(vx-2.*ip*onesq2);
		vyb.push_back(vy-2.*ip*(-1)*onesq2);
	}

}


void TPShocktube2DSolidBoundary::reflectLeftNorth(double x, double y, double pressure, double vx, double vy,  
	vector<double>& xb, vector<double>& yb, vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb) {
		
	xb.push_back(2*lbo-x);
	yb.push_back(2*nbo-y);
	pressureb.push_back(pressure);
	double ip = ip2d(vx,vy,-onesq2,onesq2);
	if((ip+epsilon)<=0) { // approaching solid boundary
		vxb.push_back(0);
		vyb.push_back(0);
	}
	else {
		vxb.push_back(vx-2.*ip*(-1)*onesq2);
		vyb.push_back(vy-2.*ip*onesq2);
	}

}


void TPShocktube2DSolidBoundary::reflectLeftSouth(double x, double y, double pressure, double vx, double vy,  
	vector<double>& xb, vector<double>& yb, vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb) {
		
	xb.push_back(2*lbo-x);
	yb.push_back(2*sbo-y);
	pressureb.push_back(pressure);
	double ip = ip2d(vx,vy,-onesq2,-onesq2);
	if((ip+epsilon)<=0) { // approaching solid boundary
		vxb.push_back(0);
		vyb.push_back(0);
	}
	else {
		vxb.push_back(vx-2.*ip*(-1)*onesq2);
		vyb.push_back(vy-2.*ip*(-1)*onesq2);
	}

}

////////////////////////////////////////////////////////////////////////////////////////
// End of Shocktube2DSolidBoundary
////////////////////////////////////////////////////////////////////////////////////////
