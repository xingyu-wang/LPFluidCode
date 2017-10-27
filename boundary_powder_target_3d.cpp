#include "boundary_powder_target_3d.h"
#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;


////////////////////////////////////////////////////////////////////////////////////////
// Start of PowderTarget3DSolidBoundary
////////////////////////////////////////////////////////////////////////////////////////

PowderTarget3DSolidBoundary::PowderTarget3DSolidBoundary():R(0.75), L(0.2), H(0.75),length(4.0) {} // L is thickness

int PowderTarget3DSolidBoundary::operator()(double x, double y, double z, double pressure, double vx, double vy, double vz,  
	vector<double>& xb, vector<double>& yb, vector<double>& zb, 
	vector<double>& pressureb, vector<double>& vxb,	vector<double>& vyb, vector<double>& vzb) {	

  double dist = sqrt(x*x + y*y);
  
  
  if (x >= R-L && x <= R && y >= 0 && y <= H && fabs(z) <= length/2.0-L) {// right wall
    reflectRight(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    return 1;
  }
  
  
  if (x <= -R+L && x >= -R && y >= 0 && y <= H && fabs(z) <= length/2.0-L) {// left wall
    reflectLeft(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    return 1;
  }
  
  if ((x > -R+L && x < R-L && y>=0 && y<=H && z> length/2.0 -L && z< length/2.0) ||( dist <=R-L  && y< 0 &&z> length/2.0 -L && z< length/2.0)) {// front wall
    reflectFront(x,y,z,pressure,vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    return 1;
  }
 
  if ((x > -R+L && x < R-L && y>=0 && y<=H && z< -length/2.0 +L && z>-length/2.0)||( dist <=R-L && y< 0 &&z<= -length/2.0 +L && z>=- length/2.0)){// back wall
    reflectBack(x,y,z,pressure,vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    return 1;
  }
 
  
  if (dist >= R-L && dist <= R && y <= 0 && fabs(z) <= length/2.0-L) {//low semicircle
    reflectSemiCircle(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    return 1;	
  }
  
  if (x >= R-L && x <= R && y >= 0 && y <= H && z > length/2.0-L && z <=length/2.0 ) {// right front wall ( 2 boundaries)
    reflectRight(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    reflectFront(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    reflectRightFront(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    return 3;
  }
  
  if (x <= -R+L && x >= -R && y >= 0 && y <= H && z > length/2.0-L && z <=length/2.0) {// left front wall
    reflectLeft(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    reflectFront(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    reflectLeftFront(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    return 3;
  }
  
  
  if  ( dist>=R-L && dist<=R && y < 0 &&  z > length/2.0-L && z <=length/2.0) {//low semicircle, front boundary
    reflectFront(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    reflectSemiCircle(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    reflectSemiCircleFront(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    return 3;	
  }
  
  
  if (x >= R-L && x <= R && y >= 0 && y <= H && z <= -length/2.0+L  && z>=-length/2.0) {// right back wall
    reflectRight(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    reflectBack(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    reflectRightBack(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    return 3;
  }

  if (x <= -R+L && x >= -R && y >= 0 && y <= H && z <= -length/2.0+L  && z>=-length/2.0) {// left back wall
    reflectLeft(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    reflectBack(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    reflectLeftBack(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    return 3;
  }
  
  
  if ( dist>=R-L && dist<=R && y < 0 && z <= -length/2.0+L  && z>=-length/2.0) {//low semicircle, back boundaries( z<0)
    reflectBack(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    reflectSemiCircle(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    reflectSemiCircleBack(x, y, z,pressure, vx, vy, vz, xb, yb, zb, pressureb, vxb, vyb, vzb);
    return 3;	
  }
  
    
  return 0; //outside the container
}


void PowderTarget3DSolidBoundary::reflectRight(double x, double y,double z, double pressure, double vx, double vy,  double vz, 
     vector<double>& xb, vector<double>& yb, vector<double>&zb, vector<double>& pressureb, 
     vector<double>& vxb, vector<double>& vyb, vector<double>& vzb) {
	
        xb.push_back(2*R-x);
	yb.push_back(y); 
	zb.push_back(z);
	pressureb.push_back(pressure);
	if(ip2d(vx,vy,1,0)<=0) { // leaving solid boundary
		vxb.push_back(0);
		vyb.push_back(0);
		vzb.push_back(0);
	}
	else {
		vxb.push_back(-vx);
		vyb.push_back(vy);
		vzb.push_back(vz);
	}

}

void PowderTarget3DSolidBoundary::reflectLeft(double x, double y, double z, double pressure, double vx, double vy,  double vz, 
     vector<double>& xb, vector<double>& yb, vector<double>&zb, vector<double>& pressureb, 
     vector<double>& vxb, vector<double>& vyb, vector<double>& vzb) {
	
        xb.push_back(-2*R-x);
	yb.push_back(y);
	zb.push_back(z);
	pressureb.push_back(pressure);
	if(ip2d(vx,vy,-1,0)<=0) { // leaving solid boundary
		vxb.push_back(0);
		vyb.push_back(0);
		vzb.push_back(0);	
	}
	else {
		vxb.push_back(-vx);
		vyb.push_back(vy);
		vzb.push_back(vz);
	}
}

void PowderTarget3DSolidBoundary::reflectFront(double x, double y, double z, double pressure, double vx, double vy,  double vz, 
     vector<double>& xb, vector<double>& yb, vector<double>&zb, vector<double>& pressureb, 
     vector<double>& vxb, vector<double>& vyb, vector<double>& vzb) {
	
        xb.push_back(x);
	yb.push_back(y);
	zb.push_back(length - z);
	pressureb.push_back(pressure);
	if(ip3d(vx,vy,vz,0,0,1)<=0) { // leaving solid boundary
		vxb.push_back(0);
		vyb.push_back(0);
		vzb.push_back(0);	
	}
	else {
		vxb.push_back(vx);
		vyb.push_back(vy);
		vzb.push_back(-vz);
	}
}

void PowderTarget3DSolidBoundary::reflectBack(double x, double y, double z, double pressure, double vx, double vy,  double vz, 
     vector<double>& xb, vector<double>& yb, vector<double>&zb, vector<double>& pressureb, 
     vector<double>& vxb, vector<double>& vyb, vector<double>& vzb) {
        xb.push_back(x);
	yb.push_back(y);
	zb.push_back(-length - z);
	pressureb.push_back(pressure);
	if(ip3d(vx,vy,vz,0,0,-1)<=0) { // leaving solid boundary
		vxb.push_back(0);
		vyb.push_back(0);
		vzb.push_back(0);	
	}
	else {
		vxb.push_back(vx);
		vyb.push_back(vy);
		vzb.push_back(-vz);
	}
}


void PowderTarget3DSolidBoundary::reflectSemiCircle(double x, double y, double z, double pressure, double vx, double vy,  double vz, 
     vector<double>& xb, vector<double>& yb, vector<double>&zb, vector<double>& pressureb, 
     vector<double>& vxb, vector<double>& vyb, vector<double>& vzb){
    double dist = sqrt(x*x+y*y);
    double factor = (2.*R-dist)/dist;
    double normal_vx = x/dist;
    double normal_vy = y/dist;
    
    xb.push_back(factor*x);
    yb.push_back(factor*y);
    zb.push_back(z);
    pressureb.push_back(pressure);
    double ip = ip2d(vx,vy,normal_vx,normal_vy);
    if(ip<=0) { // leaving solid boundary
      vxb.push_back(0);
      vyb.push_back(0);	
      vzb.push_back(0);
    }
    else { 
      vxb.push_back(vx-2.*ip*normal_vx);
      vyb.push_back(vy-2.*ip*normal_vy);
      vzb.push_back(vz);
        }
}


void PowderTarget3DSolidBoundary::reflectLeftFront(double x, double y, double z, double pressure, double vx, double vy,  double vz, 
     vector<double>& xb, vector<double>& yb, vector<double>&zb, vector<double>& pressureb, 
     vector<double>& vxb, vector<double>& vyb, vector<double>& vzb){
      	xb.push_back(-2*R-x);
	yb.push_back(y);
	zb.push_back(length - z); // 2 *length/2.0 -z
	pressureb.push_back(pressure);
	double normal_x = -1/std::sqrt(2.), normal_z = 1/std::sqrt(2.);
	double ip = ip2d(vx, vz, normal_x, normal_z);

	if(ip<=0) { // leaving solid boundary
		vxb.push_back(0);
		vyb.push_back(0);
		vzb.push_back(0);	
	}
        else { 
		vxb.push_back(vx-2.*ip*normal_x);
		vyb.push_back(vy);
		vzb.push_back(vz-2 *ip*normal_z);     
	}
}

void PowderTarget3DSolidBoundary::reflectRightFront(double x, double y, double z, double pressure, double vx, double vy,  double vz, 
     vector<double>& xb, vector<double>& yb, vector<double>&zb, vector<double>& pressureb, 
     vector<double>& vxb, vector<double>& vyb, vector<double>& vzb){
      	xb.push_back(2*R-x);
	yb.push_back(y);
	zb.push_back(length - z); // 2 *length/2.0 -z
	pressureb.push_back(pressure);
	double normal_x = 1/std::sqrt(2.), normal_z = 1/std::sqrt(2.);
	double ip = ip2d(vx, vz, normal_x, normal_z);

	if(ip<=0) { // leaving solid boundary
		vxb.push_back(0);
		vyb.push_back(0);
		vzb.push_back(0);	
	}
        else { 
		vxb.push_back(vx-2.*ip*normal_x);
		vyb.push_back(vy);
		vzb.push_back(vz-2 *ip*normal_z);     
	}
}
void PowderTarget3DSolidBoundary::reflectSemiCircleFront(double x, double y, double z, double pressure, double vx, double vy,  double vz, 
     vector<double>& xb, vector<double>& yb, vector<double>&zb, vector<double>& pressureb, 
     vector<double>& vxb, vector<double>& vyb, vector<double>& vzb){
    double dist = sqrt(x*x+y*y);
    double factor = (2.*R-dist)/dist;
    double normal_x = x/dist/std::sqrt(2.);
    double normal_y = y/dist/std::sqrt(2.);
    double normal_z = 1/std::sqrt(2.);
    
    xb.push_back(factor*x);
    yb.push_back(factor*y);
    zb.push_back(length - z);
    pressureb.push_back(pressure);
    double ip = ip3d(vx, vy, vz, normal_x, normal_y, normal_z); 

    if(ip<=0) { // leaving solid boundary
      vxb.push_back(0);
      vyb.push_back(0);	
      vzb.push_back(0);
    }
    else { 
      vxb.push_back(vx-2.*ip*normal_x);
      vyb.push_back(vy-2.*ip*normal_y);
      vzb.push_back(vz-2.*ip*normal_z);
        }

}

void PowderTarget3DSolidBoundary::reflectLeftBack(double x, double y, double z, double pressure, double vx, double vy,  double vz, 
     vector<double>& xb, vector<double>& yb, vector<double>&zb, vector<double>& pressureb, 
     vector<double>& vxb, vector<double>& vyb, vector<double>& vzb){
      	xb.push_back(-2*R-x);
	yb.push_back(y);
	zb.push_back(-length - z); // -2 *length/2.0 -z
	pressureb.push_back(pressure);
	double normal_x = -1/std::sqrt(2.), normal_z = -1/std::sqrt(2.);
	double ip = ip2d(vx, vz, normal_x, normal_z);
      
	if(ip<=0) { // leaving solid boundary
		vxb.push_back(0);
		vyb.push_back(0);
		vzb.push_back(0);	
	}
        else { 
		vxb.push_back(vx-2.*ip*normal_x);
		vyb.push_back(vy);
		vzb.push_back(vz-2 *ip*normal_z);     
	}
}


void PowderTarget3DSolidBoundary::reflectRightBack(double x, double y, double z, double pressure, double vx, double vy,  double vz, 
     vector<double>& xb, vector<double>& yb, vector<double>&zb, vector<double>& pressureb, 
     vector<double>& vxb, vector<double>& vyb, vector<double>& vzb){
      	xb.push_back(2*R-x);
	yb.push_back(y);
	zb.push_back(-length - z); // -2 *length/2.0 -z
	pressureb.push_back(pressure);
	double normal_x = 1/std::sqrt(2.), normal_z = -1/std::sqrt(2.);
	double ip = ip2d(vx, vz, normal_x, normal_z);

	if(ip<=0) { // leaving solid boundary
		vxb.push_back(0);
		vyb.push_back(0);
		vzb.push_back(0);	
	}
        else { 
		vxb.push_back(vx-2.*ip*normal_x);
		vyb.push_back(vy);
		vzb.push_back(vz-2 *ip*normal_z);     
	}
}
void PowderTarget3DSolidBoundary::reflectSemiCircleBack(double x, double y, double z, double pressure, double vx, double vy,  double vz, 
     vector<double>& xb, vector<double>& yb, vector<double>&zb, vector<double>& pressureb, 
     vector<double>& vxb, vector<double>& vyb, vector<double>& vzb){
    double dist = sqrt(x*x+y*y);
    double factor = (2.*R-dist)/dist;
    double normal_x = x/dist/std::sqrt(2.);
    double normal_y = y/dist/std::sqrt(2.);
    double normal_z = -1/std::sqrt(2.);
    
    xb.push_back(factor*x);
    yb.push_back(factor*y);
    zb.push_back(-length - z);
    pressureb.push_back(pressure);
    double ip = ip3d(vx, vy, vz, normal_x, normal_y, normal_z); 

    if(ip<=0) { // leaving solid boundary
      vxb.push_back(0);
      vyb.push_back(0);	
      vzb.push_back(0);
    }
    else { 
      vxb.push_back(vx-2.*ip*normal_x);
      vyb.push_back(vy-2.*ip*normal_y);
      vzb.push_back(vz-2.*ip*normal_z);
        }
}


////////////////////////////////////////////////////////////////////////////////////////
// End of PowderTarget3DSolidBoundary
////////////////////////////////////////////////////////////////////////////////////////
